/**
    bambam
    Copyright (C) 2015 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.    
    
    Andrew Whitwham, April 2015
**/
#include "config.h"

#include <iostream>
#include <queue>
#include <vector>

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamHeaderUpdate.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>


static int getDefaultVerbose() {return 0;}
static int getDefaultLevel() {return Z_DEFAULT_COMPRESSION;}
static int getDefaultMaxMisordered() {return 500000;}
static int getDefaultExpectedStep() {return 2;}

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>

static int getDefaultMD5() { return 0; }


// structs to hold unsorted data
struct Grouped {
    uint64_t rank;
    libmaus2::bambam::BamAlignment *algn;
    
    Grouped(uint64_t irank, libmaus2::bambam::BamAlignment *ialgn) : rank(irank), algn(ialgn) {
    }
    
    ~Grouped() {
    	delete algn;
    }
};


struct Unsorted {
    uint64_t rank;
    std::vector<Grouped *> group;
    
    Unsorted(uint64_t irank, std::vector<Grouped *> igroup) : rank(irank), group(igroup) {
    }
    
    ~Unsorted() {
    }
};


static bool rank_cmp(Unsorted *one, Unsorted *two) {
    return one->rank < two->rank;
}


static uint64_t get_rank(libmaus2::bambam::BamAlignment &algn) {
    char const *name = algn.getName();
    char const *u1 = name;
    bool ok = true;
    uint64_t rank = 0;

    while (*u1 && *u1 != '_') {
	rank *= 10;
	rank += (*u1 - '0');
	ok = ok && isdigit(*u1);
	++u1;
    }

    // unable to find rank?
    if (!ok) {
	::libmaus2::exception::LibMausException se;
	se.getStream() << "Aborting, no rank found.  Read name is " << name << std::endl;
	se.finish();
	throw se;
    }

    return rank;
}


uint64_t get_grouped_reads(libmaus2::bambam::BamAlignment &alg, libmaus2::bambam::BamAlignmentDecoder &decoder, std::vector<Grouped *> &group, bool &cont) {
    uint64_t const rank = get_rank(alg);
    uint64_t found_rank = rank;
    bool first_found  = false;
    bool second_found = false;
    
    while (found_rank == rank) {
	// store read line
	Grouped *tmp_group = new Grouped(rank, new libmaus2::bambam::BamAlignment); 
	tmp_group->algn->swap(alg);
	group.push_back(tmp_group);
	
	if (tmp_group->algn->isRead1()) first_found  = true;
	if (tmp_group->algn->isRead2()) second_found = true;
	
	if (decoder.readAlignment()) {
	    found_rank = get_rank(alg);
	} else {
	    // unable to read more reads
	    cont = false;
	    break;
	}
    }
    
    if (!(first_found && second_found)) {
	::libmaus2::exception::LibMausException se;
	se.getStream() << "Aborting, rank " << rank << "_" << (rank + 1) << " lacks first and second reads " << first_found << " " << second_found << std::endl;
	se.finish();
	throw se;
    }
    
    return rank;
}

    
void write_grouped_reads(std::vector<Grouped *> &group, libmaus2::bambam::BamBlockWriterBase::unique_ptr_type const &scribe, int const verbose) {

    for (std::vector<Grouped *>::iterator itr = group.begin(); itr != group.end(); ++itr) {
	scribe->writeAlignment(*((*itr)->algn));
	uint64_t rank = (*itr)->rank;
	
	if (verbose > 1) {
	    std::cerr << "Write rank " << rank << std::endl;
	}
	
	delete (*itr);
    }

    group.clear();
}


void store_group(std::vector<Grouped *> &group, std::vector<Unsorted *> &unsorted) {
	Unsorted *tmp_unsorted = new Unsorted(group[0]->rank, group); 
	unsorted.push_back(tmp_unsorted);
}


bool check_ranks(std::vector<Unsorted *> const &unsorted, uint64_t &wanted, int const verbose, int const step) {
    wanted = unsorted[0]->rank;
    
    if (verbose > 0) {
    	std::cerr << "check_ranks called on " << unsorted.size() << " reading groups" << std::endl;
    }

    for (std::vector<Unsorted *>::const_iterator itr = unsorted.begin(); itr != unsorted.end(); ++itr) {
    	if ((*itr)->rank == wanted) {
	    wanted += step;
	} else {
	    if (verbose > 0) {
	    	std::cerr << "Not found rank " << wanted << " continuing search." << std::endl;
	    } 

	    return false;
	}
    }
    
    return true;
}

 	
void search_for_missing_group(libmaus2::bambam::BamAlignment &alg, libmaus2::bambam::BamAlignmentDecoder &decoder,
    	    	    	       std::vector<Grouped *> &group, libmaus2::bambam::BamBlockWriterBase::unique_ptr_type 
			       const &scribe, bool &cont, uint64_t &expected, int const verbose,
			       int const misordered, int const step) {
    std::vector<Unsorted *> unsorted;
    
    // store current group
    store_group(group, unsorted);
    group.clear();
    
    bool still_missing  = true;
    bool found          = false;
    int const interval  = 10;
    int count           = 0;
    int max_count       = 0;
    
    while (still_missing) {
    	uint64_t const current_rank = get_grouped_reads(alg, decoder, group, cont);
    	store_group(group, unsorted);
    	group.clear();
	
	count++;
	
	if (current_rank == expected) found = true;
	
	if ((count >= interval || !cont) && found) { 
	
	    sort(unsorted.begin(), unsorted.end(), rank_cmp);
	    
	    if (check_ranks(unsorted, expected, verbose, step)) {
	    	still_missing = false;
		
		for (std::vector<Unsorted *>::iterator itr = unsorted.begin(); itr != unsorted.end(); ++itr) {
		    write_grouped_reads((*itr)->group, scribe, verbose);
		    delete (*itr);
		}
	    } else {
	    	max_count += count;
	    	count = 0;
		found = false;
	    }
	}
	
	if (!cont && still_missing) {
	    ::libmaus2::exception::LibMausException se;
	    se.getStream() << "Aborting, rank " << expected << "_" << (expected + 1) << " not found, end of file" << std::endl;
	    se.finish();
	    throw se;
	} else if (((count + max_count) >= misordered) && still_missing) {
	    ::libmaus2::exception::LibMausException se;
	    se.getStream() << "Aborting, number of misordered is greater than " << misordered << std::endl;
	    se.finish();
	    throw se;
	}
    }
}    	
    

int bamranksort(libmaus2::util::ArgInfo const &arginfo) {

    if (isatty(STDIN_FILENO)) {
	::libmaus2::exception::LibMausException se;
	se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
	se.finish();
	throw se;
    }

    if (isatty(STDOUT_FILENO)) {
	::libmaus2::exception::LibMausException se;
	se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
	se.finish();
	throw se;
    }

    int const verbose = arginfo.getValue<int>("verbose", getDefaultVerbose());
    int const misordered = arginfo.getValue<int>("misordered", getDefaultMaxMisordered());
    int const step = arginfo.getValue<int>("step", getDefaultExpectedStep());

    // input decoder wrapper
    libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
	    libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
		    arginfo,false // put rank
	    )
    );

    ::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
    ::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;

    libmaus2::bambam::BamAlignment &algn = dec.getAlignment();
    libmaus2::bambam::BamHeader const &header = dec.getHeader();
    
    if (verbose > 0) {
    	std::cerr << "Running bamranksort" << std::endl;
    }

    std::string const headertext(header.text);

    // add PG line to header
    std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
	    headertext,
	    "bamranksort", // ID
	    "bamranksort", // PN
	    arginfo.commandline, // CL
	    ::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
	    std::string(PACKAGE_VERSION) // VN			
    );
    
    ::libmaus2::bambam::BamHeader uphead(upheadtext);
    
    /*
     * start md5 callbacks
     */
    std::string md5filename;

    std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > cbs;
    ::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;

    if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
    {
	    if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
		    md5filename = arginfo.getUnparsedValue("md5filename","");
	    else
		    std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

	    if ( md5filename.size() )
	    {
		    ::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus2::lz::BgzfDeflateOutputCallbackMD5);
		    Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
		    cbs.push_back(Pmd5cb.get());
	    }
    }

    std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;

    if ( cbs.size() )
	    Pcbs = &cbs;
    /*
     * end md5
     */

    libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(
	    libmaus2::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, Pcbs)
    );

    // major assumption is that all alignments are grouped by read name
    // containing forward and reverse reads plus any number of secondary and
    // supplemental reads.  Isolated reads will be a cause of program termination.  
    
    // get the first read
    if (!dec.readAlignment()) {
	::libmaus2::exception::LibMausException se;
	se.getStream() << "Aborting, unable to read first alignment" << std::endl;
	se.finish();
	throw se;
    }
    
    std::vector<Grouped *> grouped_reads;
    
    uint64_t expected_rank = 0;
    bool keep_going = true;
    
    while (keep_going) {
    	uint64_t const current_rank = get_grouped_reads(algn, dec, grouped_reads, keep_going);
	
	if (current_rank == expected_rank) {
	    // we have what we came for
	    write_grouped_reads(grouped_reads, writer, verbose);
	    expected_rank += step;
	} else {
	    if (verbose > 0) {
	    	std::cerr << "Found rank " << current_rank << " expecting " << expected_rank << ".  Starting search." << std::endl;
	    } 
	
	    search_for_missing_group(algn, dec, grouped_reads, writer, keep_going, expected_rank, verbose, misordered, step);
	}
    }

    writer.reset();
    
    if (Pmd5cb) {
    	Pmd5cb->saveDigestAsFile(md5filename);
    }

    return EXIT_SUCCESS;
}


int main(int argc, char* argv[]) {
    try {
    	::libmaus2::util::ArgInfo const arginfo(argc, argv);
	
	for (uint64_t i = 0; i < arginfo.restargs.size(); ++i) {
	    if (arginfo.restargs[i] == "-v" ||
	    	arginfo.restargs[i] == "--version") {
	    	std::cerr << ::biobambam2::Licensing::license();
		return EXIT_SUCCESS;
		
	    } else if (arginfo.restargs[i] == "-h" ||
	    	    	arginfo.restargs[i] == "--help") {
		std::cerr << ::biobambam2::Licensing::license();
		std::cerr << std::endl;
		std::cerr << "Key=Value pairs:" << std::endl;
		std::cerr << std::endl;
		
		std::vector<std::pair<std::string,std::string>> V;
			
		V.push_back(std::pair<std::string,std::string> ("level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText()));
		V.push_back(std::pair<std::string,std::string> ("verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report"));
		V.push_back(std::pair<std::string,std::string> ("md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)"));
		V.push_back(std::pair<std::string,std::string> ("md5filename=<filename>", "file name for md5 check sum (default: extend output file name)"));
		V.push_back(std::pair<std::string,std::string> ("misordered=<["+::biobambam2::Licensing::formatNumber(getDefaultMaxMisordered())+"]>", "number of read pairs (including secondary/supplementary alignments) allowed to accumulate before exiting"));
		V.push_back(std::pair<std::string,std::string> ("step=<["+::biobambam2::Licensing::formatNumber(getDefaultExpectedStep())+"]>", "the increment expected between one rank and the next"));
		V.push_back(std::pair<std::string,std::string> ("outputthreads=<1>", "output helper threads (for outputformat=bam only, default: 1)"));

    	    	::biobambam2::Licensing::printMap(std::cerr,V);
		
		std::cerr << std::endl;
		return EXIT_SUCCESS;
	    }
	}
	
	return bamranksort(arginfo);
	
    } catch (std::exception const &ex) {
    	std::cerr << ex.what() << std::endl;
	return EXIT_FAILURE;
    }

}

