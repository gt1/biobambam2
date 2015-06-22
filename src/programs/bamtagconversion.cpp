/**
    biobambam
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
    
    
    Program to read the old, incorrect tags and ouput correct ones.
    
    MS:i -> ms:i
    MC:i -> mc:i (MC:Z already in spec)
    MT:i -> mt:i
    
    Andrew Whitwham, May 2015
**/

#include "config.h"

#include <iostream>

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

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>

static int getDefaultMD5() { return 0; }

int bamtagconversion(libmaus2::util::ArgInfo const &arginfo) {
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
    
    std::string const headertext(header.text);

    // add PG line to header
    std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
	    headertext,
	    "bamtagconversion", // ID
	    "bamtagconversion", // PN
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
    
    if (verbose > 0) {
    	std::cerr << "Running bamtagconversion" << std::endl;
    }
    
    libmaus2::bambam::BamAuxFilterVector MSfilter;
    libmaus2::bambam::BamAuxFilterVector MCfilter;
    libmaus2::bambam::BamAuxFilterVector MTfilter;
    MSfilter.set("MS"); MSfilter.set("ms"); 
    MCfilter.set("MC"); MCfilter.set("mc");
    MTfilter.set("MT"); MTfilter.set("mt");
   

    while (dec.readAlignment()) {
	bool found = false;
	
	if (algn.hasAuxNC("MS")) {
	    int64_t const score = algn.getAuxAsNumberNC<int64_t>("ms");
	    algn.filterOutAux(MSfilter);
	    algn.putAuxNumber("ms",'i',score);

	    if (verbose > 0) {
	    	std::cerr << "found ms " << score << " ";
	    	found = true;
	    }
	}
	
	if (algn.hasAuxNC("MC")) {
	    int64_t mate_coord;
	    
	    if (algn.getAuxAsNumberNC<int64_t>("mc", mate_coord)) {
		if (verbose > 0) {
	    	    std::cerr << "found mc " << mate_coord << " ";
	    	    found = true;
		}
	    
	    	algn.filterOutAux(MCfilter);
		algn.putAuxNumber("mc",'i',mate_coord);
	    }
	}
	
	if (algn.hasAuxNC("MT")) {
	    const char *mt = algn.getAuxStringNC("MT");
	    char *copy = strdup(mt);
	    
	    algn.filterOutAux(MTfilter);
	    algn.putAuxString("mt", copy);
	
	    if (verbose > 0) {
	    	std::cerr << "found mt " << copy << " ";
	    	found = true;
	    }
	    
	    free(copy);
	}
	
	if (found) {
	    std::cerr << std::endl;
	}
	
    	writer->writeAlignment(algn);
    }
    
    writer.reset();
    
    if (Pmd5cb) {
    	Pmd5cb->saveDigestAsFile(md5filename);
    }

    return EXIT_SUCCESS;
}


int main(int argc, char *argv[]) {
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
		V.push_back(std::pair<std::string,std::string> ("outputthreads=<1>", "output helper threads (for outputformat=bam only, default: 1)"));
		
    	    	::biobambam2::Licensing::printMap(std::cerr,V);
		
		std::cerr << std::endl;
		return EXIT_SUCCESS;
	    }
	}
	
	return bamtagconversion(arginfo);
	
    } catch (std::exception const &ex) {
    	std::cerr << ex.what() << std::endl;
	return EXIT_FAILURE;
    }

}

    	
