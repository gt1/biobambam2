/**
    bambam
    Copyright (C) 2009-2014 German Tischler
    Copyright (C) 2011-2014 Genome Research Limited

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
**/
#include "config.h"

#include <iostream>
#include <queue>

#include <libmaus/aio/CheckedOutputStream.hpp>

#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/bambam/MdNmRecalculation.hpp>

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/lz/BufferedGzipStream.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 1; }
static bool getDefaultDisableValidation() { return false; }
static std::string getDefaultInputFormat() { return "bam"; }

struct IntervalPairHistogram
{
	std::vector<std::string> names;
	std::map<std::string,uint64_t> nameids;
	std::map<uint64_t,uint64_t> hist;
	std::map< uint64_t,std::vector<std::string> > readnames;
	
	IntervalPairHistogram()
	{
	
	}
	
	uint64_t getNameId(std::string const & name)
	{
		if ( nameids.find(name) == nameids.end() )
		{
			uint64_t const id = nameids.size();
			nameids[name] = id;
			names.push_back(name);
		}
		
		return nameids.find(name)->second;
	}
		
	static std::vector< std::string > getComments(libmaus::bambam::BamAlignment * algn)
	{
		std::vector< std::string > V;
		
		if ( algn->hasAux("CO") )
		{
			char const * co = algn->getAuxString("CO");
			assert ( co );
			
			std::deque<std::string> tokens = libmaus::util::stringFunctions::tokenize<std::string>(std::string(co),std::string(";"));
			
			V = std::vector< std::string >(tokens.begin(),tokens.end());
		}
		
		return V;
	}
	
	std::ostream & printHistogram(std::ostream & out) const
	{
		std::map< uint64_t,std::vector<uint64_t> > histhist;
		for ( std::map<uint64_t,uint64_t>::const_iterator ita = hist.begin(); ita != hist.end(); ++ita )
		{
			uint64_t const ids = ita->first;
			#if 0
			uint64_t const id1 = (ita->first >> 32) & 0xFFFFFFFFULL;
			uint64_t const id2 = (ita->first >>  0) & 0xFFFFFFFFULL;
			#endif
			uint64_t const cnt = ita->second;
			
			histhist[cnt].push_back(ids);
		}
		
		for ( std::map< uint64_t,std::vector<uint64_t> >::const_reverse_iterator ita = histhist.rbegin();
			ita != histhist.rend(); ++ita )
		{
			uint64_t const cnt = ita->first;
			std::vector<uint64_t> const & V = ita->second;
			
			for ( uint64_t i = 0; i < V.size(); ++i )
			{
				uint64_t const ids = V[i];
				uint64_t const id1 = (ids >> 32) & 0xFFFFFFFFULL;
				uint64_t const id2 = (ids >>  0) & 0xFFFFFFFFULL;
				std::string const & n1 = names[id1];
				std::string const & n2 = names[id2];
				
				if ( n1 < n2 )
					out << cnt << "\t" << n1 << "\t" << n2;
				else
					out << cnt << "\t" << n2 << "\t" << n1;
				
				std::map< uint64_t,std::vector<std::string> >::const_iterator it = readnames.find(ids);
				assert ( it != readnames.end() );
				std::vector<std::string> const & rnl = it->second;
				
				for ( uint64_t j = 0; j < rnl.size(); ++j )
				{
					out << "\t" << rnl[j];
				}
				
				out << "\n";
			}
		}
		
		return out;
	}
	
	void evaluateList(std::vector<libmaus::bambam::BamAlignment *> & curlist)
	{
		std::vector<libmaus::bambam::BamAlignment *> R1;
		std::vector<libmaus::bambam::BamAlignment *> R2;
		
		for ( uint64_t i = 0; i < curlist.size(); ++i )
		{
			if ( curlist[i]->isRead1() )
				R1.push_back(curlist[i]);
			else if ( curlist[i]->isRead2() )
				R2.push_back(curlist[i]);
		}
		
		for ( uint64_t i = 0; i < R1.size(); ++i )
			for ( uint64_t j = 0; j < R2.size(); ++j )
			{
				libmaus::bambam::BamAlignment * A1 = R1[i];
				libmaus::bambam::BamAlignment * A2 = R2[j];
				
				char const * rname1 = A1->getName();
				char const * rname2 = A2->getName();
				
				assert ( strcmp(rname1,rname2) == 0 );
				
				std::vector< std::string > co1 = getComments(A1);
				std::vector< std::string > co2 = getComments(A2);
				
				for ( uint64_t k = 0; k < co1.size(); ++k )
					for ( uint64_t l = 0; l < co2.size(); ++l )
					{
						std::string const & n1 = co1[k];
						std::string const & n2 = co2[l];
						
						uint64_t const id1 = getNameId(n1);
						uint64_t const id2 = getNameId(n2);
						
						uint64_t const id = (id1 <= id2) ? ((id1<<32) | id2) : ((id2<<32) | id1);

						// std::cerr << n1 << "\t" << id1 << "\t" << n2 << "\t" << id2 << std::endl;
		
						hist[id]++;
						readnames[id].push_back(rname1);
					}
			}
		
		// std::cerr << R1.size() << "\t" << R2.size() << std::endl;
	}
};

int bamintervalcommenthist(::libmaus::util::ArgInfo const & arginfo)
{
	::libmaus::util::TempFileRemovalContainer::setup();
	
	bool const inputisstdin = (!arginfo.hasArg("I")) || (arginfo.getUnparsedValue("I","-") == "-");

	if ( isatty(STDIN_FILENO) && inputisstdin && (arginfo.getValue<std::string>("inputformat","bam") != "sam") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const disablevalidation = arginfo.getValue<int>("disablevalidation",getDefaultDisableValidation());

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	if ( disablevalidation )
		dec.disableValidation();
	::libmaus::bambam::BamHeader const & header = dec.getHeader();
	
	libmaus::bambam::BamAlignment & curalgn = dec.getAlignment();
	uint64_t c = 0;
	
	std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> algnpool;
	std::vector<libmaus::bambam::BamAlignment *> algnfreelist;
	std::vector<libmaus::bambam::BamAlignment *> curlist;
	IntervalPairHistogram hist;
	
	while ( dec.readAlignment() )
	{
		if ( 
			(curlist.size()) 
			&&
			(
				std::string(curalgn.getName()) !=
				std::string(curlist.back()->getName())
			)
		)
		{
			hist.evaluateList(curlist);
				
			// move pointers to free list
			while ( curlist.size() )
			{
				algnfreelist.push_back(curlist.back());
				curlist.pop_back();
			}
		}
		
		if ( ! algnfreelist.size() )
		{
			libmaus::bambam::BamAlignment::shared_ptr_type ptr(new libmaus::bambam::BamAlignment);
			algnpool.push_back(ptr);
			algnfreelist.push_back(algnpool.back().get());
		}
		
		assert ( algnfreelist.size() );
		libmaus::bambam::BamAlignment * addalgn = algnfreelist.back();
		algnfreelist.pop_back();
		
		curalgn.swap(*addalgn);
		curlist.push_back(addalgn);
	
		if ( (++c & (1024*1024-1)) == 0 && verbose )
			std::cerr << "[V] " << c << std::endl;		
	}
	
	if ( curlist.size() )
	{
		hist.evaluateList(curlist);
				
		// move pointers to free list
		while ( curlist.size() )
		{
			algnfreelist.push_back(curlist.back());
			curlist.pop_back();
		}	
	}
		
	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

	hist.printHistogram(std::cout);

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		
		for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
			if ( 
				arginfo.restargs[i] == "-v"
				||
				arginfo.restargs[i] == "--version"
			)
			{
				std::cerr << ::biobambam::Licensing::license();
				return EXIT_SUCCESS;
			}
			else if ( 
				arginfo.restargs[i] == "-h"
				||
				arginfo.restargs[i] == "--help"
			)
			{
				std::cerr << ::biobambam::Licensing::license();
				std::cerr << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );				

				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamintervalcommenthist(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
