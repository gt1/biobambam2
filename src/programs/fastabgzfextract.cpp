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
#include <config.h>
#include <cstdlib>
#include <iostream>
#include <libmaus2/aio/PosixFdInputStream.hpp>
#include <libmaus2/fastx/FastABgzfIndex.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/GetFileSize.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

void fastabgzfextract(libmaus2::util::ArgInfo const & arginfo)
{
	if ( ! arginfo.hasArg("reference") )
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "reference key is missing." << std::endl;
		se.finish();
		throw se;			
	}
	
	std::string const reference = arginfo.getUnparsedValue("reference","");
	
	if ( ! libmaus2::util::GetFileSize::fileExists(reference) )
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "file " << reference << " does not exist." << std::endl;
		se.finish();
		throw se;				
	}

	if ( ! libmaus2::util::GetFileSize::fileExists(reference+".idx") )
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "file " << reference << " does not exist." << std::endl;
		se.finish();
		throw se;				
	}
	
	libmaus2::aio::PosixFdInputStream PFIS(reference,128*1024);
	libmaus2::aio::InputStreamInstance indexCIS(reference+".idx");
	libmaus2::fastx::FastABgzfIndex index(indexCIS);
	libmaus2::autoarray::AutoArray<char> B;
	
	while ( std::cin )
	{
		std::string line;
		std::getline(std::cin,line);
		
		if ( line.size() )
		{
			std::deque<std::string> tokens = ::libmaus2::util::stringFunctions::tokenize(line,std::string("\t"));
			
			if ( tokens.size() != 3 )
				continue;
			
			std::string const & name = tokens[0];
			std::istringstream posistr(tokens[1]);
			std::istringstream lenistr(tokens[2]);
			uint64_t pos, len;
			
			posistr >> pos;
			lenistr >> len;
			
			int64_t const thisseqid = index.getSequenceId(name);
			
			if ( thisseqid >= 0 )
			{
				libmaus2::fastx::FastABgzfDecoder::unique_ptr_type Pstr = index.getStream(PFIS,thisseqid);
				Pstr->seekg(pos);
				
				if ( len > B.size() )
					B = libmaus2::autoarray::AutoArray<char>(len,false);
					
				Pstr->read(B.begin(),len);
				uint64_t const rlen = Pstr->gcount();
				
				std::cout.write(B.begin(),rlen);
				std::cout.put('\n');
			}
		}
	}
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo const arginfo(argc,argv);

		
		for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
			if ( 
				arginfo.restargs[i] == "-v"
				||
				arginfo.restargs[i] == "--version"
			)
			{
				std::cerr << ::biobambam2::Licensing::license();
				return EXIT_SUCCESS;
			}
			else if ( 
				arginfo.restargs[i] == "-h"
				||
				arginfo.restargs[i] == "--help"
			)
			{
				std::cerr << ::biobambam2::Licensing::license();
				std::cerr << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
			
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA.bgzf" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}


		fastabgzfextract(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
