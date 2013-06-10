/**
    bambam
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

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
#include <libmaus/bambam/BamAlignmentNameComparator.hpp>
#include <libmaus/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamEntryContainer.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus/lz/SnappyCompress.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GetObject.hpp>
#include <libmaus/util/PutObject.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }

int bamfilteraux(::libmaus::util::ArgInfo const & arginfo)
{
	::libmaus::util::TempFileRemovalContainer::setup();

	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	if ( arginfo.hasArg("keep") && arginfo.hasArg("remove") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "The keep and remove keys are mutually exclusive." << std::endl;
		se.finish();
		throw se;		
	}

	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	
	switch ( level )
	{
		case Z_NO_COMPRESSION:
		case Z_BEST_SPEED:
		case Z_BEST_COMPRESSION:
		case Z_DEFAULT_COMPRESSION:
			break;
		default:
		{
			::libmaus::exception::LibMausException se;
			se.getStream()
				<< "Unknown compression level, please use"
				<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
				<< " level=" << Z_BEST_SPEED << " (fast) or"
				<< " level=" << Z_BEST_COMPRESSION << " (best) or"
				<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
			se.finish();
			throw se;
		}
			break;
	}

	::libmaus::bambam::BamDecoder dec(std::cin,false);
	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamfilteraux", // ID
		"bamfilteraux", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	libmaus::bambam::BamHeader const uphead(upheadtext);
 	libmaus::bambam::BamWriter writer(std::cout,uphead,level);
 	libmaus::bambam::BamAuxFilterVector bafv;
 	
 	std::string const filterlist = arginfo.hasArg("keep") ? arginfo.getUnparsedValue("keep","") : arginfo.getUnparsedValue("remove","");
	std::deque<std::string> const tokens = libmaus::util::stringFunctions::tokenize<std::string>(filterlist,std::string(","));
	
	for ( uint64_t i = 0; i < tokens.size(); ++i )
		if ( tokens[i].size() == 2 )
		{
			bafv.set(tokens[i]);
		}
		else
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Invalid tag name " << tokens[i] << std::endl;
			se.finish();
			throw se; 			
		}

	libmaus::bambam::BamAlignment & algn = dec.getAlignment();
	uint64_t c = 0;

 	if ( 
 		arginfo.hasArg("remove") 
 		|| 
 		(
 			(!arginfo.hasArg("remove"))
 			&&
 			(!arginfo.hasArg("keep"))
		)
	)
 	{	
 		while ( dec.readAlignment() )
 		{
 			algn.filterOutAux(bafv);
 			algn.serialise(writer.getStream());
 			
 			if ( verbose && (++c & (1024*1024-1)) == 0 )
 				std::cerr << "[V] " << c/(1024*1024) << std::endl;
 		}
 	}
 	else
 	{
 		while ( dec.readAlignment() )
 		{
 			algn.filterAux(bafv);
 			algn.serialise(writer.getStream());

 			if ( verbose && (++c & (1024*1024-1)) == 0 )
 				std::cerr << "[V] " << c/(1024*1024) << std::endl;
 		} 		
 	}

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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );
				V.push_back ( std::pair<std::string,std::string> ( "keep=<[]>", "keep these tags (default: keep all)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "remove=<[]>","remove these tags (default: remove none)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				
				std::cerr << "The keep and remove keys are mutually exclusive. Tags are given by their two character ids. Multiple ids are separated by commas." << std::endl;
				
				return EXIT_SUCCESS;
			}
			
		return bamfilteraux(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

