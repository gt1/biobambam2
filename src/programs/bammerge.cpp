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

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/bambam/BamMergeCoordinate.hpp>
#include <libmaus/bambam/BamMergeQueryName.hpp>
#include <libmaus/bambam/BamWriter.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static std::string getDefaultSortOrder() { return "coordinate"; }

::libmaus::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus::util::ArgInfo const & arginfo,
	::libmaus::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bammerge", // ID
		"bammerge", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus::bambam::BamHeader(upheadtext));
	
	return UNIQUE_PTR_MOVE(uphead);
}

int bammerge(libmaus::util::ArgInfo const & arginfo)
{
	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	std::string const sortorder = arginfo.getValue<std::string>("SO",getDefaultSortOrder());
	
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

	std::vector<std::string> inputfilenames = arginfo.getPairValues("I");
	for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
		inputfilenames.push_back(arginfo.restargs[i]);

	if ( sortorder == "queryname" )
	{
		libmaus::bambam::BamMergeQueryName bamdec(inputfilenames /* ,true */);
		libmaus::bambam::BamAlignment const & algn = bamdec.getAlignment();
		libmaus::bambam::BamHeader const & header = bamdec.getHeader();
		::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));
		libmaus::bambam::BamWriter writer(std::cout,*uphead,level);
		libmaus::bambam::BamWriter::stream_type & bamoutstr = writer.getStream();
		while ( bamdec.readAlignment() )
			algn.serialise(bamoutstr);
	}
	else
	{
		libmaus::bambam::BamMergeCoordinate bamdec(inputfilenames /* ,true */);
		libmaus::bambam::BamAlignment const & algn = bamdec.getAlignment();
		libmaus::bambam::BamHeader const & header = bamdec.getHeader();
		::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));
		libmaus::bambam::BamWriter writer(std::cout,*uphead,level);
		libmaus::bambam::BamWriter::stream_type & bamoutstr = writer.getStream();
		while ( bamdec.readAlignment() )
			algn.serialise(bamoutstr);	
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
			
				V.push_back ( std::pair<std::string,std::string> ( "I=<[filename]>", "input file, can be set multiple times" ) );
				V.push_back ( std::pair<std::string,std::string> ( "SO=<["+getDefaultSortOrder()+"]>]", "sort order (coordinate or queryname)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bammerge(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

