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
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <bambam/Licensing.hpp>

#include <config.h>

static uint64_t getDefaultMaskNeg()
{
	return 0;
}

static int getDefaultLevel()
{
	return Z_DEFAULT_COMPRESSION;
}

int bammaskflags(::libmaus::util::ArgInfo const & arginfo)
{
	uint64_t const maskpos = arginfo.getValue<uint64_t>("maskpos",0xFFFFUL);
	uint64_t const maskneg = arginfo.getValue<uint64_t>("maskneg",getDefaultMaskNeg());
	uint64_t const mask = maskpos & (~maskneg);

	if ( mask )
	{
		std::cerr << "Keeping flags ";
		for ( uint64_t i = 1; i <= ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP; i <<= 1 )
			if ( mask & i )
				std::cerr << static_cast< ::libmaus::bambam::BamFlagBase::bam_flags >(i) << ";";
		std::cerr << std::endl;
		std::cerr << "Erasing flags ";
		for ( uint64_t i = 1; i <= ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP; i <<= 1 )
			if ( !(mask & i) )
				std::cerr << static_cast< ::libmaus::bambam::BamFlagBase::bam_flags >(i) << ";";
		std::cerr << std::endl;
	}
	else
	{
		std::cerr << "Erasing all flags." << std::endl;
	}

	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	
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

	::libmaus::bambam::BamDecoder BD(std::cin);
	::libmaus::bambam::BamHeader const & bamheader = BD.bamheader;

	std::string const headertext(bamheader.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bammaskflags", // ID
		"bammaskflags", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);

	::libmaus::bambam::BamAlignment & alignment = BD.alignment;
	::libmaus::bambam::BamWriter writer(std::cout,uphead,level);
	
	while ( BD.readAlignment() )
	{
		alignment.putFlags(alignment.getFlags() & mask);
		alignment.serialise(writer.bgzfos);
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
				std::cerr << ::bambam::Licensing::license();
				return EXIT_SUCCESS;
			}
			else if ( 
				arginfo.restargs[i] == "-h"
				||
				arginfo.restargs[i] == "--help"
			)
			{
				std::cerr << ::bambam::Licensing::license();
				std::cerr << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::bambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "maskneg=<["+::bambam::Licensing::formatNumber(getDefaultMaskNeg())+"]>", "flag masking bitmask" ) );

				::bambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bammaskflags(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
