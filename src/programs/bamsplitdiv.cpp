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
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamCat.hpp>
#include <libmaus/bambam/BamWriter.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static uint64_t getDefaultDiv() { return 1; }
static std::string getDefaultFilePrefix(::libmaus::util::ArgInfo const & arginfo) { return arginfo.getDefaultTmpFileName(); }

::libmaus::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus::util::ArgInfo const & arginfo,
	::libmaus::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamsplitmod", // ID
		"bamsplitmod", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus::bambam::BamHeader(upheadtext));
	
	return UNIQUE_PTR_MOVE(uphead);
}

int bamsplitmod(libmaus::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	uint64_t const div = arginfo.getValue<int>("div",getDefaultDiv());
	std::string const prefix = arginfo.getUnparsedValue("prefix",getDefaultFilePrefix(arginfo));
	
	if ( ! div )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "div cannot be 0." << std::endl;
		se.finish();
		throw se;
	}

	libmaus::bambam::BamDecoder bamdec(std::cin);
	libmaus::bambam::BamAlignment const & algn = bamdec.getAlignment();
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));
	
	libmaus::autoarray::AutoArray<libmaus::aio::CheckedOutputStream::unique_ptr_type> COS(div);
	libmaus::autoarray::AutoArray<libmaus::bambam::BamWriter::unique_ptr_type> writers(div);
	std::vector < std::string > filenames;
	for ( uint64_t i = 0; i < div; ++i )
	{
		std::ostringstream ostr;
		ostr << prefix << "_" << std::setw(6) << std::setfill('0') << i << std::setw(0) << ".bam";
	
		libmaus::aio::CheckedOutputStream::unique_ptr_type tCOS(new libmaus::aio::CheckedOutputStream(ostr.str()));
		COS[i] = UNIQUE_PTR_MOVE(tCOS);
		libmaus::bambam::BamWriter::unique_ptr_type twriter(new libmaus::bambam::BamWriter(*COS[i],*uphead,level));
		writers[i] = UNIQUE_PTR_MOVE(twriter);
	}

	uint64_t c = 0;
	if ( verbose )
	{
		while ( bamdec.readAlignment() )
		{
			algn.serialise ( writers [ (c++) % div ] -> getStream() );
			
			if ( ((c) & ((1ull<<20)-1)) == 0 )
				std::cerr << "[V] " << c << std::endl;
		}
		std::cerr << "[V] " << c << std::endl;
	}
	else
	{
		while ( bamdec.readAlignment() )
			algn.serialise ( writers [ (c++) % div ] -> getStream() );
	}

	for ( uint64_t i = 0; i < div; ++i )		
	{
		writers[i].reset();
		COS[i]->flush();
		COS[i].reset();
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
			
				V.push_back ( std::pair<std::string,std::string> ( "div=<["+::biobambam::Licensing::formatNumber(getDefaultDiv())+"]>", "divisor" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "prefix=<["+getDefaultFilePrefix(arginfo)+"]>", "default output file prefix" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamsplitmod(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

