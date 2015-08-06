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

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamCat.hpp>
#include <libmaus2/bambam/BamWriter.hpp>

#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static uint64_t getDefaultN() { return 64*1024; }
static std::string getDefaultFilePrefix(::libmaus2::util::ArgInfo const & arginfo) { return arginfo.getDefaultTmpFileName(); }

::libmaus2::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus2::util::ArgInfo const & arginfo,
	::libmaus2::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamsplit", // ID
		"bamsplit", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus2::bambam::BamHeader(upheadtext));
	
	return UNIQUE_PTR_MOVE(uphead);
}

int bamsplit(libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	uint64_t const n = arginfo.getValue<int>("n",getDefaultN());
	std::string const prefix = arginfo.getUnparsedValue("prefix",getDefaultFilePrefix(arginfo));

	libmaus2::bambam::BamDecoder bamdec(std::cin);
	libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
	libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));

	libmaus2::aio::OutputStreamInstance::unique_ptr_type COS;
	libmaus2::bambam::BamWriter::unique_ptr_type writer;
	
	uint64_t c = 0;
	uint64_t f = 0;
	while ( bamdec.readAlignment() )
	{
		if ( c++ % n == 0 )
		{
			writer.reset();
			if ( COS )
				COS->flush();
			COS.reset();
			
			std::ostringstream fnostr;
			fnostr << prefix << "_" << std::setw(6) << std::setfill('0') << f++ << std::setw(0) << ".bam";
			std::string const fn = fnostr.str();
			
			libmaus2::aio::OutputStreamInstance::unique_ptr_type tCOS(new libmaus2::aio::OutputStreamInstance(fn));
			COS = UNIQUE_PTR_MOVE(tCOS);
			
			libmaus2::bambam::BamWriter::unique_ptr_type twriter(new libmaus2::bambam::BamWriter(*COS,*uphead,level));
			writer = UNIQUE_PTR_MOVE(twriter);
			
			if ( verbose )
				std::cerr << "[V] opened file " << fn << std::endl;
		}
		
		algn.serialise(writer->getStream());
	}
	
	writer.reset();
	if ( COS )
		COS->flush();
	COS.reset();

	return EXIT_SUCCESS;
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
			
				V.push_back ( std::pair<std::string,std::string> ( "n=<["+::biobambam2::Licensing::formatNumber(getDefaultN())+"]>", "number of entries" ) );
				V.push_back ( std::pair<std::string,std::string> ( "prefix=<["+getDefaultFilePrefix(arginfo)+"]>", "default output file prefix" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamsplit(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

