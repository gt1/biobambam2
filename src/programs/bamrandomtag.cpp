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

#include <libmaus2/aio/OutputStreamInstance.hpp>

#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamAlignmentNameComparator.hpp>
#include <libmaus2/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamEntryContainer.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/GetObject.hpp>
#include <libmaus2/util/PutObject.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

int bamrandomtag(::libmaus2::util::ArgInfo const & arginfo)
{
	::libmaus2::util::TempFileRemovalContainer::setup();

	if ( isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());

	// tag field
	bool const havetag = arginfo.hasArg("tag");
	std::string const tag = arginfo.getUnparsedValue("tag","no tag");
	
	if ( (!havetag)
		||
			(
				tag.size() != 2 
				||
				(!isalpha(tag[0]))
				||
				(!isalnum(tag[1]))
			)
	)
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "mandatory argument tag is invalid" << std::endl;
		se.finish();
		throw se;			
	}

	
	::libmaus2::bambam::BamDecoder dec(std::cin,false);
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamrandomtag", // ID
		"bamrandomtag", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	libmaus2::bambam::BamHeader const uphead(upheadtext);

	::libmaus2::bambam::BamWriter::unique_ptr_type writer(new ::libmaus2::bambam::BamWriter(std::cout,uphead,level));
 	
	libmaus2::bambam::BamAlignment & algn = dec.getAlignment();
	uint64_t c = 0;

	libmaus2::bambam::BamAuxFilterVector tagfilter;
	tagfilter.set(tag.c_str());
	char const * taga = "A";
	char const * tagb = "C";
	srand(time(0));

	while ( dec.readAlignment() )
	{
		algn.filterOutAux(tagfilter);
		algn.putAuxString(tag.c_str(),(rand()&1) ? taga : tagb);
		algn.serialise(writer->getStream());
 			
		if ( verbose && (++c & (1024*1024-1)) == 0 )
			std::cerr << "[V] " << c/(1024*1024) << std::endl;
	}

	writer.reset();

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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
								
				return EXIT_SUCCESS;
			}
			
		return bamrandomtag(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

