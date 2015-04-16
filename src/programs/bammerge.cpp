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
#include <libmaus2/bambam/BamMergeCoordinate.hpp>
#include <libmaus2/bambam/BamMergeQueryName.hpp>
#include <libmaus2/bambam/BamWriter.hpp>

#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static std::string getDefaultSortOrder() { return "coordinate"; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

#if defined(LIBMAUS2_HAVE_IRODS)
#include <libmaus2/irods/IRodsInputStreamFactory.hpp>
#endif

::libmaus2::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus2::util::ArgInfo const & arginfo,
	::libmaus2::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bammerge", // ID
		"bammerge", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus2::bambam::BamHeader(upheadtext));
	
	return UNIQUE_PTR_MOVE(uphead);
}

int bammerge(libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	std::string const sortorder = arginfo.getValue<std::string>("SO",getDefaultSortOrder());

	std::vector<std::string> inputfilenames = arginfo.getPairValues("I");
	for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
		inputfilenames.push_back(arginfo.restargs[i]);
	std::vector<std::string> inputmetafilenames = arginfo.getPairValues("IL");
	for ( uint64_t i = 0; i < inputmetafilenames.size(); ++i )
	{
		libmaus2::aio::InputStream::unique_ptr_type iptr(libmaus2::aio::InputStreamFactoryContainer::constructUnique(inputmetafilenames[i]));
		std::istream & in = *iptr;
		while ( in )
		{
			std::string line;
			std::getline(in,line);
			if ( line.size() )
				inputfilenames.push_back(line);
		}
	}

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string md5filename;
	std::string indexfilename;

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
	libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	if ( arginfo.getValue<unsigned int>("index",getDefaultIndex()) )
	{
		if ( arginfo.hasArg("indexfilename") &&  arginfo.getUnparsedValue("indexfilename","") != "" )
			indexfilename = arginfo.getUnparsedValue("indexfilename","");
		else
			std::cerr << "[V] no filename for index given, not creating index" << std::endl;

		if ( indexfilename.size() )
		{
			libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Tindex(new libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex(tmpfileindex));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
			cbs.push_back(Pindex.get());
		}
	}
	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	if ( sortorder == "queryname" )
	{
		libmaus2::bambam::BamMergeQueryName bamdec(arginfo,inputfilenames /* ,true */);
		libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
		libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
		::libmaus2::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));
		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(
			libmaus2::bambam::BamBlockWriterBaseFactory::construct(*uphead,arginfo,Pcbs));
		
		if ( verbose )
		{
			uint64_t c = 0;
			while ( bamdec.readAlignment() )
			{
				Pwriter->writeAlignment(algn);
	
				if ( ((++c) & ((1ull<<20)-1)) == 0 )
					std::cerr << "[V] " << c << std::endl;
			}
		
			std::cerr << "[V] " << c << std::endl;
		}
		else
			while ( bamdec.readAlignment() )
				Pwriter->writeAlignment(algn);
	}
	else
	{
		libmaus2::bambam::BamMergeCoordinate bamdec(arginfo,inputfilenames /* ,true */);
		libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
		libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
		::libmaus2::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));
		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(
			libmaus2::bambam::BamBlockWriterBaseFactory::construct(*uphead,arginfo,Pcbs));
			
		if ( verbose )
		{
			uint64_t c = 0;

			while ( bamdec.readAlignment() )
			{
				Pwriter->writeAlignment(algn);
	
				if ( ((++c) & ((1ull<<20)-1)) == 0 )
					std::cerr << "[V] " << c << std::endl;
			}
		
			std::cerr << "[V] " << c << std::endl;
		}
		else
		{
			while ( bamdec.readAlignment() )
				Pwriter->writeAlignment(algn);
		}
	}

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		#if defined(LIBMAUS2_HAVE_IRODS)
                libmaus2::irods::IRodsInputStreamFactory::registerHandler();
                #endif

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
			
				V.push_back ( std::pair<std::string,std::string> ( "I=<[filename]>", "input file, can be set multiple times" ) );
				V.push_back ( std::pair<std::string,std::string> ( "SO=<["+getDefaultSortOrder()+"]>]", "sort order (coordinate or queryname)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

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

