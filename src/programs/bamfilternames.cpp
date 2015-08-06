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

#include <libmaus2/aio/CheckedOutputStream.hpp>

#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
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

int bamfilternames(::libmaus2::util::ArgInfo const & arginfo)
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

	::libmaus2::bambam::BamDecoder dec(std::cin,false);
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamfilternames", // ID
		"bamfilternames", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	libmaus2::bambam::BamHeader const uphead(upheadtext);

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
	::libmaus2::bambam::BamWriter::unique_ptr_type writer(new ::libmaus2::bambam::BamWriter(std::cout,uphead,level,Pcbs));

	::libmaus2::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type LHTsnofailure;
	
	if ( arginfo.hasArg("names") )
	{
		std::string const names = arginfo.getUnparsedValue("names",std::string());
		
		libmaus2::aio::InputStreamInstance namestr(names);
		std::vector<std::string> vnames;
		while ( namestr )
		{
			std::string line;
			std::getline(namestr,line);
			
			if ( line.size() )
				vnames.push_back(line);
		}

		::libmaus2::trie::Trie<char> trienofailure;
		trienofailure.insertContainer(vnames);
		::libmaus2::trie::LinearHashTrie<char,uint32_t>::unique_ptr_type LHTnofailure(trienofailure.toLinearHashTrie<uint32_t>());
		LHTsnofailure = ::libmaus2::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type(LHTnofailure.release());
	}
	
	libmaus2::bambam::BamAlignment & algn = dec.getAlignment();
	uint64_t c = 0;
	
	if ( ! LHTsnofailure )
	{
		while ( dec.readAlignment() )
		{
			algn.serialise(writer->getStream());
			if ( verbose && (++c & (1024*1024-1)) == 0 )
				std::cerr << "[V] " << c/(1024*1024) << std::endl;
		}
	}
	else
	{
		while ( dec.readAlignment() )
		{
			if ( LHTsnofailure->searchCompleteNoFailureZ(algn.getName()) != -1 )
				algn.serialise(writer->getStream());
			if ( verbose && (++c & (1024*1024-1)) == 0 )
				std::cerr << "[V] " << c/(1024*1024) << std::endl;
		}	
	}

	writer.reset();

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
				V.push_back ( std::pair<std::string,std::string> ( "names=<[]>", "file containing read names to be kept (default: keep all)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
								
				return EXIT_SUCCESS;
			}
			
		return bamfilternames(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

