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
#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <iomanip>

#include <config.h>

#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamParallelRewrite.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

#if 0
::libmaus::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus::util::ArgInfo const & arginfo,
	::libmaus::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamfilterflags", // ID
		"bamfilterflags", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus::bambam::BamHeader(upheadtext));
	
	return UNIQUE_PTR_MOVE(uphead);
}

struct UpdateHeader : public libmaus::bambam::BamHeaderRewriteCallback
{
	libmaus::util::ArgInfo const & arginfo;

	UpdateHeader(libmaus::util::ArgInfo const & rarginfo)
	: arginfo(rarginfo)
	{
	
	}

	::libmaus::bambam::BamHeader::unique_ptr_type operator()(::libmaus::bambam::BamHeader const & header)  const
	{
		::libmaus::bambam::BamHeader::unique_ptr_type ptr(updateHeader(arginfo,header));
		return UNIQUE_PTR_MOVE(ptr);
	}
};
#endif

int bamfilterflags(::libmaus::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude",""));

	std::cerr << "[V] excluding " << excludeflags << std::endl;

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION));
	
	uint64_t const numthreads = arginfo.getValue<uint64_t>("numthreads",1);
	uint64_t cnt = 0;
	uint64_t kept = 0;

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string md5filename;
	std::string indexfilename;

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
			md5filename = arginfo.getUnparsedValue("md5filename","");
		else
			std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

		if ( md5filename.size() )
		{
			::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
			Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
			cbs.push_back(Pmd5cb.get());
		}
	}
	libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	if ( arginfo.getValue<unsigned int>("index",getDefaultIndex()) )
	{
		if ( arginfo.hasArg("indexfilename") &&  arginfo.getUnparsedValue("indexfilename","") != "" )
			indexfilename = arginfo.getUnparsedValue("indexfilename","");
		else
			std::cerr << "[V] no filename for index given, not creating index" << std::endl;

		if ( indexfilename.size() )
		{
			libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Tindex(new libmaus::bambam::BgzfDeflateOutputCallbackBamIndex(tmpfileindex));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
			cbs.push_back(Pindex.get());
		}
	}
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */
	 

	if ( numthreads == 1 )
	{
		::libmaus::bambam::BamDecoder BD(std::cin);
		::libmaus::bambam::BamHeader const & bamheader = BD.getHeader();
		::libmaus::bambam::BamHeader::unique_ptr_type uphead(libmaus::bambam::BamHeaderUpdate::updateHeader(arginfo,bamheader,"bamfilterflags",std::string(PACKAGE_VERSION)));
		::libmaus::bambam::BamAlignment & alignment = BD.getAlignment();
		::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(std::cout,*uphead,level,Pcbs));
	
		for ( ; BD.readAlignment(); ++cnt )
		{
			if ( cnt % (1024*1024) == 0 )
				std::cerr << "[V] processed " << cnt << " kept " << kept << " removed " << (cnt-kept) << std::endl;
			if ( ! (alignment.getFlags() & excludeflags) )			
			{
				alignment.serialise(writer->getStream());
				++kept;
			}
		}

		std::cerr << "[V] " << cnt << std::endl;	
	}
	else
	{
		::libmaus::bambam::BamHeaderUpdate UH(arginfo,"bamfilterflags",std::string(PACKAGE_VERSION));
		libmaus::bambam::BamParallelRewrite BPR(std::cin,UH,std::cout,Z_DEFAULT_COMPRESSION,numthreads,4 /* blocks per thread */,Pcbs);
		libmaus::bambam::BamAlignmentDecoder & dec = BPR.getDecoder();
		libmaus::bambam::BamParallelRewrite::writer_type & writer = BPR.getWriter();

		libmaus::bambam::BamAlignment const & algn = dec.getAlignment();
		for ( ; dec.readAlignment(); ++cnt )
		{
			if ( cnt % (1024*1024) == 0 )
				std::cerr << "[V] processed " << cnt << " kept " << kept << " removed " << (cnt-kept) << std::endl;
			if ( ! (algn.getFlags() & excludeflags) )
			{
				algn.serialise(writer.getStream());	
				++kept;
			}
		}
		
		std::cerr << "[V] " << cnt << std::endl;	
	}
	
	std::cerr << "[V] kept " << kept << " removed " << cnt-kept << std::endl;

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
		libmaus::timing::RealTimeClock rtc; rtc.start();
		
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
				std::cerr << ::biobambam::Licensing::license() << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
				
				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "exclude=<[]>", "exclude alignments matching any of the given flags" ) );
				V.push_back ( std::pair<std::string,std::string> ( "numthreads=<[1]>", "number of recoding threads" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				std::cerr << "Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" << std::endl;

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
		
			
		bamfilterflags(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
