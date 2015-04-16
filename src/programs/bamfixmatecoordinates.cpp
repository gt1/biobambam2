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
#include <config.h>

#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>
#include <libmaus/bambam/EncoderBase.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/lz/BufferedGzipStream.hpp>

#include <libmaus/fastx/CompactFastEncoder.hpp>
#include <libmaus/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus/util/GetObject.hpp>

#include <libmaus/lz/Deflate.hpp>
#include <libmaus/util/CountPutObject.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/timing/RealTimeClock.hpp>

#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }
static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static unsigned int getDefaultColHashBits() { return 20; }
static uint64_t getDefaultColListSize() { return 512*1024; }

static uint64_t getMapCnt(::libmaus::bambam::CollatingBamDecoder::alignment_ptr_type const & p)
{
	if ( p && p->isMapped() )
		return 1;
	else
		return 0;
}

int bamfixmatecoordinates(::libmaus::util::ArgInfo const & arginfo)
{
	::libmaus::util::TempFileRemovalContainer::setup();
	::libmaus::timing::RealTimeClock rtc; rtc.start();
	
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	unsigned int const colhashbits = arginfo.getValue<unsigned int>("colhashbits",getDefaultColHashBits());
	unsigned int const collistsize = arginfo.getValue<unsigned int>("collistsize",getDefaultColListSize());
	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	
	std::string const tmpfilename = tmpfilenamebase + "_bamcollate";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	
	::libmaus::bambam::CollatingBamDecoder CBD(std::cin,tmpfilename,false /* put rank */,colhashbits/*hash bits*/,collistsize/*size of output list*/);
	::libmaus::bambam::BamFormatAuxiliary auxdata;
	::libmaus::bambam::BamHeader const & bamheader = CBD.getHeader();
	
	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		bamheader.text,
		"bamfixmatecoordinates", // ID
		"bamfixmatecoordinates", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(bamheader.text).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);
	
	if ( uphead.getSortOrder() != "queryname" )
		uphead.changeSortOrder("unknown");

	/*
	 * start index/md5 callbacks
	 */
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
	
	// setup bam writer
	::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(std::cout,uphead,level,Pcbs));
	
	#if 0
	::libmaus::bambam::ProgramHeaderLineSet PHLS(bamheader.text);
	std::cerr << "Last id in PG chain: " << PHLS.getLastIdInChain() << std::endl;
	#endif

	// std::cout << bamheader.text;

	typedef ::libmaus::bambam::CollatingBamDecoder::alignment_ptr_type alignment_ptr_type;
	std::pair<alignment_ptr_type,alignment_ptr_type> P;
	uint64_t const mod = 1024*1024;
	uint64_t proc = 0;
	uint64_t lastproc = 0;
	uint64_t paircnt = 0;
	
	while ( CBD.tryPair(P) )
	{
		uint64_t const mapcnt = getMapCnt(P.first) + getMapCnt(P.second);
		
		if ( mapcnt == 1 )
		{
			int32_t refid = -1;
			int32_t pos = -1;
			
			if ( P.first )
			{
				refid = P.first->getRefID();
				pos = P.first->getPos();
			}
			else
			{
				assert ( P.second );

				refid = P.second->getRefID();
				pos = P.second->getPos();
			}
			
			P.first->putRefId(refid);
			P.first->putPos(pos);
			P.first->putNextRefId(refid);
			P.first->putNextPos(pos);
			P.second->putRefId(refid);
			P.second->putPos(pos);
			P.second->putNextRefId(refid);
			P.second->putNextPos(pos);
		}
		
		if ( P.first )
		{
			P.first->serialise(writer->getStream());
			++proc;
		}
		if ( P.second )
		{
			P.second->serialise(writer->getStream());
			++proc;
		}
		if ( P.first && P.second )
		{
			paircnt++;
		}
		
		if ( verbose && (proc/mod != lastproc/mod) )
		{
			std::cerr 
				<< "Processed " << proc << " fragments, " << paircnt << " pairs, " 
				<< proc/rtc.getElapsedSeconds() << " al/s"
				<< std::endl;
			lastproc = proc;
		}
	}		

	if ( verbose )
		std::cerr 	
			<< "Processed " << proc << " fragments, " << paircnt << " pairs, " 
			<< proc/rtc.getElapsedSeconds() << " al/s"
			<< std::endl;

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
				
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhashbits=<["+::biobambam::Licensing::formatNumber(getDefaultColHashBits())+"]>", "log_2 of size of hash table used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collistsize=<["+::biobambam::Licensing::formatNumber(getDefaultColListSize())+"]>", "output list size for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
		
		return bamfixmatecoordinates(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
