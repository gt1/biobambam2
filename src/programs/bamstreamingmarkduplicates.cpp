/**
    biobambam
    Copyright (C) 2009-2014 German Tischler
    Copyright (C) 2011-2014 Genome Research Limited

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

// local
#include <config.h>

// std
#include <algorithm>
#include <iostream>
#include <vector>

// libmaus
#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamStreamingMarkDuplicates.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return true; }
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }
static int getDefaultResetDupFlag() { return 0; }

int bamstreamingmarkduplicates(libmaus::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue<uint64_t>("verbose",getDefaultVerbose());
	bool const resetdupflag = arginfo.getValue<uint64_t>("resetdupflag",getDefaultResetDupFlag());
	std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpfile",arginfo.getDefaultTmpFileName());	

	libmaus::aio::PosixFdInputStream PFIS(STDIN_FILENO);
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,true /* put rank */, 0 /* copy stream */, PFIS
		)
	);
	libmaus::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();
	libmaus::bambam::BamAlignment & algn = dec.getAlignment();

	libmaus::bambam::BamHeader const & header = dec.getHeader();

	::libmaus::bambam::BamHeader::unique_ptr_type genuphead(
		libmaus::bambam::BamHeaderUpdate::updateHeader(arginfo,dec.getHeader(),"bamstreamingmarkduplicates",std::string(PACKAGE_VERSION))
	);

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

	// construct writer
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(libmaus::bambam::BamBlockWriterBaseFactory::construct(*genuphead,arginfo,Pcbs));
	libmaus::bambam::BamBlockWriterBase & wr = *Pwriter;

	libmaus::bambam::BamStreamingMarkDuplicates BSMD(arginfo,header,wr);

	uint64_t cnt = 0;

	libmaus::timing::RealTimeClock globalrtc;
	globalrtc.start();
	libmaus::timing::RealTimeClock batchrtc;
	batchrtc.start();

	uint32_t const flagmask = ~static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);
				
	while ( dec.readAlignment() )
	{
		if ( resetdupflag )
		{
			uint32_t const flags = algn.getFlags();
			algn.putFlags(flags&flagmask);
		}
	
		BSMD.addAlignment(algn);
		
		if ( verbose && ((++cnt % (1024*1024)) == 0) )
		{
			std::cerr 
				<< "[V] " 
				<< cnt << " " 
				<< BSMD.OQ.nextout << " " 
				<< libmaus::util::MemUsage() << " "
				<< globalrtc.formatTime(globalrtc.getElapsedSeconds()) << " "
				<< batchrtc.formatTime(batchrtc.getElapsedSeconds())
				<< std::endl;
			
			batchrtc.start();
		}
	}
	
	BSMD.flush();

	// reset BAM writer
	Pwriter.reset();

	if ( Pmd5cb )
		Pmd5cb->saveDigestAsFile(md5filename);
	if ( Pindex )
		Pindex->flush(std::string(indexfilename));

	// write metrics
	BSMD.writeMetrics(arginfo);
	
	return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
	try
	{
		libmaus::util::ArgInfo const arginfo(argc,argv);

		
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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+libmaus::bambam::BamAlignmentDecoderInfo::getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				#if 0
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );

				V.push_back ( 
					std::pair<std::string,std::string> ( 
						std::string("maxreadlen=<[")+::biobambam::Licensing::formatNumber(libmaus::bambam::BamStreamingMarkDuplicates::getDefaultMaxReadLen())+"]>", 
						std::string("maximum read length (default ") + ::biobambam::Licensing::formatNumber(libmaus::bambam::BamStreamingMarkDuplicates::getDefaultMaxReadLen()) + ")" 
					)
				);
				V.push_back ( 
					std::pair<std::string,std::string> ( 
						std::string("optminpixeldif=<[")+::biobambam::Licensing::formatNumber(libmaus::bambam::BamStreamingMarkDuplicates::getDefaultOptMinPixelDif())+"]>", 
						std::string("maximum distance for optical duplicates (default ") + ::biobambam::Licensing::formatNumber(libmaus::bambam::BamStreamingMarkDuplicates::getDefaultOptMinPixelDif()) + ")" 
					)
				);
				V.push_back ( 
					std::pair<std::string,std::string> ( 
						std::string("resetdupflag=<[")+::biobambam::Licensing::formatNumber(getDefaultResetDupFlag())+"]>", 
						std::string("reset dup flag before checking for duplicates (default ") + ::biobambam::Licensing::formatNumber(getDefaultResetDupFlag()) + ")" 
					)
				);

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}


		return bamstreamingmarkduplicates(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}
}
