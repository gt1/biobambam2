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
#include <libmaus2/bambam/BamAlignmentNameComparator.hpp>
#include <libmaus2/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamEntryContainer.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/BamHeaderUpdate.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/GetObject.hpp>
#include <libmaus2/util/PutObject.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }

#include <biobambam2/BamBamConfig.hpp>

#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
#include <libmaus2/bambam/ScramDecoder.hpp>
#endif

#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static bool getDefaultDisableValidation() { return false; }
static std::string getDefaultInputFormat() { return "bam"; }

int bamfixmateinformation(::libmaus2::util::ArgInfo const & arginfo)
{
	::libmaus2::util::TempFileRemovalContainer::setup();
	
	bool const inputisstdin = (!arginfo.hasArg("I")) || (arginfo.getUnparsedValue("I","-") == "-");
	bool const outputisstdout = (!arginfo.hasArg("O")) || (arginfo.getUnparsedValue("O","-") == "-");

	if ( isatty(STDIN_FILENO) && inputisstdin && (arginfo.getValue<std::string>("inputformat","bam") != "sam") )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) && outputisstdout && (arginfo.getValue<std::string>("outputformat","bam") != "sam") )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const disablevalidation = arginfo.getValue<int>("disablevalidation",getDefaultDisableValidation());

	std::string const inputformat = arginfo.getUnparsedValue("inputformat",getDefaultInputFormat());

	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilenameout = tmpfilenamebase + "_bamfixmateinformation";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilenameout);

	// input decoder wrapper
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;
	if ( disablevalidation )
		dec.disableValidation();
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();
	
	if ( ::libmaus2::bambam::BamHeader::getSortOrderStatic(header.text) != "queryname" )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "bamfixmateinformation: file is not sorted by queryname" << std::endl;
		se.finish();
		throw se;	
	}

	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(
		libmaus2::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bamfixmateinformation",std::string(PACKAGE_VERSION))
	);

	/*
	 * start index/md5 callbacks
	 */
	std::string md5filename;

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

	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	// construct writer
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(
		libmaus2::bambam::BamBlockWriterBaseFactory::construct(*uphead,arginfo,Pcbs)
	);
	
	libmaus2::bambam::BamAlignment & curalgn = dec.getAlignment();
	libmaus2::bambam::BamAlignment prevalgn;
	bool prevalgnvalid = false;
	libmaus2::bambam::BamAuxFilterVector MQfilter;
	MQfilter.set("MQ");
	uint64_t c = 0;
	
	while ( dec.readAlignment() )
	{
		if ( curalgn.isSecondary() || curalgn.isSupplementary() )
		{
			Pwriter->writeAlignment(curalgn);
		}
		else if ( prevalgnvalid )
		{
			// different name
			if ( strcmp(curalgn.getName(),prevalgn.getName()) )
			{
				Pwriter->writeAlignment(prevalgn);
				curalgn.swap(prevalgn);
			}
			// same name
			else
			{
				libmaus2::bambam::BamAlignment::fixMateInformation(prevalgn,curalgn,MQfilter);
				Pwriter->writeAlignment(prevalgn);
				Pwriter->writeAlignment(curalgn);
				prevalgnvalid = false;
			}
		}
		else
		{
			prevalgn.swap(curalgn);
			prevalgnvalid = true;
		}
		
		if ( verbose && ( ( ++c & ((1ull<<20)-1) ) == 0 ) )
			std::cerr << "[V] " << c << std::endl;
	}
	
	if ( prevalgnvalid )
	{
		Pwriter->writeAlignment(prevalgn);
		prevalgnvalid = false;
	}

	// flush encoder so callbacks see all output data
	Pwriter.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	
	if ( verbose )
		std::cerr << "[V] " << c << std::endl;	

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
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam2::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus2::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus2::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<1>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamfixmateinformation(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
