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

#include <libmaus/aio/CheckedOutputStream.hpp>

#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamAlignmentNameComparator.hpp>
#include <libmaus/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamEntryContainer.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GetObject.hpp>
#include <libmaus/util/PutObject.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

#include <biobambam/BamBamConfig.hpp>

#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
#include <libmaus/bambam/ScramDecoder.hpp>
#endif

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static std::string getDefaultSortOrder() { return "coordinate"; }
static uint64_t getDefaultBlockSize() { return 1024; }
static bool getDefaultDisableValidation() { return false; }
static std::string getDefaultInputFormat() { return "bam"; }
static int getDefaultFixMates() { return 0; }
static int getDefaultSortThreads() { return 1; }
static int getDefaultCalMdNm() { return 0; }
static int getDefaultCalMdNmRecompIndetOnly() { return 0; }
static int getDefaultCalMdNmWarnChange() { return 0; }
static int getDefaultAddDupMarkSupport() { return 0; }

int bamsort(::libmaus::util::ArgInfo const & arginfo)
{
	::libmaus::util::TempFileRemovalContainer::setup();
	
	bool const inputisstdin = (!arginfo.hasArg("I")) || (arginfo.getUnparsedValue("I","-") == "-");
	bool const outputisstdout = (!arginfo.hasArg("O")) || (arginfo.getUnparsedValue("O","-") == "-");

	if ( isatty(STDIN_FILENO) && inputisstdin && (arginfo.getValue<std::string>("inputformat","bam") != "sam") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) && outputisstdout && (arginfo.getValue<std::string>("outputformat","bam") != "sam") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const disablevalidation = arginfo.getValue<int>("disablevalidation",getDefaultDisableValidation());

	std::string const inputformat = arginfo.getUnparsedValue("inputformat",getDefaultInputFormat());

	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilenameout = tmpfilenamebase + "_bamsort";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenameout);
	uint64_t blockmem = arginfo.getValue<uint64_t>("blockmb",getDefaultBlockSize())*1024*1024;
	std::string const sortorder = arginfo.getValue<std::string>("SO","coordinate");
	bool const fixmates = arginfo.getValue<int>("fixmates",getDefaultFixMates());
	bool const addMSMC = arginfo.getValue<int>("adddupmarksupport",getDefaultAddDupMarkSupport());
	uint64_t sortthreads = arginfo.getValue<uint64_t>("sortthreads",getDefaultSortThreads());

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	if ( disablevalidation )
		dec.disableValidation();
	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamsort", // ID
		"bamsort", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);

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
	if ( sortorder != "queryname" )
		uphead.changeSortOrder("coordinate");
	else
		uphead.changeSortOrder("queryname");

	libmaus::bambam::BamBlockWriterBase::unique_ptr_type Pout ( libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, Pcbs) );

	if ( fixmates )
	{
		if ( sortorder != "queryname" )
		{
			::libmaus::bambam::BamEntryContainer< ::libmaus::bambam::BamAlignmentPosComparator > 
				BEC(blockmem,tmpfilenameout,sortthreads);

			if ( verbose )
				std::cerr << "[V] Reading alignments from source." << std::endl;
			uint64_t incnt = 0;

			// current alignment
			libmaus::bambam::BamAlignment & curalgn = dec.getAlignment();
			// previous alignment
			libmaus::bambam::BamAlignment prevalgn;
			// previous alignment valid
			bool prevalgnvalid = false;
			// MQ field filter
			libmaus::bambam::BamAuxFilterVector MQfilter;
			libmaus::bambam::BamAuxFilterVector MSfilter;
			libmaus::bambam::BamAuxFilterVector MCfilter;
			MQfilter.set("MQ");
			MSfilter.set("MS");
			MCfilter.set("MC");
			
			while ( dec.readAlignment() )
			{
				if ( curalgn.isSecondary() || curalgn.isSupplementary() )
				{
					BEC.putAlignment(curalgn);
				}
				else if ( prevalgnvalid )
				{
					// different name
					if ( strcmp(curalgn.getName(),prevalgn.getName()) )
					{
						BEC.putAlignment(prevalgn);
						curalgn.swap(prevalgn);
					}
					// same name
					else
					{
						libmaus::bambam::BamAlignment::fixMateInformation(prevalgn,curalgn,MQfilter);

						if ( addMSMC )
						{
							libmaus::bambam::BamAlignment::addMateBaseScore(prevalgn,curalgn,MSfilter);
							libmaus::bambam::BamAlignment::addMateCoordinate(prevalgn,curalgn,MCfilter);
						}

						BEC.putAlignment(prevalgn);
						BEC.putAlignment(curalgn);
						prevalgnvalid = false;
					}
				}
				else
				{
					prevalgn.swap(curalgn);
					prevalgnvalid = true;
				}
				
				if ( verbose && ( ( ++incnt & ((1ull<<20)-1) ) == 0 ) )
					std::cerr << "[V] " << incnt << std::endl;
			}
			
			if ( prevalgnvalid )
			{
				BEC.putAlignment(prevalgn);
				prevalgnvalid = false;
			}

			if ( verbose )
				std::cerr << "[V] read " << incnt << " alignments" << std::endl;

			// BEC.createOutput(std::cout, uphead, level, verbose, Pcbs);
			BEC.createOutput(*Pout, verbose);
		}
		else
		{
			::libmaus::bambam::BamEntryContainer< ::libmaus::bambam::BamAlignmentNameComparator > 
				BEC(blockmem,tmpfilenameout,sortthreads);
			
			if ( verbose )
				std::cerr << "[V] Reading alignments from source." << std::endl;
			uint64_t incnt = 0;
			
			// current alignment
			libmaus::bambam::BamAlignment & curalgn = dec.getAlignment();
			// previous alignment
			libmaus::bambam::BamAlignment prevalgn;
			// previous alignment valid
			bool prevalgnvalid = false;
			// MQ field filter
			libmaus::bambam::BamAuxFilterVector MQfilter;
			libmaus::bambam::BamAuxFilterVector MSfilter;
			libmaus::bambam::BamAuxFilterVector MCfilter;
			MQfilter.set("MQ");
			MSfilter.set("MS");
			MCfilter.set("MC");
			
			while ( dec.readAlignment() )
			{
				if ( curalgn.isSecondary() || curalgn.isSupplementary() )
				{
					BEC.putAlignment(curalgn);
				}
				else if ( prevalgnvalid )
				{
					// different name
					if ( strcmp(curalgn.getName(),prevalgn.getName()) )
					{
						BEC.putAlignment(prevalgn);
						curalgn.swap(prevalgn);
					}
					// same name
					else
					{
						libmaus::bambam::BamAlignment::fixMateInformation(prevalgn,curalgn,MQfilter);

						if ( addMSMC )
						{
							libmaus::bambam::BamAlignment::addMateBaseScore(prevalgn,curalgn,MSfilter);
							libmaus::bambam::BamAlignment::addMateCoordinate(prevalgn,curalgn,MCfilter);
						}

						BEC.putAlignment(prevalgn);
						BEC.putAlignment(curalgn);
						prevalgnvalid = false;
					}
				}
				else
				{
					prevalgn.swap(curalgn);
					prevalgnvalid = true;
				}
				
				if ( verbose && ( ( ++incnt & ((1ull<<20)-1) ) == 0 ) )
					std::cerr << "[V] " << incnt << std::endl;
			}
			
			if ( prevalgnvalid )
			{
				BEC.putAlignment(prevalgn);
				prevalgnvalid = false;
			}
			
			if ( verbose )
				std::cerr << "[V] read " << incnt << " alignments" << std::endl;

			// BEC.createOutput(std::cout, uphead, level, verbose, Pcbs);
			BEC.createOutput(*Pout, verbose);
		}
	}
	else
	{
		if ( sortorder != "queryname" )
		{
			bool const calmdnm = arginfo.getValue<unsigned int>("calmdnm",getDefaultCalMdNm());			
			if ( calmdnm && (! arginfo.hasArg("calmdnmreference")) )
			{
				libmaus::exception::LibMausException lme;
				lme.getStream() << "calmdnm is set but required calmdnmreference is not, aborting." << std::endl;
				lme.finish();
				throw lme;
			}
			std::string const calmdnmreference = arginfo.getUnparsedValue("calmdnmreference","");
			bool const calmdnmrecompindetonly = arginfo.getValue<unsigned int>("calmdnmrecompindetonly",getDefaultCalMdNmRecompIndetOnly());
			bool const calmdnmwarnchange = arginfo.getValue<unsigned int>("calmdnmwarnchange",getDefaultCalMdNmWarnChange());
			
			::libmaus::bambam::BamEntryContainer< ::libmaus::bambam::BamAlignmentPosComparator > BEC(blockmem,tmpfilenameout,sortthreads);

			if ( verbose )
				std::cerr << "[V] Reading alignments from source." << std::endl;
			uint64_t incnt = 0;
			
			while ( dec.readAlignment() )
			{
				BEC.putAlignment(dec.getAlignment());
				incnt++;
				if ( verbose && (incnt % (1024*1024) == 0) )
					std::cerr << "[V] " << incnt/(1024*1024) << "M" << std::endl;
			}

			if ( verbose )
				std::cerr << "[V] read " << incnt << " alignments" << std::endl;


			if ( calmdnm )
			{
				libmaus::bambam::MdNmRecalculation mdnmrecalc(calmdnmreference,false /* do not validate again */,calmdnmrecompindetonly,calmdnmwarnchange,64*1024);
				BEC.createOutput(*Pout, verbose, &mdnmrecalc);
			}
			else
			{
				BEC.createOutput(*Pout, verbose, 0);			
			}
		}
		else
		{
			::libmaus::bambam::BamEntryContainer< ::libmaus::bambam::BamAlignmentNameComparator > BEC(blockmem,tmpfilenameout,sortthreads);
			
			if ( verbose )
				std::cerr << "[V] Reading alignments from source." << std::endl;
			uint64_t incnt = 0;
			
			while ( dec.readAlignment() )
			{
				BEC.putAlignment(dec.getAlignment());
				incnt++;
				if ( verbose && (incnt % (1024*1024) == 0) )
					std::cerr << "[V] " << incnt/(1024*1024) << "M" << std::endl;
			}
			
			if ( verbose )
				std::cerr << "[V] read " << incnt << " alignments" << std::endl;

			// BEC.createOutput(std::cout, uphead, level, verbose, Pcbs);
			BEC.createOutput(*Pout, verbose);
		}
	}

	// flush encoder so callbacks see all output data
	Pout.reset();

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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "SO=<["+getDefaultSortOrder()+"]>", "sorting order (coordinate or queryname)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "blockmb=<["+::biobambam::Licensing::formatNumber(getDefaultBlockSize())+"]>", "size of internal memory buffer used for sorting in MiB" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("fixmates=<[")+::biobambam::Licensing::formatNumber(getDefaultFixMates())+"]>", "fix mate information (for name collated input only, disabled by default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("calmdnm=<[")+::biobambam::Licensing::formatNumber(getDefaultCalMdNm())+"]>", "calculate MD and NM aux fields (for coordinate sorted output only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("calmdnmreference=<[]>"), "reference for calculating MD and NM aux fields (calmdnm=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("calmdnmrecompindetonly=<[")+::biobambam::Licensing::formatNumber(getDefaultCalMdNm())+"]>", "only recalculate MD and NM in the presence of indeterminate bases (calmdnm=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("calmdnmwarnchange=<[")+::biobambam::Licensing::formatNumber(getDefaultCalMdNmWarnChange())+"]>", "warn when changing existing MD/NM fields (calmdnm=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("adddupmarksupport=<[")+::biobambam::Licensing::formatNumber(getDefaultAddDupMarkSupport())+"]>", "add info for streaming duplicate marking (for name collated input only, ignored for fixmate=0, disabled by default)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamsort(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

