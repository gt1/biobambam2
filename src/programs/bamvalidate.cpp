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
#include "config.h"

#include <iostream>

#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

int getDefaultMD5()
{
	return 0;
}
int getDefaultIndex()
{
	return 0;
}

int getDefaultLevel()
{
	return Z_DEFAULT_COMPRESSION;
}

bool getDefaultPassThrough()
{
	return false;
}

int getDefaultVerbose()
{
	return false;
}

char const * getDefaultInputFormat()
{
	return "bam";
}

template<bool passthrough>
int bamvalidateTemplate(::libmaus::util::ArgInfo const & arginfo)
{
	libmaus::timing::RealTimeClock rtc; rtc.start();	
	bool const verbose = arginfo.getValue("verbose",getDefaultVerbose());
	
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamHeader const & header = dec.getHeader();	
	::libmaus::bambam::BamAlignment const & algn = dec.getAlignment();

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		header.text,
		"bamvalidate", // ID
		"bamvalidate", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(header.text).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);

	/*
	 * start index/md5 callbacks and alignment writer
	 */
	std::string md5filename;
	std::string indexfilename;

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type Pout;

	if ( passthrough )
	{
		std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
		std::string const tmpfileindex = tmpfilenamebase + "_index";
		::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);
		
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


		libmaus::bambam::BamBlockWriterBase::unique_ptr_type Tout ( 
			libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, Pcbs) 
		);
		Pout = UNIQUE_PTR_MOVE(Tout);
	}
	
	libmaus::autoarray::AutoArray<char> lastvalidname(256); // max valid read name is 255 bytes
	uint64_t alsok = 0;

	try
	{
		while ( dec.readAlignment() )
		{
			if ( passthrough )
				Pout->writeAlignment(algn);

			uint64_t const lname = algn.getLReadName();
			char const * name = algn.getName();
			std::copy(name,name+lname+1,lastvalidname.begin());
			
			alsok += 1;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << "[E] name of last valid alignment was " << lastvalidname.begin() << std::endl;
		std::cerr << "[E] read " << alsok << " valid alignments" << std::endl;
		throw;
	}
	
	Pout.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
	}

	if ( verbose )
		std::cerr << "[V] checked " << alsok << " alignments in " << rtc.formatTime(rtc.getElapsedSeconds()) 
			<< " (" << alsok / rtc.getElapsedSeconds() << " al/s)" << std::endl;
	
	return EXIT_SUCCESS;
}

int bamvalidate(::libmaus::util::ArgInfo const & arginfo)
{
	bool const passthrough = arginfo.getValue<unsigned int>("passthrough",getDefaultPassThrough());
	
	if ( passthrough )
		return bamvalidateTemplate<true>(arginfo);
	else
		return bamvalidateTemplate<false>(arginfo);
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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print stats at the end of a successfull run" ) );
				V.push_back ( std::pair<std::string,std::string> ( "passthrough=<["+::biobambam::Licensing::formatNumber(getDefaultPassThrough())+"]>", "write alignments to standard output (default: do not pass through)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory (passthrough=1, index=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0, passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0, passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ", passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1, passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset, passthrough=1 only)" ) );


				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}


		return bamvalidate(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
