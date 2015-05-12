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

#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

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

bool getDefaultBaseQualHist()
{
	return false;
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
int bamvalidateTemplate(::libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::timing::RealTimeClock rtc; rtc.start();	
	bool const verbose = arginfo.getValue("verbose",getDefaultVerbose());
	bool const basequalhist = arginfo.getValue("basequalhist",getDefaultBaseQualHist());
	
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();	
	::libmaus2::bambam::BamAlignment const & algn = dec.getAlignment();

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		header.text,
		"bamvalidate", // ID
		"bamvalidate", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(header.text).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus2::bambam::BamHeader uphead(upheadtext);

	/*
	 * start index/md5 callbacks and alignment writer
	 */
	std::string md5filename;
	std::string indexfilename;

	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pout;

	if ( passthrough )
	{
		std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
		std::string const tmpfileindex = tmpfilenamebase + "_index";
		::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);
		
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


		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Tout ( 
			libmaus2::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, Pcbs) 
		);
		Pout = UNIQUE_PTR_MOVE(Tout);
	}
	
	libmaus2::autoarray::AutoArray<char> lastvalidname(256); // max valid read name is 255 bytes
	uint64_t alsok = 0;

	::libmaus2::autoarray::AutoArray<char> qual;
	libmaus2::autoarray::AutoArray<uint64_t> H(static_cast<uint64_t>(std::numeric_limits<uint8_t>::max())+1);
	std::fill(H.begin(),H.end(),0ull);
	
	try
	{
		while ( dec.readAlignment() )
		{
			if ( passthrough )
				Pout->writeAlignment(algn);

			if ( basequalhist )
			{
				uint64_t const l = algn.getLseq();
				uint8_t const * Qc = libmaus2::bambam::BamAlignmentDecoderBase::getQual(algn.D.begin());					
				uint8_t const * const Qe = Qc + l;
				
				while ( Qc != Qe )
					H[*(Qc++)]++;
			}

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
			
	if ( basequalhist )
	{
		uint64_t const s = std::accumulate(H.begin(),H.end(),0ull);
		
		uint64_t a = 0;		
		uint64_t minq = std::numeric_limits<uint64_t>::max();
		uint64_t maxq = 0;
		
		for ( uint64_t i = 0; i < H.size(); ++i )
			if ( H[i] )
			{
				minq = std::min(minq,i);
				maxq = std::max(maxq,i);
			
				a += H[i];
				
				std::cerr 
					<< "[H]\t" << i << "\t";
				
				if ( ( static_cast<uint64_t>(i+33) < static_cast<uint64_t>(std::numeric_limits<char>::max()) && isprint(i+33)) )
					std::cerr << static_cast<char>(i+33);
				
				std::cerr << "\t"
					<< H[i] << "\t" 
					<< (H[i] / static_cast<double>(s)) << "\t"
					<< (a / static_cast<double>(s))
					<< std::endl;
			}
			
		if ( s )
		{
			std::cerr << "[H]\tmin\t" << minq << "\t";
			if ( ( static_cast<uint64_t>(minq+33) < static_cast<uint64_t>(std::numeric_limits<char>::max()) && isprint(minq+33)) )
				std::cerr << static_cast<char>(minq+33);
			std::cerr << std::endl;
			std::cerr << "[H]\tmax\t" << maxq << "\t";
			if ( ( static_cast<uint64_t>(maxq+33) < static_cast<uint64_t>(std::numeric_limits<char>::max()) && isprint(maxq+33)) )
				std::cerr << static_cast<char>(maxq+33);
			std::cerr << std::endl;
		}
	}
	
	return EXIT_SUCCESS;
}

int bamvalidate(::libmaus2::util::ArgInfo const & arginfo)
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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print stats at the end of a successfull run" ) );
				V.push_back ( std::pair<std::string,std::string> ( "basequalhist=<["+::biobambam2::Licensing::formatNumber(getDefaultBaseQualHist())+"]>", "print base quality histogram at end of a successfull run" ) );
				V.push_back ( std::pair<std::string,std::string> ( "passthrough=<["+::biobambam2::Licensing::formatNumber(getDefaultPassThrough())+"]>", "write alignments to standard output (default: do not pass through)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory (passthrough=1, index=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0, passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0, passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus2::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus2::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ", passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1, passthrough=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset, passthrough=1 only)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

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
