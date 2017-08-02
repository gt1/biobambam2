/**
    bambam
    Copyright (C) 2009-2017 German Tischler
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

#include <libmaus2/bambam/BamPeeker.hpp>
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/StrCmpNum.hpp>

#include <libmaus2/util/ArgInfo.hpp>

#include <biobambam2/Licensing.hpp>
#include <biobambam2/AttachRank.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static std::string getDefaultInputFormat() { return "bam"; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

int bamfillquery(::libmaus2::util::ArgInfo const & arginfo)
{
	if ( arginfo.getNumRestArgs() < 1 )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "[E] required first argument (FastA containing query sequences) is missing" << std::endl;
		se.finish();
		throw se;
	}

	// int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());

	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,
			false, // do not put rank
			0, /* copy stream */
			std::cin /* standard input */
		)
	);
	::libmaus2::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();

	libmaus2::bambam::BamPeeker BP(dec);
	::libmaus2::bambam::BamHeader const & header = BP.dec.getHeader();
	libmaus2::fastx::FastaPeeker FP(arginfo.getUnparsedRestArg(0));

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamfillquery", // ID
		"bamfillquery", // PN
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
		if ( libmaus2::bambam::BamBlockWriterBaseFactory::getMD5FileName(arginfo) != std::string() )
			md5filename = libmaus2::bambam::BamBlockWriterBaseFactory::getMD5FileName(arginfo);
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
		if ( libmaus2::bambam::BamBlockWriterBaseFactory::getIndexFileName(arginfo) != std::string() )
			indexfilename = libmaus2::bambam::BamBlockWriterBaseFactory::getIndexFileName(arginfo);
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

	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Uout ( libmaus2::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, Pcbs) );
	libmaus2::bambam::BamBlockWriterBase & Pout = *Uout;

	libmaus2::bambam::BamAlignment algn;
	uint64_t c = 0;

	libmaus2::fastx::FastAReader::pattern_type pat;
	std::string prevname;
	bool prevnamevalid = false;
	std::string prevbamname;
	bool prevbamnamevalid = false;

	while ( FP.peekNext(pat) )
	{
		FP.getNext(pat);
		std::string const name = pat.getShortStringId();

		if ( prevnamevalid )
		{
			int const r = libmaus2::bambam::StrCmpNum::strcmpnum(prevname.c_str(),name.c_str());

			if ( r >= 0 )
			{
				::libmaus2::exception::LibMausException se;
				se.getStream() << "[E] FastA file is not sorted by query name" << std::endl;
				se.getStream() << "[E] " << prevname << std::endl;
				se.getStream() << "[E] " << name << std::endl;
				se.finish();
				throw se;
			}
		}

		if ( verbose )
			std::cerr << name << std::endl;

		// check order in alignment file
		if ( prevbamnamevalid )
		{
			if ( BP.peekNext(algn) )
			{
				if ( libmaus2::bambam::StrCmpNum::strcmpnum(algn.getName(),prevbamname.c_str()) < 0 )
				{
					::libmaus2::exception::LibMausException se;
					se.getStream() << "[E] Alignment stream is not sorted by query name" << std::endl;
					se.getStream() << "[E] " << algn.getName() << std::endl;
					se.finish();
					throw se;
				}
			}
		}

		if ( BP.peekNext(algn) && libmaus2::bambam::StrCmpNum::strcmpnum(algn.getName(),name.c_str()) < 0 )
		{
			::libmaus2::exception::LibMausException se;
			se.getStream() << "[E] Alignment stream contains query name not cotained in FastA input or files are not properly sorted" << std::endl;
			se.getStream() << "[E] " << algn.getName() << std::endl;
			se.finish();
			throw se;
		}

		while ( BP.peekNext(algn) && algn.getName() == name )
		{
			BP.getNext(algn);

			std::string seq;
			if ( algn.isReverse() )
			{
				seq = libmaus2::fastx::reverseComplementUnmapped(libmaus2::fastx::remapString(libmaus2::fastx::mapString(pat.spattern)));
			}
			else
			{
				seq = libmaus2::fastx::remapString(libmaus2::fastx::mapString(pat.spattern));
			}

			if ( algn.isMapped() )
			{
				uint64_t const fc = algn.getFrontHardClipping();
				uint64_t const bc = algn.getBackHardClipping();
				uint64_t const l = seq.size() - (fc+bc);

				assert ( fc + bc <= seq.size() );

				seq = seq.substr(fc,l);
			}

			algn.replaceSequence(
				seq,
				std::string(seq.size(),'H')
			);

			Pout.writeAlignment(algn);

			prevbamnamevalid = true;
			prevbamname = algn.getName();

			++c;
		}

		prevname = name;
		prevnamevalid = true;
	}

	Uout.reset();

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
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus2::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus2::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;

				std::cerr << "The keep and remove keys are mutually exclusive. Tags are given by their two character ids. Multiple ids are separated by commas." << std::endl;

				return EXIT_SUCCESS;
			}

		return bamfillquery(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
