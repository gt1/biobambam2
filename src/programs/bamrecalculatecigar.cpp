/**
    bambam
    Copyright (C) 2009-2016 German Tischler
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
#include <libmaus2/bambam/BamCat.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>

#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

::libmaus2::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus2::util::ArgInfo const & arginfo,
	::libmaus2::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamrecalculatecigar", // ID
		"bamrecalculatecigar", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN
	);
	// construct new header
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus2::bambam::BamHeader(upheadtext));

	return UNIQUE_PTR_MOVE(uphead);
}

#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>

int bamrecalculatecigar(libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	// input decoder wrapper
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);

	libmaus2::bambam::BamAlignmentDecoder & bamdec = decwrapper->getDecoder();
	libmaus2::bambam::BamAlignment & algn = bamdec.getAlignment();
	libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,header));

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

	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus2::bambam::BamBlockWriterBaseFactory::construct(*uphead, arginfo, Pcbs)
	);

	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigopin;
	libmaus2::autoarray::AutoArray<char> readdata;
	libmaus2::bambam::BamAlignment::D_array_type T;

	if ( ! arginfo.hasArg("reference") )
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "reference key is missing." << std::endl;
		se.finish();
		throw se;
	}

	std::string const reference = arginfo.getUnparsedValue("reference","");

	if ( ! libmaus2::util::GetFileSize::fileExists(reference) )
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "file " << reference << " does not exist." << std::endl;
		se.finish();
		throw se;
	}

	libmaus2::fastx::FastAIndex::unique_ptr_type FAindex(libmaus2::fastx::FastAIndex::load(reference + ".fai"));
	libmaus2::aio::InputStreamInstance FAISI(reference);

	uint64_t c = 0;
	libmaus2::autoarray::AutoArray<char> ref;
	int64_t refloaded = -1;

	while ( bamdec.readAlignment() )
	{
		if ( algn.isMapped() )
		{
			assert ( algn.getRefID() >= 0 );
			if ( algn.getRefID() != refloaded )
			{
				if ( algn.getRefID() < refloaded )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "bamrecalculatecigar: file is not sorted by coordinate" << std::endl;
					lme.finish();
					throw lme;
				}

				ref = FAindex->readSequence(FAISI,algn.getRefID());
				refloaded = algn.getRefID();
			}

			assert (
				algn.getPos() + algn.getReferenceLength() <= ref.size()
			);

			uint64_t const numcig = libmaus2::bambam::BamAlignmentDecoderBase::recalculateCigar(
				algn.D.begin(),
				ref.begin() + algn.getPos(),
				cigopin,
				readdata
			);
			algn.replaceCigarString(cigopin,numcig,T);
		}

		writer->writeAlignment(algn);

		if ( ((++c) & ((1ull<<20)-1)) == 0 && verbose )
			std::cerr << "[V] " << c << std::endl;
	}

	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

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

				V.push_back ( std::pair<std::string,std::string> ( "I=<filename>", "input file, can be set multiple times" ) );
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

		return bamrecalculatecigar(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
