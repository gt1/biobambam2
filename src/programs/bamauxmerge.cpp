/**
    biobambam2
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

#include <libmaus2/aio/OutputStreamInstance.hpp>

#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamAlignmentNameComparator.hpp>
#include <libmaus2/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamEntryContainer.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/GetObject.hpp>
#include <libmaus2/util/PutObject.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamPeeker.hpp>

#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

int bamauxmerge(::libmaus2::util::ArgInfo const & arginfo)
{
	::libmaus2::util::TempFileRemovalContainer::setup();

	libmaus2::util::ArgInfo arg_a = arginfo;
	arg_a.replaceKey("I",arginfo.restargs.at(0));

	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper_a(libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arg_a,false /* put rank */));
	::libmaus2::bambam::BamAlignmentDecoder * ppdec_a = &(decwrapper_a->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec_a = *ppdec_a;
	::libmaus2::bambam::BamHeader const & header_a = dec_a.getHeader();
	//::libmaus2::bambam::BamAlignment const & algn_a = dec_a.getAlignment();
	::libmaus2::bambam::BamPeeker peek_a(dec_a);

	libmaus2::util::ArgInfo arg_b = arginfo;
	arg_b.replaceKey("I",arginfo.restargs.at(1));

	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper_b(libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arg_b,false /* put rank */));
	::libmaus2::bambam::BamAlignmentDecoder * ppdec_b = &(decwrapper_b->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec_b = *ppdec_b;
	::libmaus2::bambam::BamHeader const & header_b = dec_b.getHeader();
	// ::libmaus2::bambam::BamAlignment const & algn_b = dec_b.getAlignment();
	::libmaus2::bambam::BamPeeker peek_b(dec_b);

	std::vector<libmaus2::bambam::HeaderLine> VSQ = libmaus2::bambam::HeaderLine::extractLinesByType(header_b.text,"SQ");
	std::ostringstream outhdrostr;
	outhdrostr << header_a.text;

	for ( uint64_t i = 0; i < VSQ.size(); ++i )
		outhdrostr << VSQ[i].line << "\n";

	std::string const headertext(outhdrostr.str());

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamauxmerge", // ID
		"bamauxmerge", // PN
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

	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus2::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, Pcbs)
	);

	libmaus2::bambam::BamAuxSortingBuffer sortbuf;

	libmaus2::bambam::BamAlignment algn_a;
	libmaus2::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > auxlist_a;
	libmaus2::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > auxlist_b;
	libmaus2::bambam::BamAuxFilterVector auxvec;

	while ( peek_a.getNext(algn_a) )
	{
		std::string const name_a = algn_a.getName();
		algn_a.sortAux(sortbuf);
		uint64_t const numaux_a = algn_a.enumerateAuxTags(auxlist_a);

		libmaus2::bambam::BamAlignment algn_b;
		while ( peek_b.peekNext(algn_b) && algn_b.getName() == name_a )
		{
			peek_b.getNext(algn_b);
			algn_b.sortAux(sortbuf);
			uint64_t const numaux_b = algn_b.enumerateAuxTags(auxlist_b);

			for ( uint64_t i = 0; i < numaux_a; ++i )
				auxvec.set(auxlist_a[i].first,auxlist_a[i].second);
			for ( uint64_t i = 0; i < numaux_b; ++i )
				auxvec.clear(auxlist_b[i].first,auxlist_b[i].second);

			algn_b.copyAuxTags(algn_a,auxvec);

			writer->writeAlignment(algn_b);
		}

		for ( uint64_t i = 0; i < numaux_a; ++i )
			auxvec.clear(auxlist_a[i].first,auxlist_a[i].second);
	}

	if ( peek_b.peekNext(algn_a) )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] mapped file contains read with name " << algn_a.getName() << " which is not in unmapped file" << std::endl;
		lme.finish();
		throw lme;
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
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;

				return EXIT_SUCCESS;
			}

		return bamauxmerge(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
