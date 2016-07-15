/**
    biobambam2
    Copyright (C) 2009-2016 German Tischler
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
#include <libmaus2/aio/OutputStreamInstance.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/math/iabs.hpp>
#include <libmaus2/util/GrowingFreeList.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/bambam/BamHeaderUpdate.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/StrCmpNum.hpp>
#include <libmaus2/util/FreeList.hpp>
#include <libmaus2/util/SimpleQueue.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

#include <config.h>

static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

static int getDefaultVerbose()
{
	return 1;
}

static int getDefaultDisableValidation()
{
	return false;
}

static std::string getDefaultInputFormat()
{
	return "bam";
}

static int getDefaultNumThreads()
{
	return 1;
}


struct BamAlignmentFreeListDefaultTypeInfo
{
	typedef libmaus2::bambam::BamAlignment element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}
};

struct BamAlignmentFreeListDefaultAllocator
{
	typedef libmaus2::bambam::BamAlignment element_type;

	BamAlignmentFreeListDefaultAllocator()
	{
	}
	~BamAlignmentFreeListDefaultAllocator()
	{
	}

	element_type::shared_ptr_type operator()() const
	{
		element_type::shared_ptr_type ptr(new element_type);
		return ptr;
	}
};

void handleQueue(
	libmaus2::util::SimpleQueue < libmaus2::bambam::BamAlignment::shared_ptr_type > & Q,
	libmaus2::util::FreeList < libmaus2::bambam::BamAlignment, BamAlignmentFreeListDefaultAllocator, BamAlignmentFreeListDefaultTypeInfo > & FL,
	libmaus2::bambam::BamBlockWriterBase & wr,
	libmaus2::bambam::BamAuxFilterVector const & auxvec,
	uint64_t const numthreads
)
{
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < Q.size(); ++i )
	{
		libmaus2::bambam::BamAlignment::shared_ptr_type P = Q[i];

		if ( P->isMateUnmap() )
			P->filterOutAux(auxvec);
	}

	while ( ! Q.empty() )
	{
		libmaus2::bambam::BamAlignment::shared_ptr_type P = Q.pop_front();
		wr.writeAlignment(*P);
		FL.put(P);
	}
}

int bamfiltermc(libmaus2::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue("verbose",getDefaultVerbose());

	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();
	::libmaus2::bambam::BamAlignment & algn = dec.getAlignment();
	std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpfile",arginfo.getDefaultTmpFileName());
	uint64_t const numthreads = arginfo.getValueUnsignedNumeric<uint64_t>("numthreads",getDefaultNumThreads());

	/*
	 * start index/md5 callbacks
	 */
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

	::libmaus2::bambam::BamHeader::unique_ptr_type genuphead(
		libmaus2::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bamfiltermc",std::string(PACKAGE_VERSION))
	);
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(libmaus2::bambam::BamBlockWriterBaseFactory::construct(*genuphead,arginfo,Pcbs));
	libmaus2::bambam::BamBlockWriterBase & wr = *Pwriter;

	// freelist size
	uint64_t const flsize = 16*1024;
	libmaus2::util::FreeList < libmaus2::bambam::BamAlignment, BamAlignmentFreeListDefaultAllocator, BamAlignmentFreeListDefaultTypeInfo > FL(flsize);
	libmaus2::util::SimpleQueue < libmaus2::bambam::BamAlignment::shared_ptr_type > Q;

	libmaus2::bambam::BamAuxFilterVector auxvec;
	auxvec.set('M','C');

	uint64_t alcnt = 0;
	while ( dec.readAlignment() )
	{
		if ( FL.empty() )
			handleQueue(Q,FL,wr,auxvec,numthreads);
		assert ( ! FL.empty() );

		libmaus2::bambam::BamAlignment::shared_ptr_type P = FL.get();
		P->swap(algn);

		Q.push_back(P);

		if ( verbose && ((++alcnt % (1024*1024)) == 0) )
			std::cerr << "[V] " << alcnt << std::endl;
	}

	handleQueue(Q,FL,wr,auxvec,numthreads);

	// reset BAM writer
	Pwriter.reset();

	if ( Pmd5cb )
		Pmd5cb->saveDigestAsFile(md5filename);
	if ( Pindex )
		Pindex->flush(std::string(indexfilename));

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

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

				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "numthreads=<[1]>", "number of worker threads (default: 1)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return bamfiltermc(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
