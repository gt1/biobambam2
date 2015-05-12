/**
    bambam
    Copyright (C) 2009-2015 German Tischler
    Copyright (C) 2011-2015 Genome Research Limited

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

#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>

#include <libmaus2/bambam/parallel/BlockSortControl.hpp>
#include <libmaus2/bambam/parallel/BlockMergeControl.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBufferPosComparator.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBufferQueryNameComparator.hpp>

#include <libmaus2/digest/Digests.hpp>

#if defined(LIBMAUS2_HAVE_SHA2_ASSEMBLY)
#include <libmaus2/digest/DigestFactory_SHA2_ASM.hpp>
#endif

#if defined(LIBMAUS2_HAVE_SMMINTRIN_H) && defined(HAVE_SSE4)
#include <libmaus2/digest/DigestFactory_CRC32C_SSE42.hpp>
#endif

#include <libmaus2/digest/DigestFactoryContainer.hpp>

#include <libmaus2/parallel/NumCpus.hpp>

#include <libmaus2/util/ArgInfo.hpp>
#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultTempLevel() { return Z_BEST_SPEED; }
// static std::string getDefaultSortOrder() { return "coordinate"; }
static std::string getDefaultInputFormat() { return "bam"; }
static std::string getDefaultSeqChksumHash() { return "crc32prod"; }
static std::string getDefaultDigest() { return "md5"; }

template<typename heap_element_type>
static int 
	mergeBlocks
(
	libmaus2::parallel::SimpleThreadPool & STP,
	std::ostream & out,
	libmaus2::autoarray::AutoArray<char> & sheader,
	std::vector<libmaus2::bambam::parallel::GenericInputControlStreamInfo> const & BI,
	libmaus2::bitio::BitVector::unique_ptr_type const & Pdupvec,
	int level,
	uint64_t inputblocksize,
	uint64_t inputblocksperfile,
	uint64_t mergebuffersize,
	uint64_t mergebuffers,
	uint64_t complistsize,
	std::string const & hash,
	std::string const & headerchecksumstr,
	std::string const & digesttype,
	std::string const & indextmpfileprefix,
	std::string const & indexfilename,
	std::string const & digestfilename,
	libmaus2::bambam::parallel::BlockMergeControlTypeBase::block_merge_output_format_t const oformat,
	bool const bamindex,
	std::string const & sortordername
)
{	
	typedef libmaus2::digest::DigestInterface digest_interface_type;
	typedef digest_interface_type::unique_ptr_type digest_interface_pointer_type;
	digest_interface_pointer_type Pdigest(::libmaus2::digest::DigestFactoryContainer::construct(digesttype));

	bool const computerefidintervals = (oformat == libmaus2::bambam::parallel::BlockMergeControlTypeBase::output_format_cram) && (sortordername == "coordinate");
	
	libmaus2::bambam::parallel::BlockMergeControl<heap_element_type> BMC(
		STP, // libmaus2::parallel::SimpleThreadPool &
		out, // std::ostream & 
		sheader, // libmaus2::autoarray::AutoArray<char> const & 
		BI, // std::vector<libmaus2::bambam::parallel::GenericInputControlStreamInfo> const &
		Pdupvec.get(), // libmaus2::bitio::BitVector::unique_ptr_type const &
		level, // int
		inputblocksize, // uint64_t
		inputblocksperfile /* uint64_t, blocks per channel */,
		mergebuffersize /* uint64_t, merge buffer size */,
		mergebuffers /* uint64_t, number of merge buffers */,
		complistsize /* uint64_t, number of bgzf preload blocks */,
		hash, // std::string
		indextmpfileprefix,
		Pdigest.get(),
		oformat,
		bamindex,
		computerefidintervals
	);
	BMC.addPending();			
	BMC.waitWritingFinished();
	std::ostringstream finalheaderchecksumstr;
	BMC.printChecksumsForBamHeader(finalheaderchecksumstr);
	// std::cerr << finalheaderchecksumstr.str();
	if ( indexfilename.size() )
		BMC.writeBamIndex(indexfilename);
	std::string const digeststr = BMC.getFileDigest();
	if ( digestfilename.size() )
	{
		libmaus2::aio::PosixFdOutputStream PFOS(digestfilename);
		PFOS << digeststr << std::endl;
		PFOS.flush();
	}

	std::cerr << "[D]\t" << digesttype << "\t" << digeststr << std::endl;

	if ( finalheaderchecksumstr.str() != headerchecksumstr )
	{
		std::cerr << "[W] checksum mismatch between input and output" << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cerr << "[V] checksum ok" << std::endl;
		return EXIT_SUCCESS;
	}	
}

template<typename order_type, typename heap_element_type, bool create_dup_mark_info>
int bamsormadupTemplate(
	::libmaus2::util::ArgInfo const & arginfo,
	bool const bamindex,
	bool const fixmates,
	bool const dupmarksupport,
	std::string const & sortordername
)
{
	int returncode = EXIT_SUCCESS;

	libmaus2::timing::RealTimeClock progrtc; progrtc.start();
	// typedef libmaus2::bambam::parallel::FragmentAlignmentBufferPosComparator order_type;
	// typedef libmaus2::bambam::parallel::FragmentAlignmentBufferNameComparator order_type;
	
	libmaus2::timing::RealTimeClock rtc;
	
	rtc.start();
	uint64_t const numlogcpus = arginfo.getValue<int>("threads",libmaus2::parallel::NumCpus::getNumLogicalProcessors());
	libmaus2::aio::PosixFdInputStream in(STDIN_FILENO,256*1024);
	std::string const tmpfilebase = arginfo.getUnparsedValue("tmpfile",arginfo.getDefaultTmpFileName());
	int const templevel = arginfo.getValue<int>("templevel",getDefaultTempLevel());
	std::string const seqchksumhash = arginfo.getValue<std::string>("seqchksumhash",getDefaultSeqChksumHash());
	std::string const digest = arginfo.getValue<std::string>("digest",getDefaultDigest());
	std::string const indexfilename = arginfo.getUnparsedValue("indexfilename",std::string());
	std::string const digestfilename = arginfo.getUnparsedValue("digestfilename",std::string());

	std::string const sinputformat = arginfo.getUnparsedValue("inputformat",getDefaultInputFormat());
	libmaus2::bambam::parallel::BlockSortControlBase::block_sort_control_input_enum inform = libmaus2::bambam::parallel::BlockSortControlBase::block_sort_control_input_bam;
	
	if ( sinputformat == "bam" )
	{
		inform = libmaus2::bambam::parallel::BlockSortControlBase::block_sort_control_input_bam;
	}
	else if ( sinputformat == "sam" )
	{
		inform = libmaus2::bambam::parallel::BlockSortControlBase::block_sort_control_input_sam;			
	}
	else
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "Unknown input format " << sinputformat << std::endl;
		lme.finish();
		throw lme;				
	}
				
	libmaus2::parallel::SimpleThreadPool STP(numlogcpus);
	typename libmaus2::bambam::parallel::BlockSortControl<order_type,create_dup_mark_info>::unique_ptr_type VC(
		new libmaus2::bambam::parallel::BlockSortControl<order_type,create_dup_mark_info>(
			inform,STP,in,templevel,tmpfilebase,seqchksumhash,
			fixmates,
			dupmarksupport
		)
	);
	VC->enqueReadPackage();
	VC->waitDecodingFinished();
	// VC->printChecksums(std::cerr);
	std::ostringstream headerchecksumstr;
	VC->printChecksumsForBamHeader(headerchecksumstr);
	// std::cerr << headerchecksumstr.str();
	// VC->printSizes(std::cerr);
	// VC->printPackageFreeListSizes(std::cerr);
	#if defined(AUTOARRAY_TRACE)
	libmaus2::autoarray::autoArrayPrintTraces(std::cerr);
	#endif
	VC->freeBuffers();

	/**
	 * set up metrics stream
	 **/
	::libmaus2::aio::CheckedOutputStream::unique_ptr_type pM;
	std::ostream * pmetricstr = 0;
	
	if ( arginfo.hasArg("M") && (arginfo.getValue<std::string>("M","") != "") )
	{
		::libmaus2::aio::CheckedOutputStream::unique_ptr_type tpM(
                                new ::libmaus2::aio::CheckedOutputStream(arginfo.getValue<std::string>("M",std::string("M")))
                        );
		pM = UNIQUE_PTR_MOVE(tpM);
		pmetricstr = pM.get();
	}
	else
	{
		pmetricstr = & std::cerr;
	}

	std::ostream & metricsstr = *pmetricstr;

	if ( create_dup_mark_info )
	{
		VC->flushReadEndsLists();
		VC->mergeReadEndsLists(metricsstr,"bamsormadup");
	}
	
	metricsstr.flush();
	pM.reset();

	std::vector<libmaus2::bambam::parallel::GenericInputControlStreamInfo> const BI = VC->getBlockInfo();
	libmaus2::bitio::BitVector::unique_ptr_type Pdupvec(VC->releaseDupBitVector());
	libmaus2::bambam::BamHeader::unique_ptr_type Pheader(VC->getHeader());
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(libmaus2::bambam::BamHeaderUpdate::updateHeader(arginfo,*Pheader,"bamsormadup",PACKAGE_VERSION));
	uphead->changeSortOrder(sortordername /* "coordinate" */);
	uphead->text = uphead->filterOutChecksum(uphead->text);
	uphead->text += headerchecksumstr.str();

	libmaus2::bambam::parallel::BlockMergeControlTypeBase::block_merge_output_format_t oformat = libmaus2::bambam::parallel::BlockMergeControlTypeBase::output_format_bam;
		
	if ( arginfo.getUnparsedValue("outputformat","bam") == "sam" )
		oformat = libmaus2::bambam::parallel::BlockMergeControlTypeBase::output_format_sam;
	if ( arginfo.getUnparsedValue("outputformat","bam") == "cram" )
		oformat = libmaus2::bambam::parallel::BlockMergeControlTypeBase::output_format_cram;

	std::string const reference = arginfo.getUnparsedValue("reference",std::string());
	if ( oformat == libmaus2::bambam::parallel::BlockMergeControlTypeBase::output_format_cram )
	{
		try
		{
			uphead->checkSequenceChecksums(reference);
			
			if ( ! uphead->checkSequenceChecksumsCached(false /* throw */) )
			{
				char const * refcache = getenv("REF_CACHE");
				
				if ( (! refcache) || (!*refcache) )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "Sequence cache is missing sequences but REF_CACHE is not set" << std::endl;
					lme.finish();
					throw lme;	
				}
				
				// try to fill cache
				uphead->getSequenceURSet(true);
			}

			if ( ! uphead->checkSequenceChecksumsCached(true /* throw */) )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Sequence cache is missing sequences" << std::endl;
				lme.finish();
				throw lme;
			}
		}
		catch(...)
		{
			STP.terminate();
			STP.join();
			throw;
		}
	}

	std::ostringstream hostr;
	uphead->serialise(hostr);
	std::string const hostrstr = hostr.str();
	libmaus2::autoarray::AutoArray<char> sheader(hostrstr.size(),false);
	std::copy(hostrstr.begin(),hostrstr.end(),sheader.begin());		
	VC.reset();
				
	std::cerr << "[V] blocks generated in time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	
	rtc.start();
	uint64_t const inputblocksize = 1024*1024;
	uint64_t const inputblocksperfile = 8;
	uint64_t const mergebuffermemory = arginfo.getValueUnsignedNumeric("mergebuffermemory", 1024ull * 1024ull * 1024ull);
	uint64_t const complistsize = 32;
	int const level = arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION);


	uint64_t const mergebuffers = 
		(oformat == libmaus2::bambam::parallel::BlockMergeControlTypeBase::output_format_cram)
		?
		(2*STP.getNumThreads())
		:
		4
	;
	uint64_t const mergebuffersize = std::max(mergebuffermemory / mergebuffers,static_cast<uint64_t>(1ull));

	try
	{
		mergeBlocks<heap_element_type>(
			STP,std::cout,sheader,BI,Pdupvec,level,inputblocksize,inputblocksperfile,mergebuffersize,mergebuffers,complistsize,seqchksumhash,headerchecksumstr.str(),
			digest,
			tmpfilebase+"_index",
			indexfilename,
			digestfilename,
			oformat,
			bamindex,
			sortordername
		);
	}
	catch(...)
	{
		STP.terminate();
		STP.join();
		throw;
	}
	
	std::cerr << "[V] blocks merged in time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

	STP.terminate();
	STP.join();
	
	std::cerr << "[V] run time " << progrtc.formatTime(progrtc.getElapsedSeconds()) << " (" << progrtc.getElapsedSeconds() << " s)" << "\t" << libmaus2::util::MemUsage() << std::endl;

	return returncode;
}

int bamsormadup(::libmaus2::util::ArgInfo const & arginfo)
{
	std::string const so = arginfo.getUnparsedValue("SO","coordinate");
	
	if ( so == "queryname" )
		return bamsormadupTemplate<
				libmaus2::bambam::parallel::FragmentAlignmentBufferQueryNameComparator,
				libmaus2::bambam::parallel::GenericInputControlMergeHeapEntryQueryName,
				false /* create dup mark info */>(arginfo,false /* bam index */,false /* fix mates */,false /* dup mark support */,"queryname");
	else
		return bamsormadupTemplate<
				libmaus2::bambam::parallel::FragmentAlignmentBufferPosComparator,
				libmaus2::bambam::parallel::GenericInputControlMergeHeapEntryCoordinate,
				true /* create dup mark info */>(arginfo,true /* bam index */,true /* fix mates */,true /* dup mark support */,"coordinate");
}

int main(int argc, char * argv[])
{
	try
	{
		#if defined(LIBMAUS2_HAVE_SHA2_ASSEMBLY)
		libmaus2::digest::DigestFactoryContainer::addFactories(libmaus2::digest::DigestFactory_SHA2_ASM());
		#endif
		#if defined(LIBMAUS2_HAVE_SMMINTRIN_H) && defined(HAVE_SSE4)
		libmaus2::digest::DigestFactoryContainer::addFactories(libmaus2::digest::DigestFactory_CRC32C_SSE42());
		#endif
		
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
				V.push_back ( std::pair<std::string,std::string> ( "templevel=<["+::biobambam2::Licensing::formatNumber(getDefaultTempLevel())+"]>", "compression setting for temporary files (see level for options)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "threads=<["+::biobambam2::Licensing::formatNumber(libmaus2::parallel::NumCpus::getNumLogicalProcessors())+"]>", "number of threads" ) );

				// V.push_back ( std::pair<std::string,std::string> ( "SO=<["+getDefaultSortOrder()+"]>", "sorting order (coordinate or queryname)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("tmpfile=<[")+arginfo.getDefaultTmpFileName()+"]>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (sam,bam)") ) );
				V.push_back ( std::pair<std::string,std::string> ( "M=<filename>", "metrics file, stderr if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("seqchksumhash=<[")+getDefaultSeqChksumHash()+"]>", "seqchksum digest function: " + libmaus2::bambam::ChecksumsFactory::getSupportedHashVariantsList()) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("digest=<[")+getDefaultDigest()+"]>", "hash digest computed for output stream (md5, sha512)") );
				V.push_back ( std::pair<std::string,std::string> ( std::string("digestfilename=<[]>"), "name of file for storing hash digest computed for output stream (not stored by default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("indexfilename=<[]>"), "name of file for storing BAM index (not stored by default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("SO=<coordinate|queryname>"), "output sort order (coordinate by default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[bam]>"), std::string("output format (sam,bam,cram)") ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("reference=<[]>"), std::string("reference FastA for writing cram") ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return bamsormadup(arginfo);			
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

