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

#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamAlignmentHeapComparator.hpp>
#include <libmaus/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamHeader.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/lz/BgzfInflateBase.hpp>
#include <libmaus/lz/SimpleCompressedConcatInputStream.hpp>
#include <libmaus/lz/SimpleCompressedConcatInputStreamFragment.hpp>
#include <libmaus/lz/SimpleCompressedInputBlockConcatBlock.hpp>
#include <libmaus/lz/SimpleCompressedInputBlockConcat.hpp>
#include <libmaus/lz/SimpleCompressedOutputStream.hpp>
#include <libmaus/lz/SimpleCompressedStreamInterval.hpp>
#include <libmaus/lz/SimpleCompressedStreamNamedInterval.hpp>
#include <libmaus/lz/SnappyCompressorObjectFactory.hpp>
#include <libmaus/lz/SnappyDecompressorObjectFactory.hpp>
#include <libmaus/lz/ZlibCompressorObjectFactory.hpp>
#include <libmaus/lz/ZlibDecompressorObjectFactory.hpp>
#include <libmaus/parallel/LockedBool.hpp>
#include <libmaus/parallel/LockedQueue.hpp>
#include <libmaus/parallel/NumCpus.hpp>
#include <libmaus/parallel/SimpleThreadPool.hpp>
#include <libmaus/parallel/SimpleThreadPoolWorkPackageFreeList.hpp>
#include <libmaus/sorting/ParallelStableSort.hpp>
#include <libmaus/util/GetObject.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <config.h>

int getDefaultLevel()
{
	return Z_DEFAULT_COMPRESSION;
}

std::string getDefaultSortOrder()
{
	return "coordinate";
}

bool getDefaultVerbose()
{
	return true;
}

uint64_t getDefaultMem()
{
	return 2ull * 1024ull * 1024ull * 1024ull;
}

uint64_t getDefaultProcessBuffers()
{
	return 4;
}

struct BamThreadPoolDecodeBamParseQueueInfo
{
	uint64_t packageid;
	std::pair<uint64_t,uint64_t> blockmeta;
	uint64_t baseid;
	uint64_t blockid;
	
	BamThreadPoolDecodeBamParseQueueInfo()
	: packageid(0), blockmeta(0,0), baseid(0), blockid(0)
	{
	
	}
	BamThreadPoolDecodeBamParseQueueInfo(
		uint64_t rpackageid,
		std::pair<uint64_t,uint64_t> rblockmeta,
		uint64_t rbaseid,
		uint64_t rblockid
	)
	: packageid(rpackageid), blockmeta(rblockmeta), baseid(rbaseid), blockid(rblockid)
	{
	
	}
	
	bool operator<(BamThreadPoolDecodeBamParseQueueInfo const & o) const
	{
		return blockid > o.blockid;
	}
};

struct BamThreadPoolDecodeBamProcessQueueInfo
{
	uint64_t packageid;
	uint64_t baseid;
	uint64_t blockid;
	
	BamThreadPoolDecodeBamProcessQueueInfo()
	: packageid(0), baseid(0), blockid(0)
	{
	
	}
	BamThreadPoolDecodeBamProcessQueueInfo(
		uint64_t rpackageid,
		uint64_t rbaseid,
		uint64_t rblockid
	)
	: packageid(rpackageid), baseid(rbaseid), blockid(rblockid)
	{
	
	}
	
	bool operator<(BamThreadPoolDecodeBamProcessQueueInfo const & o) const
	{
		return blockid > o.blockid;
	}
};

struct BamProcessBuffer
{
	typedef BamProcessBuffer this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus::util::shared_ptr<this_type>::type shared_ptr_type;
	
	typedef uint64_t pointer_type;

	libmaus::autoarray::AutoArray<pointer_type> B8;
	uint64_t const bspace;

	uint8_t * const ca;
	uint8_t * cc;
	
	pointer_type * const pa;
	pointer_type * pc;
	
	int64_t bufferid;
	
	void reset()
	{
		cc = ca;
		pc = pa;
	}
	
	BamProcessBuffer(uint64_t const bytesize)
	: B8(
		(bytesize + sizeof(pointer_type) - 1)/sizeof(pointer_type),
		false
	  ),
	  bspace(sizeof(pointer_type) * B8.size()),
	  ca(reinterpret_cast<uint8_t *>(B8.begin())), cc(ca),
	  pa(B8.end()), pc(pa),
	  bufferid(-1)
	{
	}
	
	uint64_t getFill() const
	{
		return pa-pc;
	}
	
	bool put(uint8_t const * p, uint64_t const s)
	{
		uint64_t const ptrmult = 2;
		uint64_t const spaceused_data = cc-ca;
		uint64_t const spaceused_ptr = (pa-pc)*sizeof(pointer_type)*ptrmult;
		uint64_t const spaceused_len = sizeof(uint32_t);
		uint64_t const spaceused = (spaceused_data+spaceused_ptr+spaceused_len);
		uint64_t const spaceav = bspace - spaceused;
		
		if ( spaceav >= s + sizeof(uint32_t) + ptrmult*sizeof(pointer_type) )
		{
			*(--pc) = cc-ca;
			for ( unsigned int i = 0; i < sizeof(uint32_t); ++i )
				*(cc++) = (s >> (i*8)) & 0xFF;
			std::copy(p,p+s,cc);
			cc += s;
			return true;
		}
		else
		{
			return false;
		}
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeContextBase;

struct BamThreadPoolDecodeContextBaseConstantsBase
{
	enum bamthreadpooldecodecontextbase_dispatcher_ids
	{
		bamthreadpooldecodecontextbase_dispatcher_id_read = 0,
		bamthreadpooldecodecontextbase_dispatcher_id_decompress = 1,
		bamthreadpooldecodecontextbase_dispatcher_id_bamparse = 2,
		bamthreadpooldecodecontextbase_dispatcher_id_bamprocess = 3,
		bamthreadpooldecodecontextbase_dispatcher_id_bamsort = 4,
		bamthreadpooldecodecontextbase_dispatcher_id_bamwrite = 5
	};
	
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_read = 10;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_decompress = 10;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_bamparse = 10;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_bamprocess = 10;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_bamsort = 0;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_bamwrite = 0;
};

template<typename _order_type>
struct BamThreadPoolDecodeReadPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolDecodeReadPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolDecodeContextBase<order_type> * contextbase;

	BamThreadPoolDecodeReadPackage() : contextbase(0) {}
	BamThreadPoolDecodeReadPackage(
		uint64_t const rpackageid,
		BamThreadPoolDecodeContextBase<order_type> * rcontextbase
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_read,
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_read,
		rpackageid
	), contextbase(rcontextbase)
	{
	
	}

	virtual char const * getPackageName() const
	{
		return "BamThreadPoolDecodeReadPackage";
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeDecompressPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolDecodeDecompressPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolDecodeContextBase<order_type> * contextbase;

	std::pair<uint64_t,uint64_t> blockmeta; // block size compressed and uncompressed
	uint64_t baseid;
	uint64_t blockid;

	BamThreadPoolDecodeDecompressPackage() : contextbase(0) {}
	BamThreadPoolDecodeDecompressPackage(
		uint64_t const rpackageid,
		BamThreadPoolDecodeContextBase<order_type> * rcontextbase,
		std::pair<uint64_t,uint64_t> rblockmeta,
		uint64_t rbaseid,
		uint64_t rblockid
		
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_decompress,
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_decompress,
		rpackageid
	), contextbase(rcontextbase), blockmeta(rblockmeta), baseid(rbaseid), blockid(rblockid)
	{
	
	}

	virtual char const * getPackageName() const
	{
		return "BamThreadPoolDecodeDecompressPackage";
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamParsePackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolDecodeBamParsePackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolDecodeContextBase<order_type> * contextbase;

	std::pair<uint64_t,uint64_t> blockmeta; // block size compressed and uncompressed
	uint64_t baseid;
	uint64_t blockid;

	BamThreadPoolDecodeBamParsePackage() : contextbase(0) {}
	BamThreadPoolDecodeBamParsePackage(
		uint64_t const rpackageid,
		BamThreadPoolDecodeContextBase<order_type> * rcontextbase,
		std::pair<uint64_t,uint64_t> rblockmeta,
		uint64_t rbaseid,
		uint64_t rblockid		
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_bamparse,
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_bamparse,
		rpackageid
	), contextbase(rcontextbase), blockmeta(rblockmeta), baseid(rbaseid), blockid(rblockid)
	{
	
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolDecodeBamParsePackage";
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamProcessPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolDecodeBamProcessPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolDecodeContextBase<order_type> * contextbase;
	uint64_t baseid;
	uint64_t blockid;

	BamThreadPoolDecodeBamProcessPackage() : contextbase(0) {}
	BamThreadPoolDecodeBamProcessPackage(
		uint64_t const rpackageid,
		BamThreadPoolDecodeContextBase<order_type> * rcontextbase,
		uint64_t rbaseid,
		uint64_t rblockid		
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_bamprocess,
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_bamprocess,
		rpackageid
	), contextbase(rcontextbase), baseid(rbaseid), blockid(rblockid)
	{
	
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolDecodeBamProcessPackage";
	}
};


template<typename _order_type>
struct BamSortInfo
{
	typedef _order_type order_type;
	typedef BamSortInfo this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;
	
	typedef libmaus::sorting::ParallelStableSort::ParallelSortControl<
		BamProcessBuffer::pointer_type *,order_type
	> sort_type;
	
	BamProcessBuffer * processBuffer;
	order_type BAPC;
	typename sort_type::unique_ptr_type sortControl;
	
	BamSortInfo(
		BamProcessBuffer * rprocessBuffer, uint64_t const rnumthreads
	) : processBuffer(rprocessBuffer), BAPC(processBuffer->ca),
	    sortControl(new sort_type(
	    	processBuffer->pc, // current
	    	processBuffer->pa, // end
	    	processBuffer->pc-(processBuffer->pa-processBuffer->pc),
	    	processBuffer->pc,
	    	BAPC,
	    	rnumthreads,
	    	true /* copy back */
	    ))
	{
		assert ( sortControl->context.ae-sortControl->context.aa == sortControl->context.be-sortControl->context.ba );

		#if 0
		for ( uint64_t * pc = processBuffer->pc; pc != processBuffer->pa; ++pc )
		{
			uint8_t const * cp = processBuffer->ca + *pc + 4;
			std::cerr << libmaus::bambam::BamAlignmentDecoderBase::getReadName(cp) << std::endl;
		}
		#endif
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamSortPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolDecodeBamSortPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	enum sort_package_type {
		sort_package_type_base,
		sort_package_type_merge_level,
		sort_package_type_merge_sublevel,
	};

	BamThreadPoolDecodeContextBase<order_type> * contextbase;
	uint64_t baseid;
	uint64_t blockid;
	typedef BamSortInfo<order_type> bam_sort_info_type;
	typename bam_sort_info_type::shared_ptr_type sortinfo;
	sort_package_type packagetype;	
	
	uint64_t sort_base_id;
	uint64_t sort_merge_id;
	uint64_t sort_submerge_id;

	BamThreadPoolDecodeBamSortPackage() : contextbase(0), packagetype(sort_package_type_base) {}
	BamThreadPoolDecodeBamSortPackage(
		uint64_t const rpackageid,
		BamThreadPoolDecodeContextBase<order_type> * rcontextbase,
		uint64_t rbaseid,
		uint64_t rblockid,
		typename bam_sort_info_type::shared_ptr_type rsortinfo,
		sort_package_type rpackagetype,
		uint64_t rsort_base_id,
		uint64_t rsort_merge_id,
		uint64_t rsort_submerge_id
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_bamsort,
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_bamsort,
		rpackageid
	), contextbase(rcontextbase), baseid(rbaseid), blockid(rblockid), sortinfo(rsortinfo), packagetype(rpackagetype),
	   sort_base_id(rsort_base_id), sort_merge_id(rsort_merge_id), sort_submerge_id(rsort_submerge_id)
	{
		assert ( 
			packagetype == sort_package_type_base ||
			packagetype == sort_package_type_merge_level ||
			packagetype == sort_package_type_merge_sublevel
		);
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolDecodeBamSortPackage";
	}
};

struct BamBlockWriteInfo
{
	typedef BamBlockWriteInfo this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus::util::shared_ptr<this_type>::type shared_ptr_type;
	
	BamProcessBuffer * processBuffer;
	std::vector < std::pair<BamProcessBuffer::pointer_type const *, BamProcessBuffer::pointer_type const *> > packets;
	libmaus::parallel::SynchronousCounter<uint64_t> packetsWritten;
	
	BamBlockWriteInfo() : processBuffer(0), packetsWritten(0) {}
	BamBlockWriteInfo(
		BamProcessBuffer * rprocessBuffer,
		std::vector < std::pair<BamProcessBuffer::pointer_type const *, BamProcessBuffer::pointer_type const *> > const & rpackets
	) : processBuffer(rprocessBuffer), packets(rpackets), packetsWritten(0) {}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamWritePackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolDecodeBamWritePackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolDecodeContextBase<order_type> * contextbase;
	uint64_t baseid;
	uint64_t blockid;
	BamBlockWriteInfo::shared_ptr_type writeinfo;
	uint64_t write_block_id;

	BamThreadPoolDecodeBamWritePackage() : contextbase(0) {}
	BamThreadPoolDecodeBamWritePackage(
		uint64_t const rpackageid,
		BamThreadPoolDecodeContextBase<order_type> * rcontextbase,
		uint64_t rbaseid,
		uint64_t rblockid,
		BamBlockWriteInfo::shared_ptr_type rwriteinfo,
		uint64_t rwrite_block_id
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_bamwrite,
		BamThreadPoolDecodeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_bamwrite,
		rpackageid
	), contextbase(rcontextbase), baseid(rbaseid), blockid(rblockid), writeinfo(rwriteinfo),
	   write_block_id(rwrite_block_id)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolDecodeBamWritePackage";
	}
	
	bool operator<(BamThreadPoolDecodeBamWritePackage const & o) const
	{
		return baseid > o.baseid;
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamWritePackageHeapComparator
{
	typedef _order_type order_type;
	
	bool operator()(BamThreadPoolDecodeBamWritePackage<order_type> const * A, BamThreadPoolDecodeBamWritePackage<order_type> const * B)
	{
		//return A->baseid > B->baseid;
		return A->blockid > B->blockid;
	}
};

struct MergeInfo
{
	std::vector< std::vector<uint64_t> > tmpfileblockcnts;
	std::vector< std::string > tmpfilenames;
	std::vector< std::vector< libmaus::lz::SimpleCompressedStreamInterval > > tmpfileblocks;
	std::vector< std::vector< libmaus::lz::SimpleCompressedStreamNamedInterval > > tmpfilenamedblocks;
	std::vector< uint64_t > tmpfileblockcntsums;
	libmaus::bambam::BamHeader::shared_ptr_type sheader;
	
	bool empty() const
	{
		return
			std::accumulate(
				tmpfileblockcntsums.begin(),
				tmpfileblockcntsums.end(),
				0ull
			) == 0;
	}
	
	MergeInfo() {}
	MergeInfo(
		std::vector< std::vector<uint64_t> > rtmpfileblockcnts,
		std::vector< std::string > rtmpfilenames,
		std::vector< std::vector< libmaus::lz::SimpleCompressedStreamInterval > > rtmpfileblocks,
		std::vector< std::vector< libmaus::lz::SimpleCompressedStreamNamedInterval > > rtmpfilenamedblocks,
		std::vector< uint64_t > rtmpfileblockcntsums,
		libmaus::bambam::BamHeader const & header
	)
	: tmpfileblockcnts(rtmpfileblockcnts),
	  tmpfilenames(rtmpfilenames),
	  tmpfileblocks(rtmpfileblocks),
	  tmpfilenamedblocks(rtmpfilenamedblocks),
	  tmpfileblockcntsums(rtmpfileblockcntsums),
	  sheader(header.sclone())
	{
	
	}
};

#include <libmaus/aio/IsKnownLocalFileSystem.hpp>

#if defined(__linux__)
#include <libmaus/aio/LinuxStreamingPosixFdOutputStream.hpp>
#include <libmaus/aio/PosixFdOutputStreamBuffer.hpp>
#endif


template<typename _order_type>
struct BamThreadPoolDecodeContextBase : public BamThreadPoolDecodeContextBaseConstantsBase
{
	typedef _order_type order_type;

	static uint64_t const alperpackalign = 1024;

	libmaus::parallel::PosixSpinLock cerrlock;

	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolDecodeReadPackage<order_type> > readFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolDecodeDecompressPackage<order_type> > decompressFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolDecodeBamParsePackage<order_type> > bamParseFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolDecodeBamProcessPackage<order_type> > bamProcessFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolDecodeBamSortPackage<order_type> > bamSortFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolDecodeBamWritePackage<order_type> > bamWriteFreeList;
	
	libmaus::parallel::PosixSpinLock inputLock;
	libmaus::timing::RealTimeClock inputRtc;
	std::istream & in;
	libmaus::parallel::SynchronousCounter<uint64_t> readCnt;
	libmaus::parallel::SynchronousCounter<uint64_t> readCompCnt;
	libmaus::parallel::LockedBool readComplete;
	libmaus::parallel::PosixSemaphore readSem;
	libmaus::autoarray::AutoArray< libmaus::lz::BgzfInflateBase::unique_ptr_type > inflateBases;
	libmaus::autoarray::AutoArray< char > inflateDecompressSpace;
	libmaus::parallel::SynchronousQueue<uint64_t> inflateBasesFreeList;
	libmaus::parallel::SynchronousCounter<uint64_t> nextInputBlockId;
	
	libmaus::parallel::PosixSpinLock decompressLock;
	uint64_t decompressCnt;
	libmaus::parallel::LockedBool decompressComplete;
	
	libmaus::parallel::PosixSpinLock bamParseQueueLock;
	std::priority_queue<BamThreadPoolDecodeBamParseQueueInfo> bamparseQueue;
	std::priority_queue<BamThreadPoolDecodeBamParseQueueInfo> bamparseStall;
	uint64_t bamParseNextBlock;

	libmaus::parallel::PosixSpinLock bamParseLock;
	uint64_t bamParseCnt;
	libmaus::parallel::LockedBool bamParseComplete;

	libmaus::parallel::LockedBool haveheader;
	::libmaus::bambam::BamHeader::BamHeaderParserState bamheaderparsestate;
	libmaus::bambam::BamHeader header;
	
	enum bam_parser_state_type {
		bam_parser_state_read_blocklength,
		bam_parser_state_read_blockdata
	};
	
	bam_parser_state_type bamParserState;
	unsigned int bamBlockSizeRead;
	uint32_t bamBlockSizeValue;
	uint32_t bamBlockDataRead;
	libmaus::bambam::BamAlignment bamGatherBuffer;
	
	libmaus::autoarray::AutoArray<BamProcessBuffer::unique_ptr_type> processBuffers;
	libmaus::parallel::LockedQueue<uint64_t> processBuffersFreeList;
	uint64_t nextProcessBufferIdIn;
	libmaus::parallel::PosixSpinLock nextProcessBufferIdOutLock;
	uint64_t nextProcessBufferIdOut;

	libmaus::parallel::LockedBool bamProcessComplete;

	libmaus::parallel::PosixSpinLock buffersSortedLock;
	uint64_t buffersSorted;
	libmaus::parallel::LockedBool bamSortComplete;

	std::vector<std::string> tmpfilenames;
	libmaus::autoarray::AutoArray<libmaus::aio::CheckedOutputStream::unique_ptr_type> tmpfilesCOS;
	#if defined(__linux__)
	libmaus::autoarray::AutoArray<libmaus::aio::LinuxStreamingPosixFdOutputStream::unique_ptr_type> tmpfilesLSPFDOSB;
	#endif
	libmaus::lz::CompressorObjectFactory & compressorFactory;
	libmaus::autoarray::AutoArray<libmaus::lz::SimpleCompressedOutputStream<std::ostream>::unique_ptr_type> compressedTmpFiles;
	libmaus::parallel::PosixSpinLock tmpfileblockslock;
	std::vector< std::vector< libmaus::lz::SimpleCompressedStreamInterval > > tmpfileblocks;
	std::vector< std::vector<uint64_t> > tmpfileblockcnts;
	std::vector< uint64_t > tmpfileblockcntsums;
	libmaus::parallel::PosixSpinLock tmpfileblockcntsumsacclock;
	uint64_t tmpfileblockcntsumsacc;

	MergeInfo getMergeInfo() const
	{
		std::vector< std::vector< libmaus::lz::SimpleCompressedStreamNamedInterval > > tmpfilenamedblocks;
		
		for ( uint64_t i = 0; i < tmpfileblocks.size(); ++i )
		{
			std::vector< libmaus::lz::SimpleCompressedStreamNamedInterval > subblocks(
				tmpfileblocks[i].size()
			);
			for ( uint64_t j = 0; j < tmpfileblocks[i].size(); ++j )
			{
				subblocks[j] = libmaus::lz::SimpleCompressedStreamNamedInterval(
					tmpfileblocks[i][j].start,
					tmpfileblocks[i][j].end,
					tmpfilenames[j]
				);
			}
			
			tmpfilenamedblocks.push_back(subblocks);
		}
			
		return MergeInfo(tmpfileblockcnts,tmpfilenames,tmpfileblocks,tmpfilenamedblocks,tmpfileblockcntsums,header);
	}

	libmaus::parallel::PosixSpinLock writesPendingLock;
	std::vector< uint64_t > writesNext;
	std::vector< 
		std::priority_queue<
			BamThreadPoolDecodeBamWritePackage<order_type> *, 
			std::vector<BamThreadPoolDecodeBamWritePackage<order_type> *>, 
			BamThreadPoolDecodeBamWritePackageHeapComparator<order_type>
		> 
	> writesPending;

	libmaus::parallel::PosixSpinLock buffersWrittenLock;
	uint64_t buffersWritten;
	libmaus::parallel::LockedBool bamWriteComplete;

	libmaus::parallel::SimpleThreadPool & TP;
	
	bool const verbose;

	typedef libmaus::lz::SimpleCompressedOutputStream<std::ostream> comp_stream_type;
	typedef comp_stream_type::unique_ptr_type comp_stream_ptr_type;
	
	BamThreadPoolDecodeContextBase(
		std::istream & rin,
		uint64_t const numInflateBases,
		uint64_t const numProcessBuffers,
		uint64_t const processBufferSize,
		std::string const & tmpfilenamebase,
		uint64_t const numthreads,
		libmaus::parallel::SimpleThreadPool & rTP,
		bool const rverbose,
		libmaus::lz::CompressorObjectFactory & rcompressorFactory
	)
	:
	  inputLock(),
	  in(rin),
	  readCnt(0),
	  readCompCnt(0),
	  readComplete(false), 
	  inflateBases(numInflateBases),
	  inflateDecompressSpace(numInflateBases * libmaus::lz::BgzfConstants::getBgzfMaxBlockSize()),
	  inflateBasesFreeList(),
	  nextInputBlockId(0),
	  decompressLock(),
	  decompressCnt(0),
	  decompressComplete(false),
	  bamParseQueueLock(),
	  bamparseQueue(),
	  bamparseStall(),
	  bamParseNextBlock(0),
	  bamParseLock(),
	  bamParseCnt(0),
	  bamParseComplete(false),
	  haveheader(false),
	  bamheaderparsestate(),
	  header(),
	  bamParserState(bam_parser_state_read_blocklength),
	  bamBlockSizeRead(0),
	  bamBlockSizeValue(0),
	  bamBlockDataRead(0),
	  processBuffers(numProcessBuffers),
	  processBuffersFreeList(),
	  nextProcessBufferIdIn(0),
	  nextProcessBufferIdOut(0),
	  bamProcessComplete(false),
	  buffersSorted(0),
	  bamSortComplete(false),
	  tmpfilenames(numthreads),
	  tmpfilesCOS(numthreads),
	  #if defined(__linux__)
	  tmpfilesLSPFDOSB(numthreads),
	  #endif
	  compressorFactory(rcompressorFactory),
	  compressedTmpFiles(numthreads),
	  tmpfileblockcntsumsacc(0),
	  writesNext(numthreads,0),
	  writesPending(numthreads),
	  buffersWritten(0),
	  bamWriteComplete(false),
	  TP(rTP),
	  verbose(rverbose)
	{
		for ( uint64_t i = 0; i < inflateBases.size(); ++i )
		{
			libmaus::lz::BgzfInflateBase::unique_ptr_type tptr(new libmaus::lz::BgzfInflateBase);
			inflateBases[i] = UNIQUE_PTR_MOVE(tptr);
			inflateBasesFreeList.enque(i);
		}
		
		for ( uint64_t i = 0; i < processBuffers.size(); ++i )
		{
			BamProcessBuffer::unique_ptr_type tptr(new BamProcessBuffer(processBufferSize));
			processBuffers[i] = UNIQUE_PTR_MOVE(tptr);
			processBuffersFreeList.push_back(i);
		}
		
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			std::ostringstream ostr;
			ostr << tmpfilenamebase << "_" << std::setw(6) << std::setfill('0') << i << std::setw(0) << ".sorttmp";
			std::string const fn = ostr.str();
			tmpfilenames[i] = fn;
			libmaus::util::TempFileRemovalContainer::addTempFile(fn);
			
			bool const local = libmaus::aio::IsKnownLocalFileSystem::isKnownLocalFileSystemCreate(fn);
			std::ostream * out = 0;

			#if defined(__linux__)
			if ( local )
			{
				libmaus::aio::LinuxStreamingPosixFdOutputStream::unique_ptr_type tptr
				(
					new libmaus::aio::LinuxStreamingPosixFdOutputStream(fn,8*1024*1024)
				);
				tmpfilesLSPFDOSB[i] = UNIQUE_PTR_MOVE(tptr);
				out = tmpfilesLSPFDOSB[i].get();
			}
			else
			#endif
			{
				libmaus::aio::CheckedOutputStream::unique_ptr_type tptr(new libmaus::aio::CheckedOutputStream(fn));
				tmpfilesCOS[i] = UNIQUE_PTR_MOVE(tptr);
				out = tmpfilesCOS[i].get();
			}
			
			comp_stream_ptr_type tcptr(new comp_stream_type(*out,compressorFactory));
			compressedTmpFiles[i] = UNIQUE_PTR_MOVE(tcptr);
		}
	}
	
	char * getDecompressSpace(uint64_t const rbaseid)
	{
		return inflateDecompressSpace.begin() + rbaseid * libmaus::lz::BgzfConstants::getBgzfMaxBlockSize();
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeReadPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolDecodeReadPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolDecodeReadPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolDecodeReadPackage<order_type> & RP = *dynamic_cast<BamThreadPoolDecodeReadPackage<order_type> *>(P);
		BamThreadPoolDecodeContextBase<order_type> & contextbase = *(RP.contextbase);
		
		// handle empty input file
		if ( contextbase.readCnt.get() == 0 )
		{
			libmaus::parallel::ScopePosixSpinLock slock(contextbase.inputLock);
			
			if ( contextbase.in.peek() == std::istream::traits_type::eof() )
			{
				contextbase.readComplete.set(true);
				contextbase.decompressComplete.set(true);
				contextbase.bamParseComplete.set(true);
				contextbase.bamProcessComplete.set(true);
				contextbase.bamSortComplete.set(true);
				tpi.terminate();
			}
		}

		if ( ! contextbase.readComplete.get() )
		{
			assert ( contextbase.in.peek() != std::istream::traits_type::eof() );
			
			uint64_t baseid;

			while (  
				(!contextbase.readComplete.get())
				&&
				contextbase.inflateBasesFreeList.trydeque(baseid) 
			)
			{
				uint64_t const nextInputBlockId = (contextbase.nextInputBlockId++);

				#if 0
				uint64_t readCnt;
				uint64_t readCompCnt;
				#endif
				std::pair<uint64_t,uint64_t> blockmeta;
				
				{			
					libmaus::parallel::ScopePosixSpinLock slock(contextbase.inputLock);
					if ( contextbase.readComplete.get() )
						break;
					blockmeta = contextbase.inflateBases[baseid]->readBlock(contextbase.in);

					#if 0
					readCnt = 
					#endif
						contextbase.readCnt++;
					
					contextbase.readCompCnt += blockmeta.first;
					
					#if 0
					readCompCnt = contextbase.readCompCnt.get();
					#endif

					if ( contextbase.in.peek() == std::istream::traits_type::eof() )
					{
						contextbase.readComplete.set(true);
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "read complete." << std::endl;
						contextbase.cerrlock.unlock();
						#endif
					}
				}
				
				#if 0
				if ( readCnt % 16384 == 0 )
				{
					contextbase.cerrlock.lock();
					std::cerr << "[D] input rate " << (readCompCnt / contextbase.inputRtc.getElapsedSeconds())/(1024.0*1024.0) << std::endl;
					contextbase.TP.printStateHistogram(std::cerr);
					contextbase.cerrlock.unlock();
				}
				#endif
				
				BamThreadPoolDecodeDecompressPackage<order_type> * pBTPDDP = RP.contextbase->decompressFreeList.getPackage();
				*pBTPDDP = BamThreadPoolDecodeDecompressPackage<order_type>(0,RP.contextbase,blockmeta,baseid,nextInputBlockId);

				tpi.enque(pBTPDDP);
			}

			contextbase.readSem.post();
		}
		
		RP.contextbase->readFreeList.returnPackage(dynamic_cast<BamThreadPoolDecodeReadPackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeDecompressPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolDecodeDecompressPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolDecodeDecompressPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolDecodeDecompressPackage<order_type> & RP = *dynamic_cast<BamThreadPoolDecodeDecompressPackage<order_type> *>(P);
		BamThreadPoolDecodeContextBase<order_type> & contextbase = *(RP.contextbase);

		char * const decompressSpace = contextbase.getDecompressSpace(RP.baseid);
		contextbase.inflateBases[RP.baseid]->decompressBlock(decompressSpace,RP.blockmeta);

		libmaus::parallel::ScopePosixSpinLock ldecompressLock(contextbase.decompressLock);
		contextbase.decompressCnt += 1;
		if ( contextbase.readComplete.get() && (contextbase.decompressCnt == contextbase.readCnt) )
		{
			contextbase.decompressComplete.set(true);

			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "decompress complete." << std::endl;
			contextbase.cerrlock.unlock();
			#endif
		}

		{
			// generate parse info object
			BamThreadPoolDecodeBamParseQueueInfo qinfo(0,RP.blockmeta,RP.baseid,RP.blockid);

			// enque the parse info object
			libmaus::parallel::ScopePosixSpinLock lbamParseQueueLock(contextbase.bamParseQueueLock);
			contextbase.bamparseQueue.push(qinfo);
			
			// check whether this is the next block for parsing
			if ( contextbase.bamParseNextBlock == contextbase.bamparseQueue.top().blockid )
			{
				BamThreadPoolDecodeBamParseQueueInfo const bqinfo = contextbase.bamparseQueue.top();
				contextbase.bamparseQueue.pop();
				
				BamThreadPoolDecodeBamParsePackage<order_type> * pBTPDBPP = RP.contextbase->bamParseFreeList.getPackage();
				*pBTPDBPP = BamThreadPoolDecodeBamParsePackage<order_type>(
					bqinfo.packageid,
					RP.contextbase,
					bqinfo.blockmeta,
					bqinfo.baseid,
					bqinfo.blockid
				);


				tpi.enque(pBTPDBPP);
			}
		}		

		RP.contextbase->decompressFreeList.returnPackage(dynamic_cast<BamThreadPoolDecodeDecompressPackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamParsePackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	template<unsigned int len>
	static uint64_t decodeLittleEndianInteger(uint8_t const * pa)
	{
		#if defined(LIBMAUS_HAVE_i386)
		if ( len == 4 )
		{
			return *reinterpret_cast<uint32_t const *>(pa);
		}
		#else
		uint64_t v = 0;
		for ( unsigned int shift = 0; shift < 8*len; shift += 8 )
			v |= static_cast<uint64_t>(*(pa++)) << shift;
		return v;
		#endif
	}

	virtual ~BamThreadPoolDecodeBamParsePackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolDecodeBamParsePackage<order_type> *>(P) != 0 );
		
		BamThreadPoolDecodeBamParsePackage<order_type> & RP = *dynamic_cast<BamThreadPoolDecodeBamParsePackage<order_type> *>(P);
		BamThreadPoolDecodeContextBase<order_type> & contextbase = *(RP.contextbase);
		
		char * const decompressSpace = contextbase.getDecompressSpace(RP.baseid);

		uint8_t const * pa = reinterpret_cast<uint8_t const *>(decompressSpace);
		uint8_t const * pc = pa + RP.blockmeta.second;

		if ( (! contextbase.haveheader.get()) && (pa != pc) )
		{			
			::libmaus::util::GetObject<uint8_t const *> G(pa);
			std::pair<bool,uint64_t> const P = ::libmaus::bambam::BamHeader::parseHeader(G,contextbase.bamheaderparsestate,pc-pa);

			// header complete?
			if ( P.first )
			{
				contextbase.header.init(contextbase.bamheaderparsestate);
				contextbase.haveheader.set(true);
				pa += P.second;
				
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << contextbase.header.text;
				contextbase.cerrlock.unlock();
				#endif
			}
			else
			{
				pa += P.second;
			}
		}

		bool bufferAvailable = !(contextbase.processBuffersFreeList.empty()); 
		int64_t freeBufferId = bufferAvailable ? contextbase.processBuffersFreeList.dequeFront() : -1;
		std::vector<uint64_t> finishedBuffers;
		
		bool stall = freeBufferId < 0;
		bool running = (pa != pc) && (!stall);
		
		while ( running )
		{
			assert ( contextbase.haveheader.get() );
			BamProcessBuffer * bamProcessBuffer = contextbase.processBuffers[freeBufferId].get();
			
			switch ( contextbase.bamParserState )
			{
				case BamThreadPoolDecodeContextBase<order_type>::bam_parser_state_read_blocklength:
				{
					if ( contextbase.bamBlockSizeRead == 0 )
					{
						uint32_t lblocksize;
						while ( (pc-pa >= 4) && ((pc-pa) >= (4+(lblocksize=decodeLittleEndianInteger<4>(pa)) )) )
						{
							if ( bamProcessBuffer->put(pa+4,lblocksize) )
							{
								pa += 4 + lblocksize;
							}
							else
							{
								finishedBuffers.push_back(freeBufferId);
								
								bufferAvailable = !(contextbase.processBuffersFreeList.empty());
								freeBufferId = bufferAvailable ? contextbase.processBuffersFreeList.dequeFront() : -1;
								
								if ( freeBufferId < 0 )
								{
									running = false;
									stall = true;
									break;
								}
								else
								{
									bamProcessBuffer = contextbase.processBuffers[freeBufferId].get();
								}
							}
						}						
					}
					while ( running && (pa != pc) && (contextbase.bamBlockSizeRead < 4) )
					{
						contextbase.bamBlockSizeValue |= static_cast<uint32_t>(*(pa++)) << (8*((contextbase.bamBlockSizeRead++)));
					}
					if ( running && (contextbase.bamBlockSizeRead == 4) )
					{
						contextbase.bamParserState = contextbase.bam_parser_state_read_blockdata;
						contextbase.bamBlockDataRead = 0;

						if ( contextbase.bamBlockSizeValue > contextbase.bamGatherBuffer.D.size() )
							contextbase.bamGatherBuffer.D = libmaus::bambam::BamAlignment::D_array_type(contextbase.bamBlockSizeValue);
					}
					break;
				}
				case BamThreadPoolDecodeContextBase<order_type>::bam_parser_state_read_blockdata:
				{
					assert ( pa != pc );
				
					uint64_t const skip = 
						std::min(
							static_cast<uint64_t>(pc-pa),
							static_cast<uint64_t>(contextbase.bamBlockSizeValue - contextbase.bamBlockDataRead)
						);
							
					std::copy(pa,pa+skip,contextbase.bamGatherBuffer.D.begin()+contextbase.bamBlockDataRead);
						
					contextbase.bamBlockDataRead += skip;
					pa += skip;

					if ( contextbase.bamBlockDataRead == contextbase.bamBlockSizeValue )
					{
						contextbase.bamGatherBuffer.blocksize = contextbase.bamBlockSizeValue;

						if ( bamProcessBuffer->put(
							contextbase.bamGatherBuffer.D.begin(),
							contextbase.bamGatherBuffer.blocksize)
						)
						{
							contextbase.bamBlockSizeRead = 0;
							contextbase.bamBlockSizeValue = 0;
							contextbase.bamParserState = contextbase.bam_parser_state_read_blocklength;
						}
						else
						{
							// queue finished buffer
							finishedBuffers.push_back(freeBufferId);

							// rewind to process this data again
							pa -= skip;
							contextbase.bamBlockDataRead -= skip;

							// check for next free buffer								
							bufferAvailable = !(contextbase.processBuffersFreeList.empty());
							freeBufferId = bufferAvailable ? contextbase.processBuffersFreeList.dequeFront() : -1;

							// no more free buffers, stall
							if ( freeBufferId < 0 )
							{
								stall = true;
								running = false;
							}
							else
								bamProcessBuffer = contextbase.processBuffers[freeBufferId].get();							
						}
					}
					break;
				}
			}
			
			running = running && (pa != pc) && (!stall);
		}

		if ( ! stall )
		{
			// add inflate object to free list if it is no longer required
			contextbase.inflateBasesFreeList.enque(RP.baseid);

			// enque next read package when stream is ready for reading
			if ( contextbase.readSem.trywait() )
			{
				BamThreadPoolDecodeReadPackage<order_type> * pBTPDRP = RP.contextbase->readFreeList.getPackage();
				*pBTPDRP = BamThreadPoolDecodeReadPackage<order_type>(0,RP.contextbase);
				tpi.enque(pBTPDRP);
			}
		}

		if ( ! stall )
		{
			
			libmaus::parallel::ScopePosixSpinLock lbamParseLock(contextbase.bamParseLock);
			contextbase.bamParseCnt += 1;
			if ( contextbase.decompressComplete.get() && (contextbase.bamParseCnt == contextbase.decompressCnt) )
			{
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "bamParse complete." << std::endl;
				contextbase.cerrlock.unlock();
				#endif
				contextbase.bamParseComplete.set(true);

				// queue final finished buffer
				finishedBuffers.push_back(freeBufferId);
			}
			else
			{
				// reinsert unfinished buffer
				contextbase.processBuffersFreeList.push_front(freeBufferId);	
			}
		}

		if ( stall )
		{
			// move remaining data to start of buffer
			memmove(reinterpret_cast<uint8_t *>(decompressSpace),pa,pc-pa);

			RP.blockmeta.second = pc-pa;

			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "stalling, rest " << pc-pa << std::endl;
			contextbase.cerrlock.unlock();
			#endif
		
			// stall package
			BamThreadPoolDecodeBamParseQueueInfo qinfo(
				RP.packageid,
				RP.blockmeta,
				RP.baseid,
				RP.blockid
			);

			// enque the parse info object in the stall queue
			libmaus::parallel::ScopePosixSpinLock lbamParseQueueLock(contextbase.bamParseQueueLock);
			contextbase.bamparseStall.push(qinfo);
		}

		// queue finished buffers
		if ( finishedBuffers.size() )
		{
			for ( uint64_t i = 0; i < finishedBuffers.size(); ++i )
			{
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "queueing process block." << std::endl;
				contextbase.cerrlock.unlock();
				#endif
			
				BamThreadPoolDecodeBamProcessQueueInfo qinfo(
					0,
					finishedBuffers[i],
					contextbase.nextProcessBufferIdIn++	
				);

				BamThreadPoolDecodeBamProcessPackage<order_type> * pprocpack = RP.contextbase->bamProcessFreeList.getPackage();
				*pprocpack = BamThreadPoolDecodeBamProcessPackage<order_type>(
					qinfo.packageid,
					RP.contextbase,
					qinfo.baseid,
					qinfo.blockid
				);
				
				tpi.enque(pprocpack);
			}			
		}
		
		if ( ! stall )
		{
			// process next bam parse object if any is available
			{
				libmaus::parallel::ScopePosixSpinLock lbamParseQueueLock(contextbase.bamParseQueueLock);
				contextbase.bamParseNextBlock += 1;

				if ( contextbase.bamParseNextBlock == contextbase.bamparseQueue.top().blockid )
				{
					BamThreadPoolDecodeBamParseQueueInfo const bqinfo = contextbase.bamparseQueue.top();
					contextbase.bamparseQueue.pop();

					BamThreadPoolDecodeBamParsePackage<order_type> * pBTPDBPP = RP.contextbase->bamParseFreeList.getPackage();
					*pBTPDBPP = BamThreadPoolDecodeBamParsePackage<order_type>(
						bqinfo.packageid,
						RP.contextbase,
						bqinfo.blockmeta,
						bqinfo.baseid,
						bqinfo.blockid
					);
					
					tpi.enque(pBTPDBPP);
				}
			}
		}

		RP.contextbase->bamParseFreeList.returnPackage(dynamic_cast<BamThreadPoolDecodeBamParsePackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamProcessPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolDecodeBamProcessPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolDecodeBamProcessPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolDecodeBamProcessPackage<order_type> & RP = *dynamic_cast<BamThreadPoolDecodeBamProcessPackage<order_type> *>(P);
		BamThreadPoolDecodeContextBase<order_type> & contextbase = *(RP.contextbase);

		BamProcessBuffer * processBuffer = contextbase.processBuffers[RP.baseid].get();
		std::reverse(processBuffer->pc,processBuffer->pa);
		
		typename BamSortInfo<order_type>::shared_ptr_type sortinfo(new BamSortInfo<order_type>(processBuffer,contextbase.TP.threads.size()));

		for ( uint64_t i = 0; i < contextbase.TP.threads.size(); ++i )
		{
			BamThreadPoolDecodeBamSortPackage<order_type> * spack = RP.contextbase->bamSortFreeList.getPackage();
			*spack = 
				BamThreadPoolDecodeBamSortPackage<order_type>(
					0, RP.contextbase, RP.baseid, RP.blockid, sortinfo, 
					BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_base,
					i,0,0
				);
			tpi.enque(spack);
		}

		// increment number of processed buffers
		{
		libmaus::parallel::ScopePosixSpinLock lnextProcessBufferIdOutLock(contextbase.nextProcessBufferIdOutLock);
		contextbase.nextProcessBufferIdOut += 1;
		}

		// check whether processing of package type is complete
		{
			libmaus::parallel::ScopePosixSpinLock lnextProcessBufferIdOutLock(contextbase.nextProcessBufferIdOutLock);
			if ( 
				contextbase.bamParseComplete.get() && (contextbase.nextProcessBufferIdOut == contextbase.nextProcessBufferIdIn) 
			)
			{
				contextbase.bamProcessComplete.set(true);
				#if 0
				std::cerr << "bamProcess complete" << std::endl;
				#endif
			}
		}

		// return package
		RP.contextbase->bamProcessFreeList.returnPackage(dynamic_cast<BamThreadPoolDecodeBamProcessPackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamSortPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolDecodeBamSortPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolDecodeBamSortPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolDecodeBamSortPackage<order_type> & RP = *dynamic_cast<BamThreadPoolDecodeBamSortPackage<order_type> *>(P);
		BamThreadPoolDecodeContextBase<order_type> & contextbase = *(RP.contextbase);
		typename BamSortInfo<order_type>::sort_type & sortControl = *(RP.sortinfo->sortControl);
		
		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "BamThreadPoolDecodeBamSortPackageDispatcher " 
			<< "packagetype=" << static_cast<int>(RP.packagetype) 
			<< ",sort_base_id=" << RP.sort_base_id
			<< ",sort_merge_id=" << RP.sort_merge_id
			<< ",sort_submerge_id=" << RP.sort_submerge_id
			<< std::endl;
		contextbase.cerrlock.unlock();
		#endif
		
		switch ( RP.packagetype )
		{
			case BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_base:
			{
				libmaus::sorting::ParallelStableSort::BaseSortRequestSet<uint64_t *,typename BamSortInfo<order_type>::order_type> & baseSortRequests =
					sortControl.baseSortRequests;

				baseSortRequests.baseSortRequests[RP.sort_base_id].dispatch();

				uint64_t const finished = ++(baseSortRequests.requestsFinished);

				if ( finished == baseSortRequests.baseSortRequests.size() )
				{
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "BamThreadPoolDecodeBamSortPackageDispatcher done" << std::endl;
					contextbase.cerrlock.unlock();	
					#endif

					BamThreadPoolDecodeBamSortPackage<order_type> * spack = RP.contextbase->bamSortFreeList.getPackage();
					*spack = BamThreadPoolDecodeBamSortPackage<order_type>(0, RP.contextbase, RP.baseid, RP.blockid, 
						RP.sortinfo, BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_merge_level,0,0,0);
					tpi.enque(spack);
				}
				break;
			}
			case BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_merge_level:
			{
				if ( sortControl.mergeLevels.levelsFinished.get() == sortControl.mergeLevels.levels.size() )
				{
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "sorted buffer " << contextbase.processBuffers[RP.baseid]->pa-contextbase.processBuffers[RP.baseid]->pc << std::endl;
					contextbase.cerrlock.unlock();
					#endif

					BamProcessBuffer * processBuffer = contextbase.processBuffers[RP.baseid].get();
					order_type BAPC(processBuffer->ca);
					
					BamProcessBuffer::pointer_type * pa =
						sortControl.needCopyBack ?
						( processBuffer->pc - (processBuffer->pa-processBuffer->pc) )
						:
						processBuffer->pc;
					BamProcessBuffer::pointer_type * pe =
						sortControl.needCopyBack ? processBuffer->pc : processBuffer->pa;

					uint64_t const numthreads = contextbase.TP.threads.size();
					std::vector < std::pair<BamProcessBuffer::pointer_type const *, BamProcessBuffer::pointer_type const *> > writepackets(numthreads);
					uint64_t const numalgn = pe-pa;
					uint64_t const prealperpack = (numalgn + numthreads - 1)/numthreads;
					uint64_t const alperpack =
						((prealperpack + contextbase.alperpackalign - 1)/contextbase.alperpackalign) * contextbase.alperpackalign;
					uint64_t allow = 0;

					for ( uint64_t i = 0; i < numthreads; ++i )
					{
						uint64_t const alhigh = std::min(allow+alperpack,numalgn);
						
						writepackets[i] = std::pair<BamProcessBuffer::pointer_type const *, BamProcessBuffer::pointer_type const *>(
							pa + allow, pa + alhigh
						);

						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "writepacket " << allow << " " << alhigh << std::endl;
						contextbase.cerrlock.unlock();
						#endif
						
						allow = alhigh;
					}

					BamBlockWriteInfo::shared_ptr_type blockWriteInfo(new BamBlockWriteInfo(processBuffer,writepackets));

					{
						libmaus::parallel::ScopePosixSpinLock ltmpfileblockslock(contextbase.tmpfileblockslock);

						while ( ! (RP.blockid < contextbase.tmpfileblocks.size() ) )
							contextbase.tmpfileblocks.push_back(
								std::vector< 
									libmaus::lz::SimpleCompressedStreamInterval
								>(
									numthreads
								)
							);
						
						while ( ! (RP.blockid < contextbase.tmpfileblockcnts.size() ) )
						{
							contextbase.tmpfileblockcnts.push_back(
								std::vector<uint64_t>(numthreads)
							);						
						}
						
						while ( ! (RP.blockid < contextbase.tmpfileblockcntsums.size() ) )
						{
							contextbase.tmpfileblockcntsums.push_back(0);
						}
					}

					libmaus::parallel::ScopePosixSpinLock lwritesPendingLock(contextbase.writesPendingLock);
					for ( uint64_t i = 0; i < numthreads; ++i )
					{
						BamThreadPoolDecodeBamWritePackage<order_type> * pBTPDBPP = RP.contextbase->bamWriteFreeList.getPackage();
						*pBTPDBPP = BamThreadPoolDecodeBamWritePackage<order_type>(
							0,RP.contextbase,RP.baseid,RP.blockid,
							blockWriteInfo,
							i
						);
						contextbase.writesPending[i].push(pBTPDBPP);
						
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "queuing " << RP.blockid << "," << i << " as pending." << std::endl;
						contextbase.cerrlock.unlock();
						#endif
					}
					for ( uint64_t i = 0; i < numthreads; ++i )
					{
						assert ( ! contextbase.writesPending[i].empty() );

						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "checking queue for thread " << i << " top is " << contextbase.writesPending[i].top()->blockid << std::endl;
						contextbase.cerrlock.unlock();
						#endif
						
						if ( 
							contextbase.writesPending[i].top()->blockid ==
							contextbase.writesNext[i]
						)
						{
							BamThreadPoolDecodeBamWritePackage<order_type> * pBTPDBPP = 
								contextbase.writesPending[i].top();
							contextbase.writesPending[i].pop();
							
							tpi.enque(pBTPDBPP);
						}
					}

					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "needCopyBack " << sortControl.needCopyBack << std::endl;
					contextbase.cerrlock.unlock();
					
					contextbase.cerrlock.lock();	
					uint8_t const * prev = 0;
					bool ok = true;
					for ( BamProcessBuffer::pointer_type * pc = pa ; pc != pe; ++pc )
					{
						uint8_t const * c = processBuffer->ca + (*pc)+4;
						#if 0
						std::cerr 
							<< libmaus::bambam::BamAlignmentDecoderBase::getReadName(c) 
							<< "\t"
							<< libmaus::bambam::BamAlignmentDecoderBase::getRefID(c)
							<< "\t"
							<< libmaus::bambam::BamAlignmentDecoderBase::getPos(c)
							<< std::endl;
						#endif
							
						if ( prev )
						{
							ok = ok && ( ! (BAPC( *pc , *(pc-1) )) );
							// assert ( ! (BAPC( *pc , *(pc-1) )) );
							// std::cerr << BAPC( *pc , *(pc-1) ) << std::endl;
						}
							
						prev = c;
					}
					std::cerr << "order: " << (ok?"ok":"failed") << std::endl;
					contextbase.cerrlock.unlock();
					#endif
					
					#if 0
					// return buffer to free list
					contextbase.processBuffers[RP.baseid]->reset();
					contextbase.processBuffersFreeList.push_back(RP.baseid);
					
					// reinsert stalled parse package
					{
						libmaus::parallel::ScopePosixSpinLock lbamParseQueueLock(contextbase.bamParseQueueLock);
						
						if ( contextbase.bamparseStall.size() )
						{
							BamThreadPoolDecodeBamParseQueueInfo const bqinfo = contextbase.bamparseStall.top();
							contextbase.bamparseStall.pop();

							BamThreadPoolDecodeBamParsePackage * pBTPDBPP = RP.contextbase->bamParseFreeList.getPackage();
							*pBTPDBPP = BamThreadPoolDecodeBamParsePackage(
								bqinfo.packageid,
								RP.contextbase,
								bqinfo.blockmeta,
								bqinfo.baseid,
								bqinfo.blockid
							);

							#if 0
							contextbase.cerrlock.lock();
							std::cerr << "reinserting stalled block " << bqinfo.blockid << std::endl;
							contextbase.cerrlock.unlock();
							#endif
				
							tpi.enque(pBTPDBPP);			
						}
					}
					#endif

					libmaus::parallel::ScopePosixSpinLock lnextProcessBufferIdOutLock(contextbase.nextProcessBufferIdOutLock);
					libmaus::parallel::ScopePosixSpinLock lbuffersSortedLock(contextbase.buffersSortedLock);					
					contextbase.buffersSorted += 1;
					if ( 
						contextbase.bamProcessComplete.get()
						&&
						contextbase.nextProcessBufferIdOut == contextbase.buffersSorted
					)
					{
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "bamSort complete" << std::endl;
						contextbase.cerrlock.unlock();
						#endif

						contextbase.bamSortComplete.set(true);
						// tpi.terminate();
					}
				}
				else
				{
					// search split points
					sortControl.mergeLevels.levels[RP.sort_merge_id].dispatch();
					
					// call merging
					for ( uint64_t i = 0; i < sortControl.mergeLevels.levels[RP.sort_merge_id].mergeRequests.size(); ++i )
					{
						BamThreadPoolDecodeBamSortPackage<order_type> * spack = RP.contextbase->bamSortFreeList.getPackage();
						*spack = BamThreadPoolDecodeBamSortPackage<order_type>(0, RP.contextbase, RP.baseid, RP.blockid, 
							RP.sortinfo, 
							BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_merge_sublevel,
							0,RP.sort_merge_id,i);
						tpi.enque(spack);	
					}
				}
				break;
			}
			case BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_merge_sublevel:
			{
				sortControl.mergeLevels.levels[RP.sort_merge_id].mergeRequests[RP.sort_submerge_id].dispatch();

				uint64_t const finished = ++(sortControl.mergeLevels.levels[RP.sort_merge_id].requestsFinished);
				
				if ( finished == sortControl.mergeLevels.levels[RP.sort_merge_id].mergeRequests.size() )
				{
					sortControl.mergeLevels.levelsFinished++;

					BamThreadPoolDecodeBamSortPackage<order_type> * spack = RP.contextbase->bamSortFreeList.getPackage();
					*spack = BamThreadPoolDecodeBamSortPackage<order_type>(0, RP.contextbase, RP.baseid, RP.blockid, 
						RP.sortinfo, 
						BamThreadPoolDecodeBamSortPackage<order_type>::sort_package_type_merge_level,
						0,sortControl.mergeLevels.levelsFinished.get(),0
					);
					tpi.enque(spack);	
				}
				
				break;
			}
		}

		RP.sortinfo.reset();
		RP.contextbase->bamSortFreeList.returnPackage(dynamic_cast<BamThreadPoolDecodeBamSortPackage<order_type> *>(P));		
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeBamWritePackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolDecodeBamWritePackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolDecodeBamWritePackage<order_type> *>(P) != 0 );
		
		BamThreadPoolDecodeBamWritePackage<order_type> & RP = *dynamic_cast<BamThreadPoolDecodeBamWritePackage<order_type> *>(P);
		BamThreadPoolDecodeContextBase<order_type> & contextbase = *(RP.contextbase);
		BamBlockWriteInfo * writeinfo = RP.writeinfo.get();		
		BamProcessBuffer * buffer = writeinfo->processBuffer;
		
		uint8_t const * ca = buffer->ca;
		BamProcessBuffer::pointer_type const * pa = writeinfo->packets[RP.write_block_id].first;
		BamProcessBuffer::pointer_type const * pe = writeinfo->packets[RP.write_block_id].second;
		libmaus::lz::SimpleCompressedOutputStream<std::ostream> * compout =
			contextbase.compressedTmpFiles[RP.write_block_id].get();
	
		#if 0
		order_type BAPC(buffer->ca);
		BamProcessBuffer::pointer_type const * prev = 0;
		bool gok = true;
		#endif
		
		#if 0
		if ( pa != pe )
		{
			::libmaus::bambam::BamFormatAuxiliary aux;
			
			contextbase.cerrlock.lock();
			std::cerr << "first " 
				<< 
					libmaus::bambam::BamAlignmentDecoderBase::formatAlignment(
						ca + *pa + 4,
						(static_cast<uint32_t>((ca + *pa)[0]) << 0) |
						(static_cast<uint32_t>((ca + *pa)[1]) << 8) |
						(static_cast<uint32_t>((ca + *pa)[2]) << 16) |
						(static_cast<uint32_t>((ca + *pa)[3]) << 24),
						contextbase.header,
						aux
					)
				<< std::endl;
			contextbase.cerrlock.unlock();
		}
		#endif
		
		std::pair<uint64_t,uint64_t> const preoff = compout->getOffset();
		#if 0
		assert ( preoff.second == 0 );
		#endif
		for ( BamProcessBuffer::pointer_type const * pc = pa; pc != pe; ++pc )
		{
			uint64_t const off = *pc;
			uint8_t const * c = ca + off;
			
			uint32_t len = 0;
			for ( unsigned int i = 0; i < 4; ++i )
				len |= (static_cast<uint32_t>(c[i]) << (i*8));
		
			#if 0	
			if ( (pc-pa) % contextbase.alperpackalign == 0 )
			{
				compout->flush();
				assert ( compout->getOffset().second == 0 );
			}
			#endif
			compout->write(reinterpret_cast<char const *>(c),len+4);
			
			#if 0
			if ( prev )
			{
				bool const ok = !BAPC(*pc,*prev);
				if ( ! ok )
				{
					contextbase.cerrlock.lock();
					std::cerr << "comp failed" << std::endl;
					contextbase.cerrlock.unlock();
				}
				
				gok = gok && ok;
			}
			
			prev = pc;
			#endif
		}
		std::pair<uint64_t,uint64_t> const postoff = compout->getOffset();
		compout->flush();

		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "interval comp " << (gok?"ok":"failed") << std::endl;
		contextbase.cerrlock.unlock();
		
		contextbase.cerrlock.lock();
		std::cerr << "BamWrite " 
			<< RP.blockid << "," << RP.write_block_id 
			<< "," << (writeinfo->packets[RP.write_block_id].second-writeinfo->packets[RP.write_block_id].first)
			<< std::endl;
		contextbase.cerrlock.unlock();
		#endif

		{
			libmaus::parallel::ScopePosixSpinLock ltmpfileblockslock(contextbase.tmpfileblockslock);
			contextbase.tmpfileblocks[RP.blockid][RP.write_block_id] =
				libmaus::lz::SimpleCompressedStreamInterval(preoff,postoff);
			contextbase.tmpfileblockcnts[RP.blockid][RP.write_block_id] = (pe-pa);
			contextbase.tmpfileblockcntsums[RP.blockid] += (pe-pa);
		}

		// enque next block for stream writing (if any is present)
		{
			libmaus::parallel::ScopePosixSpinLock lwritesPendingLock(contextbase.writesPendingLock);
			contextbase.writesNext[RP.write_block_id]++;

			if ( 
				contextbase.writesPending[RP.write_block_id].size() 
				&&
				contextbase.writesPending[RP.write_block_id].top()->blockid ==
				contextbase.writesNext[RP.write_block_id]
			)
			{
				BamThreadPoolDecodeBamWritePackage<order_type> * pBTPDBPP = 
					contextbase.writesPending[RP.write_block_id].top();
				contextbase.writesPending[RP.write_block_id].pop();
				
				tpi.enque(pBTPDBPP);
			}
		}

		uint64_t const packetsWritten = ++(writeinfo->packetsWritten);
		
		if ( packetsWritten == writeinfo->packets.size() )
		{
			{
				libmaus::parallel::ScopePosixSpinLock ltmpfileblockcntsumsacclock(contextbase.tmpfileblockcntsumsacclock);
				contextbase.tmpfileblockcntsumsacc += contextbase.tmpfileblockcntsums[RP.blockid];
				
				if ( contextbase.verbose )
				{
					contextbase.cerrlock.lock();
					std::cerr << "[V] " << contextbase.tmpfileblockcntsumsacc << std::endl;
					contextbase.cerrlock.unlock();		
				}
			}
		
			// return buffer to free list
			contextbase.processBuffers[RP.baseid]->reset();
			contextbase.processBuffersFreeList.push_back(RP.baseid);
					
			// reinsert stalled parse package
			{
				libmaus::parallel::ScopePosixSpinLock lbamParseQueueLock(contextbase.bamParseQueueLock);
				
				if ( contextbase.bamparseStall.size() )
				{
					BamThreadPoolDecodeBamParseQueueInfo const bqinfo = contextbase.bamparseStall.top();
					contextbase.bamparseStall.pop();

					BamThreadPoolDecodeBamParsePackage<order_type> * pBTPDBPP = RP.contextbase->bamParseFreeList.getPackage();
					*pBTPDBPP = BamThreadPoolDecodeBamParsePackage<order_type>(
						bqinfo.packageid,
						RP.contextbase,
						bqinfo.blockmeta,
						bqinfo.baseid,
						bqinfo.blockid
					);

					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "reinserting stalled block " << bqinfo.blockid << std::endl;
					contextbase.cerrlock.unlock();
					#endif
		
					tpi.enque(pBTPDBPP);			
				}
			}			

			libmaus::parallel::ScopePosixSpinLock lbuffersWrittenLock(contextbase.buffersWrittenLock);
			libmaus::parallel::ScopePosixSpinLock lbuffersSortedLock(contextbase.buffersSortedLock);
			contextbase.buffersWritten += 1;
			
			if ( 
				contextbase.bamSortComplete.get()
				&&
				contextbase.buffersWritten == contextbase.buffersSorted 
			)
			{
				for ( uint64_t i = 0; i < contextbase.compressedTmpFiles.size(); ++i )
				{
					contextbase.compressedTmpFiles[i]->flush();
					contextbase.compressedTmpFiles[i].reset();
					if ( contextbase.tmpfilesCOS[i] )
					{
						contextbase.tmpfilesCOS[i]->flush();
						contextbase.tmpfilesCOS[i].reset();
					}
					#if defined(__linux__)
					if ( contextbase.tmpfilesLSPFDOSB[i] )
					{
						contextbase.tmpfilesLSPFDOSB[i]->flush();
						contextbase.tmpfilesLSPFDOSB[i].reset();
					}
					#endif
				}

				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "bamWrite done." << std::endl;
				contextbase.cerrlock.unlock();
				#endif
				
				contextbase.bamWriteComplete.set(true);
				tpi.terminate();
			}
		}

		// return package to free list
		RP.contextbase->bamWriteFreeList.returnPackage(dynamic_cast<BamThreadPoolDecodeBamWritePackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolDecodeContext : public BamThreadPoolDecodeContextBase<_order_type>
{
	typedef _order_type order_type;

	libmaus::parallel::SimpleThreadPool & TP;

	BamThreadPoolDecodeContext(
		std::istream & rin, 
		uint64_t const numInflateBases,
		uint64_t const numProcessBuffers,
		uint64_t const processBufferSize,
		std::string const & tmpfilenamebase,
		uint64_t const numthreads,
		libmaus::parallel::SimpleThreadPool & rTP,
		bool const verbose,
		libmaus::lz::CompressorObjectFactory & rcompressorFactory
	)
	: BamThreadPoolDecodeContextBase<order_type>(
		rin,numInflateBases,numProcessBuffers,
		processBufferSize,tmpfilenamebase,numthreads,rTP,verbose,
		rcompressorFactory
	), TP(rTP) 
	{
	
	}
	
	void startup()
	{
		BamThreadPoolDecodeContextBase<order_type>::inputRtc.start();
	
		assert ( BamThreadPoolDecodeContextBase<order_type>::inflateBases.size() );

		BamThreadPoolDecodeReadPackage<order_type> * pBTPDRP = BamThreadPoolDecodeContextBase<order_type>::readFreeList.getPackage();
		*pBTPDRP = BamThreadPoolDecodeReadPackage<order_type>(0,this);
		TP.enque(pBTPDRP);
	}
};


template<typename _order_type>
MergeInfo produceSortedBlocks(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::lz::CompressorObjectFactory & rcompressorFactory
)
{
	typedef _order_type order_type;

	std::string const tmpfilenamebase = 
		arginfo.getUnparsedValue("tmpfileprefix",arginfo.getDefaultTmpFileName());
		
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());

	uint64_t const numthreads = arginfo.getValue<unsigned int>("numthreads", libmaus::parallel::NumCpus::getNumLogicalProcessors());

	// thread pool	
	libmaus::parallel::SimpleThreadPool TP(numthreads);

	// package dispatchers
	BamThreadPoolDecodeReadPackageDispatcher<order_type> readdispatcher;
	TP.registerDispatcher(BamThreadPoolDecodeContext<order_type>::bamthreadpooldecodecontextbase_dispatcher_id_read,&readdispatcher);
	BamThreadPoolDecodeDecompressPackageDispatcher<order_type> decompressdispatcher;
	TP.registerDispatcher(BamThreadPoolDecodeContext<order_type>::bamthreadpooldecodecontextbase_dispatcher_id_decompress,&decompressdispatcher);
	BamThreadPoolDecodeBamParsePackageDispatcher<order_type> bamparsedispatcher;
	TP.registerDispatcher(BamThreadPoolDecodeContext<order_type>::bamthreadpooldecodecontextbase_dispatcher_id_bamparse,&bamparsedispatcher);
	BamThreadPoolDecodeBamProcessPackageDispatcher<order_type> bamprocessdispatcher;
	TP.registerDispatcher(BamThreadPoolDecodeContext<order_type>::bamthreadpooldecodecontextbase_dispatcher_id_bamprocess,&bamprocessdispatcher);
	BamThreadPoolDecodeBamSortPackageDispatcher<order_type> bamsortdispatcher;
	TP.registerDispatcher(BamThreadPoolDecodeContext<order_type>::bamthreadpooldecodecontextbase_dispatcher_id_bamsort,&bamsortdispatcher);
	BamThreadPoolDecodeBamWritePackageDispatcher<order_type> bamwritedispatcher;
	TP.registerDispatcher(BamThreadPoolDecodeContext<order_type>::bamthreadpooldecodecontextbase_dispatcher_id_bamwrite,&bamwritedispatcher);

	uint64_t numProcessBuffers = arginfo.getValueUnsignedNumeric<uint64_t>("processbuffers",getDefaultProcessBuffers());; // 2*numthreads;
	uint64_t processBufferMemory = arginfo.getValueUnsignedNumeric<uint64_t>("mem",getDefaultMem());
	uint64_t processBufferSize = (processBufferMemory + numProcessBuffers-1)/numProcessBuffers;
	assert ( processBufferSize <= std::numeric_limits<BamProcessBuffer::pointer_type>::max() );

	libmaus::aio::PosixFdInputStream PFIS(STDIN_FILENO,64*1024);
	BamThreadPoolDecodeContext<order_type> context(PFIS,16*numthreads /* inflate bases */,numProcessBuffers,processBufferSize,tmpfilenamebase,numthreads,TP,verbose,rcompressorFactory);
	context.startup();
	
	TP.join();
	
	if ( verbose )
	{
		std::cerr << "[V] produced " << context.tmpfileblockcntsums.size() << " blocks: " << std::endl;
		for ( uint64_t i = 0; i < context.tmpfileblockcntsums.size(); ++i )
			std::cerr << "[V] block[" << i << "]=" << context.tmpfileblockcntsums[i] << std::endl;
	}

	MergeInfo const mergeinfo = context.getMergeInfo();

	return mergeinfo;
}

template<typename _order_type>
void mergeSortedBlocksNoThreadPool(
	libmaus::util::ArgInfo const & arginfo, MergeInfo const & mergeinfo, std::string const & neworder,
	libmaus::lz::DecompressorObjectFactory & decfact
)
{
	typedef _order_type order_type;

	uint64_t const numthreads = arginfo.getValue<unsigned int>("numthreads", libmaus::parallel::NumCpus::getNumLogicalProcessors());

	libmaus::autoarray::AutoArray<libmaus::aio::CheckedInputStream::unique_ptr_type> inputfiles(mergeinfo.tmpfilenames.size());
	for ( uint64_t i = 0; i < inputfiles.size(); ++i )
	{
		libmaus::aio::CheckedInputStream::unique_ptr_type tptr(new libmaus::aio::CheckedInputStream(mergeinfo.tmpfilenames[i]));
		inputfiles[i] = UNIQUE_PTR_MOVE(tptr);
	}
	
	libmaus::autoarray::AutoArray<libmaus::lz::SimpleCompressedConcatInputStream<std::istream>::unique_ptr_type> concfiles(mergeinfo.tmpfileblocks.size());
	for ( uint64_t i = 0; i < mergeinfo.tmpfileblocks.size(); ++i )
	{
		std::vector<libmaus::lz::SimpleCompressedConcatInputStreamFragment<std::istream> > fragments;
		
		for ( uint64_t j = 0; j < mergeinfo.tmpfileblocks[i].size(); ++j )
			fragments.push_back(
				libmaus::lz::SimpleCompressedConcatInputStreamFragment<std::istream>(
					mergeinfo.tmpfileblocks[i][j],
					inputfiles[j].get()
				)
			);
			
		libmaus::lz::SimpleCompressedConcatInputStream<std::istream>::unique_ptr_type tptr(
			new libmaus::lz::SimpleCompressedConcatInputStream<std::istream>(fragments,decfact)
		);
		concfiles[i] = UNIQUE_PTR_MOVE(tptr);
	}
	
	std::vector<uint64_t> tmpoutcnts = mergeinfo.tmpfileblockcntsums;
	order_type BAPC(0);
	::libmaus::autoarray::AutoArray< ::libmaus::bambam::BamAlignment > algns(mergeinfo.tmpfileblocks.size());
	::libmaus::bambam::BamAlignmentHeapComparator<order_type> heapcmp(BAPC,algns.begin());

	::std::priority_queue< 
		uint64_t, 
		std::vector<uint64_t>, 
		::libmaus::bambam::BamAlignmentHeapComparator<order_type> > Q(heapcmp);
	for ( uint64_t i = 0; i < tmpoutcnts.size(); ++i )
		if ( tmpoutcnts[i]-- )
		{
			::libmaus::bambam::BamDecoder::readAlignmentGz(*concfiles[i],algns[i],0 /* no header for validation */,false /* no validation */);
			Q.push(i);
		}
		
	uint64_t outcnt = 0;
	bool const verbose = true;
	mergeinfo.sheader->changeSortOrder(neworder);
	int const level = arginfo.getValue<int>("level",getDefaultLevel());;
	libmaus::bambam::BamParallelWriter writer(std::cout,std::max(1,static_cast<int>(numthreads-1)),*(mergeinfo.sheader),level);
					
	while ( Q.size() )
	{
		uint64_t const t = Q.top(); Q.pop();				
		algns[t].serialise(writer.getStream());
		// writer.writeAlgn(algns[t]);
					
		if ( verbose && (++outcnt % (1024*1024) == 0) )
			std::cerr << "[V] " << outcnt/(1024*1024) << "M" << std::endl;

		if ( tmpoutcnts[t]-- )
		{
			::libmaus::bambam::BamDecoder::readAlignmentGz(*concfiles[t],algns[t],0 /* bamheader */, false /* do not validate */);
			Q.push(t);
		}
	}	
}

struct BamThreadPoolMergeContextBaseConstantsBase
{
	enum bamthreadpoolmergecontextbase_dispatcher_ids
	{
		bamthreadpooldecodecontextbase_dispatcher_id_read = 0,
		bamthreadpooldecodecontextbase_dispatcher_id_decompress = 1,
		bamthreadpooldecodecontextbase_dispatcher_id_process = 2,
		bamthreadpooldecodecontextbase_dispatcher_id_merge = 3,
		bamthreadpooldecodecontextbase_dispatcher_id_compress = 4,
		bamthreadpooldecodecontextbase_dispatcher_id_write = 5
	};
	
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_read = 0;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_decompress = 0;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_process = 0;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_merge = 0;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_compress = 0;
	static unsigned int const bamthreadpooldecodecontextbase_dispatcher_priority_write = 0;
};

template<typename _order_type>
struct BamThreadPoolMergeContextBase;

template<typename _order_type>
struct BamThreadPoolMergeReadPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeReadPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolMergeContextBase<order_type> * contextbase;
	uint64_t blockid;

	BamThreadPoolMergeReadPackage() : contextbase(0) {}
	BamThreadPoolMergeReadPackage(
		BamThreadPoolMergeContextBase<order_type> * rcontextbase,
		uint64_t rblockid
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_read,
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_read,
		0 /* package id */
	), contextbase(rcontextbase), blockid(rblockid)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolMergeReadPackage";
	}	
};

template<typename _order_type>
struct BamThreadPoolMergeDecompressPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeDecompressPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolMergeContextBase<order_type> * contextbase;
	uint64_t blockid;
	uint64_t baseid;
	uint64_t seqid;

	BamThreadPoolMergeDecompressPackage() : contextbase(0) {}
	BamThreadPoolMergeDecompressPackage(
		BamThreadPoolMergeContextBase<order_type> * rcontextbase,
		uint64_t rblockid,
		uint64_t rbaseid,
		uint64_t rseqid
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_decompress,
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_decompress,
		0 /* package id */
	), contextbase(rcontextbase), blockid(rblockid), baseid(rbaseid), seqid(rseqid)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolMergeDecompressPackage";
	}	
};

template<typename _order_type>
struct BamThreadPoolMergeProcessPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeProcessPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolMergeContextBase<order_type> * contextbase;
	uint64_t blockid;
	uint64_t baseid;
	uint64_t seqid;

	BamThreadPoolMergeProcessPackage() : contextbase(0) {}
	BamThreadPoolMergeProcessPackage(
		BamThreadPoolMergeContextBase<order_type> * rcontextbase,
		uint64_t rblockid,
		uint64_t rbaseid,
		uint64_t rseqid
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_process,
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_process,
		0 /* package id */
	), contextbase(rcontextbase), blockid(rblockid), baseid(rbaseid), seqid(rseqid)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolMergeProcessPackage";
	}	
};

template<typename _order_type>
struct BamThreadPoolMergeMergePackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeMergePackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolMergeContextBase<order_type> * contextbase;

	BamThreadPoolMergeMergePackage() : contextbase(0) {}
	BamThreadPoolMergeMergePackage(BamThreadPoolMergeContextBase<order_type> * rcontextbase)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_merge,
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_merge,
		0 /* package id */
	), contextbase(rcontextbase)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolMergeMergePackage";
	}	
};

template<typename _order_type>
struct BamThreadPoolMergeCompressPackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeCompressPackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolMergeContextBase<order_type> * contextbase;
	uint64_t baseid;
	uint64_t seqid;

	BamThreadPoolMergeCompressPackage() : contextbase(0), baseid(0), seqid(0) {}
	BamThreadPoolMergeCompressPackage(
		BamThreadPoolMergeContextBase<order_type> * rcontextbase,
		uint64_t const rbaseid,
		uint64_t const rseqid
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_compress,
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_compress,
		0 /* package id */
	), contextbase(rcontextbase), baseid(rbaseid), seqid(rseqid)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolMergeCompressPackage";
	}	
};

struct BamThreadPoolMergeWritePendingInfo
{
	uint64_t baseid;
	uint64_t seqid;
	libmaus::lz::BgzfDeflateZStreamBaseFlushInfo flushinfo;
	
	BamThreadPoolMergeWritePendingInfo() : baseid(0), seqid(0), flushinfo() {}
	BamThreadPoolMergeWritePendingInfo(
		uint64_t const rbaseid,
		uint64_t const rseqid,
		libmaus::lz::BgzfDeflateZStreamBaseFlushInfo const & rflushinfo)
	: baseid(rbaseid), seqid(rseqid), flushinfo(rflushinfo) {}
};

struct BamThreadPoolMergeWritePendingInfoHeapComparator
{
	bool operator()(
		BamThreadPoolMergeWritePendingInfo const & A,
		BamThreadPoolMergeWritePendingInfo const & B
	)
	{
		return A.seqid > B.seqid;
	}
};

template<typename _order_type>
struct BamThreadPoolMergeWritePackage : public ::libmaus::parallel::SimpleThreadWorkPackage
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeWritePackage<order_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	BamThreadPoolMergeContextBase<order_type> * contextbase;
	BamThreadPoolMergeWritePendingInfo info;

	BamThreadPoolMergeWritePackage() : contextbase(0), info() {}
	BamThreadPoolMergeWritePackage(
		BamThreadPoolMergeContextBase<order_type> * rcontextbase,
		BamThreadPoolMergeWritePendingInfo const & rinfo
	)
	: ::libmaus::parallel::SimpleThreadWorkPackage(
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_priority_write,
		BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_write,
		0 /* package id */
	), contextbase(rcontextbase), info(rinfo)
	{
	}
	virtual char const * getPackageName() const
	{
		return "BamThreadPoolMergeWritePackage";
	}	
};

struct BamThreadPoolMergeProcessBufferInfo
{
	typedef BamThreadPoolMergeProcessBufferInfo this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t blockid;
	uint64_t seqid;
	uint64_t bufferid;
	bool eof;
	BamProcessBuffer * buffer;
	
	BamProcessBuffer::pointer_type * pc;
	BamProcessBuffer::pointer_type * pe;

	BamThreadPoolMergeProcessBufferInfo() 
	: blockid(0), seqid(0), bufferid(0), eof(false) {}
	BamThreadPoolMergeProcessBufferInfo(
		uint64_t rblockid, 
		uint64_t rseqid, 
		uint64_t rbufferid, 
		bool reof,
		BamProcessBuffer * rbuffer
	) 
	: blockid(rblockid), seqid(rseqid), bufferid(rbufferid), eof(reof), buffer(rbuffer),
	  pc(buffer->pc), pe(buffer->pa)
	{
		std::reverse(buffer->pc,buffer->pa);		
	}
	
	template<typename _order_type>
	bool checkOrder(libmaus::bambam::BamHeader const & header) const
	{
		BamProcessBuffer::pointer_type * prev = 0;
		bool ok = true;
		
		::libmaus::bambam::BamFormatAuxiliary aux;
		
		for ( BamProcessBuffer::pointer_type * tc = pc; tc != pe; ++tc )
		{
			if ( prev )
			{
				uint8_t const * ca = buffer->ca + (*prev) + 4;
				uint8_t const * cb = buffer->ca + (*tc) + 4;
				uint32_t const la =
					(static_cast<uint32_t>(ca[0]) << 0) |
					(static_cast<uint32_t>(ca[1]) << 8) |
					(static_cast<uint32_t>(ca[2]) << 16) |
					(static_cast<uint32_t>(ca[3]) << 24);
				uint32_t const lb =
					(static_cast<uint32_t>(cb[0]) << 0) |
					(static_cast<uint32_t>(cb[1]) << 8) |
					(static_cast<uint32_t>(cb[2]) << 16) |
					(static_cast<uint32_t>(cb[3]) << 24);
				
				int const r = _order_type::compareInt(ca,cb);
				
				if ( r > 0 )
				{
					std::cerr << "*order broken* " << std::endl;
					std::cerr << libmaus::bambam::BamAlignmentDecoderBase::formatAlignment(
						ca,la,header,aux) << std::endl;
					std::cerr << libmaus::bambam::BamAlignmentDecoderBase::formatAlignment(
						cb,lb,header,aux) << std::endl;
				}
			}
			
			prev = tc;
		}
		
		return ok;
	}
};

struct BamThreadPoolMergeProcessBufferInfoHeapComparator
{
	bool operator()(BamThreadPoolMergeProcessBufferInfo const & A, BamThreadPoolMergeProcessBufferInfo const & B)
	{
		return A.seqid > B.seqid;
	}
};

template<typename _order_type>
struct BamThreadPoolMergeProcessPackageHeapComparator
{
	typedef _order_type order_type;

	bool operator()(BamThreadPoolMergeProcessPackage<order_type> const & A, BamThreadPoolMergeProcessPackage<order_type> const & B) const
	{
		return A.seqid > B.seqid;
	}
	bool operator()(BamThreadPoolMergeProcessPackage<order_type> const * A, BamThreadPoolMergeProcessPackage<order_type> const * B) const
	{
		return A->seqid > B->seqid;
	}
};


struct BamThreadPoolMergeParserState
{
	enum bam_thread_pool_merge_parser_state
	{
		bam_thread_pool_merge_parser_state_reading_blocklen = 0,
		bam_thread_pool_merge_parser_state_reading_block = 1
	};
	
	bam_thread_pool_merge_parser_state state;
	unsigned int blocklenred;
	uint32_t blocklen;
	uint32_t blockred;
	libmaus::autoarray::AutoArray<uint8_t> gatherbuffer;
	
	BamThreadPoolMergeParserState()
	: state(bam_thread_pool_merge_parser_state_reading_blocklen), blocklenred(0), blocklen(0), gatherbuffer()
	{
	
	}
};

template<typename _order_type>
struct BamThreadPoolMergeHeapComparator
{
	typedef _order_type order_type;
	
	BamThreadPoolMergeHeapComparator() {}
	
	bool operator()(
		std::pair<uint64_t,uint8_t const *> const & A,
		std::pair<uint64_t,uint8_t const *> const & B
	)
	{
		int const r = order_type::compareInt(A.second+4,B.second+4);
		
		if ( r )
			return r > 0;
		else
			return A.first > B.first;
	}
};

struct LockedBitVector
{
	libmaus::parallel::PosixSpinLock lock;
	std::vector<bool> B;
	
	LockedBitVector(uint64_t const n) : lock(), B(n) {}
	~LockedBitVector() {}
	
	bool get(uint64_t const i)
	{	
		bool b;
		
		{
			libmaus::parallel::ScopePosixSpinLock slock(lock);
			b = B[i];
		}
		
		return b;
	}
	
	void set(uint64_t const i, bool const b)
	{
		libmaus::parallel::ScopePosixSpinLock slock(lock);
		B[i] = b;
	}
};

template<typename _order_type>
struct BamThreadPoolMergeContextBase : public BamThreadPoolMergeContextBaseConstantsBase
{
	typedef _order_type order_type;
	typedef BamThreadPoolMergeContextBase<order_type> this_type;

	libmaus::parallel::PosixSpinLock cerrlock;

	libmaus::parallel::SimpleThreadPool & TP;

	MergeInfo const & mergeinfo;
	uint64_t const numblocks;

	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolMergeReadPackage<order_type> > readFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolMergeDecompressPackage<order_type> > decompressFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolMergeProcessPackage<order_type> > processFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolMergeMergePackage<order_type> > mergeFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolMergeCompressPackage<order_type> > compressFreeList;
	libmaus::parallel::SimpleThreadPoolWorkPackageFreeList<BamThreadPoolMergeWritePackage<order_type> > writeFreeList;

	typedef libmaus::lz::SimpleCompressedInputBlockConcatBlock read_block_type;
	typedef read_block_type::unique_ptr_type read_block_ptr_type;
	
	libmaus::autoarray::AutoArray<libmaus::parallel::PosixSpinLock> readQueueLocks;
	libmaus::autoarray::AutoArray<libmaus::parallel::PosixMutex> readLocks;
	libmaus::autoarray::AutoArray<libmaus::parallel::SynchronousQueue<uint64_t>::unique_ptr_type> readBlocksFreeLists;
	libmaus::autoarray::AutoArray<read_block_ptr_type> readBlocks;
	libmaus::autoarray::AutoArray<uint64_t> readSeqIds;
	libmaus::autoarray::AutoArray<libmaus::lz::SimpleCompressedInputBlockConcat::unique_ptr_type> blockReaders;
	libmaus::parallel::SynchronousCounter<uint64_t> blockReadersFinished;
	LockedBitVector blockReadersFinishedVector;
	libmaus::parallel::LockedBool readFinished;
	
	typedef std::priority_queue<
		BamThreadPoolMergeProcessPackage<order_type> *,
		std::vector< BamThreadPoolMergeProcessPackage<order_type> * >,
		BamThreadPoolMergeProcessPackageHeapComparator<order_type>
	> process_pending_type;
	libmaus::parallel::PosixSpinLock processPendingLock;
	libmaus::autoarray::AutoArray<process_pending_type> processPending;
	libmaus::parallel::PosixSpinLock processPendingNextLock;
	libmaus::autoarray::AutoArray<uint64_t> processPendingNext;
	libmaus::autoarray::AutoArray<BamThreadPoolMergeParserState> parserStates;
	libmaus::autoarray::AutoArray<BamProcessBuffer::unique_ptr_type> processBuffers;
	libmaus::autoarray::AutoArray<libmaus::parallel::LockedQueue<uint64_t>::unique_ptr_type> processBuffersFreeLists;
	libmaus::autoarray::AutoArray<
		typename libmaus::parallel::SynchronousQueue < BamThreadPoolMergeProcessPackage<order_type> * >::unique_ptr_type 
	> processStallList;
	libmaus::parallel::SynchronousCounter<uint64_t> blockParsersFinished;
	libmaus::parallel::LockedBool parsingFinished;
	
	libmaus::parallel::PosixSpinLock bufferInfoLock;
	std::vector<
		std::priority_queue<
			BamThreadPoolMergeProcessBufferInfo,
			std::vector<BamThreadPoolMergeProcessBufferInfo>,
			BamThreadPoolMergeProcessBufferInfoHeapComparator
		>
	> bufferInfoPending;
	std::vector<uint64_t> bufferInfoSeq;
	std::vector<uint64_t> bufferInfoNext;
	libmaus::parallel::SynchronousCounter<uint64_t> bufferInfoInitComplete;
	
	libmaus::parallel::PosixSpinLock missingLock;
	libmaus::parallel::LockedBool missingSet;
	uint64_t missing;
	
	libmaus::parallel::PosixSpinLock alPerBlockLock;
	std::vector< uint64_t > alPerBlockFinished;
	uint64_t alFinished;
	uint64_t alTotal;
	std::vector < BamThreadPoolMergeProcessBufferInfo > activeMergeBuffers;

	std::priority_queue<
		std::pair<uint64_t,uint8_t const *>,
		std::vector< 
			std::pair<uint64_t,uint8_t const *> 
		>,
		BamThreadPoolMergeHeapComparator<order_type>
	> mergeQ;
	
	void mergeQueuePush(uint64_t const blockid, uint8_t const * data)
	{
		#if 0
		uint32_t const blocksize =
			(static_cast<uint32_t>(data[0]) << 0 ) |
			(static_cast<uint32_t>(data[1]) << 8 ) |
			(static_cast<uint32_t>(data[2]) << 16) |
			(static_cast<uint32_t>(data[3]) << 24);
		uint8_t const * aldat = data + 4;

		if ( std::string(libmaus::bambam::BamAlignmentDecoderBase::getReadName(aldat)) == "MS6_12707:1:1101:1781:15425" )
		{
			::libmaus::bambam::BamFormatAuxiliary aux;
			std::cerr << "[D]*\t" << libmaus::bambam::BamAlignmentDecoderBase::formatAlignment(
				aldat,blocksize,
				*(mergeinfo.sheader),
				aux
			) << std::endl;
			
			std::cerr << libmaus::util::StackTrace().toString();
		}
		#endif

		mergeQ.push(std::pair<uint64_t,uint8_t const *>(blockid,data));
	}

	uint64_t outputMerged;
	
	
	libmaus::autoarray::AutoArray< libmaus::lz::BgzfDeflateBase::unique_ptr_type > deflateBuffers;
	libmaus::parallel::LockedQueue<uint64_t> deflateBuffersFreeList;
	uint64_t deflateBufferSeq;
	libmaus::parallel::PosixSpinLock deflateBufferNextLock;
	uint64_t deflateBufferNext;
	libmaus::parallel::PosixSpinLock writePendingLock;
	std::priority_queue<
		BamThreadPoolMergeWritePendingInfo,
		std::vector< BamThreadPoolMergeWritePendingInfo >,
		BamThreadPoolMergeWritePendingInfoHeapComparator
	> writePending;
	libmaus::parallel::SynchronousQueue< BamThreadPoolMergeMergePackage<order_type> * > mergeStallList;
	libmaus::autoarray::AutoArray<uint8_t> mergeWritePending;
	uint64_t mergeWritePendingFill;
	
	libmaus::parallel::LockedBool mergeFinished;
	
	std::ostream & out;
	
	#if 0
	libmaus::bambam::BamAlignment thisAlgn;
	uint64_t thisAlgnBlock;
	libmaus::bambam::BamAlignment prevAlgn;
	uint64_t prevAlgnBlock;
	#endif

	libmaus::parallel::SynchronousCounter<uint64_t> mergeWriteIn;
	uint64_t mergeWriteOut;
	libmaus::parallel::PosixSpinLock mergeWriteOutLock;
	
	bool const verbose;
		
	BamThreadPoolMergeContextBase(
		libmaus::parallel::SimpleThreadPool & rTP,
		MergeInfo const & rmergeinfo,
		libmaus::lz::DecompressorObjectFactory & decfact,
		uint64_t const inputBlocksPerBlock,
		uint64_t const processBuffersPerBlock,
		uint64_t const processBuffersSize,
		uint64_t const numDeflateBuffers,
		int const level,
		std::ostream & rout,
		bool const rverbose
	)
	: 
		cerrlock(),
		TP(rTP),
		mergeinfo(rmergeinfo),
		numblocks(mergeinfo.tmpfilenamedblocks.size()),
		readFreeList(),
		decompressFreeList(),
		readQueueLocks(numblocks),
		readLocks(numblocks),
		readBlocksFreeLists(numblocks),
		readBlocks(inputBlocksPerBlock * numblocks),
		readSeqIds(numblocks),
		blockReaders(numblocks),
		blockReadersFinished(0),
		blockReadersFinishedVector(numblocks),
		readFinished(false),
		processPending(numblocks),
		processPendingNextLock(),
		processPendingNext(numblocks),
		parserStates(numblocks),
		processBuffers(numblocks * processBuffersPerBlock),
		processBuffersFreeLists(numblocks),
		processStallList(numblocks),
		blockParsersFinished(0),
		parsingFinished(false),
		bufferInfoLock(),
		bufferInfoPending(numblocks),
		bufferInfoSeq(numblocks,0),
		bufferInfoNext(numblocks,0),
		bufferInfoInitComplete(0),
		missingLock(),
		missingSet(false),
		missing(0),
		alPerBlockFinished(numblocks,0),
		alFinished(0),
		alTotal(std::accumulate(mergeinfo.tmpfileblockcntsums.begin(),mergeinfo.tmpfileblockcntsums.end(),0ull)),
		activeMergeBuffers(numblocks),
		mergeQ(),
		outputMerged(0),
		deflateBuffers(numDeflateBuffers),
		deflateBuffersFreeList(),
		deflateBufferSeq(0),
		deflateBufferNext(0),
		writePendingLock(),
		writePending(),
		mergeStallList(),
		mergeWritePendingFill(0),
		mergeFinished(false),
		out(rout),
		mergeWriteIn(0),
		mergeWriteOut(0),
		mergeWriteOutLock(),
		verbose(rverbose)
	{
		assert ( inputBlocksPerBlock );
	
		for ( uint64_t i = 0; i < readBlocksFreeLists.size(); ++i )
		{
			libmaus::parallel::SynchronousQueue<uint64_t>::unique_ptr_type tptr(
				new libmaus::parallel::SynchronousQueue<uint64_t>()
			);
			readBlocksFreeLists[i] = UNIQUE_PTR_MOVE(tptr);
		}
		for ( uint64_t i = 0; i < readBlocks.size(); ++i )
		{
			read_block_ptr_type tptr(new read_block_type(decfact));
			readBlocks[i] = UNIQUE_PTR_MOVE(tptr);
			readBlocksFreeLists[i/inputBlocksPerBlock]->enque(i);
		}
		for ( uint64_t i = 0; i < blockReaders.size(); ++i )
		{
			libmaus::lz::SimpleCompressedInputBlockConcat::unique_ptr_type tptr(
				new libmaus::lz::SimpleCompressedInputBlockConcat(mergeinfo.tmpfilenamedblocks[i])
			);
			blockReaders[i] = UNIQUE_PTR_MOVE(tptr);
		}
		for ( uint64_t i = 0; i < numblocks; ++i )
		{
			blockReadersFinishedVector.set(i,false);
		}
		for ( uint64_t i = 0; i < processBuffers.size(); ++i )
		{
			BamProcessBuffer::unique_ptr_type tptr(new BamProcessBuffer(processBuffersSize));
			processBuffers[i] = UNIQUE_PTR_MOVE(tptr);
		}
		for ( uint64_t i = 0; i < processBuffersFreeLists.size(); ++i )
		{
			libmaus::parallel::LockedQueue<uint64_t>::unique_ptr_type tptr(new libmaus::parallel::LockedQueue<uint64_t>);
			processBuffersFreeLists[i] = UNIQUE_PTR_MOVE(tptr);
		}
		for ( uint64_t i = 0; i < processBuffers.size(); ++i )
			processBuffersFreeLists[i / processBuffersPerBlock]->push_back(i);

		for ( uint64_t i = 0; i < processStallList.size(); ++i )
		{
			typename libmaus::parallel::SynchronousQueue < BamThreadPoolMergeProcessPackage<order_type> * >::unique_ptr_type 	
				tptr(new libmaus::parallel::SynchronousQueue < BamThreadPoolMergeProcessPackage<order_type> * >);
			processStallList[i] = UNIQUE_PTR_MOVE(tptr);
		}
		for ( uint64_t i = 0; i < deflateBuffers.size(); ++i )
		{
			libmaus::lz::BgzfDeflateBase::unique_ptr_type tptr(
				new libmaus::lz::BgzfDeflateBase(level,true)
			);
			deflateBuffers[i] = UNIQUE_PTR_MOVE(tptr);
			deflateBuffersFreeList.push_back(i);
		}
	}
};

/*
 * block reading
 */
template<typename _order_type>
struct BamThreadPoolMergeReadPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolMergeReadPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolMergeReadPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolMergeReadPackage<order_type> & RP = *dynamic_cast<BamThreadPoolMergeReadPackage<order_type> *>(P);
		BamThreadPoolMergeContextBase<order_type> & contextbase = *(RP.contextbase);
		
		// block id for reading
		uint64_t const blockid = RP.blockid;

		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "read " << blockid << std::endl;
		contextbase.cerrlock.unlock();
		#endif

		bool finishedLoop = (contextbase.blockReadersFinishedVector.get(blockid));
		
		libmaus::parallel::ScopePosixMutex SPM(contextbase.readLocks[blockid]);
		
		while ( !finishedLoop )
		{
			// see if we can get a package from the free list
			uint64_t baseid = 0;
			bool baseidok = false;
			{
				libmaus::parallel::ScopePosixSpinLock llock(contextbase.readQueueLocks[blockid]);
				baseidok = contextbase.readBlocksFreeLists[blockid]->peek(baseid);
			}

			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "read " << blockid << " baseidok=" << baseidok << " list size " << contextbase.readBlocksFreeLists[blockid]->getFillState() << std::endl;
			contextbase.cerrlock.unlock();
			#endif
			
			// if we have space
			if ( baseidok )
			{
				// read package
				contextbase.blockReaders[blockid]->readBlock(*(contextbase.readBlocks[baseid]));
				
				// last one?
				if ( contextbase.readBlocks[baseid]->eof )
				{
					// reading for this block is finished, so quit read loop after this one
					finishedLoop = true;
					// set finished flag for this block
					contextbase.blockReadersFinishedVector.set(blockid,true);
					// increment number of finished readers
					uint64_t const numfinished = ++(contextbase.blockReadersFinished);
					
					// all blocks finished?
					if ( numfinished == contextbase.numblocks )
					{
						// reading finished
						contextbase.readFinished.set(true);
					}
				}
				
				// deque base block
				uint64_t const rbaseid = contextbase.readBlocksFreeLists[blockid]->deque();
				assert ( baseid == rbaseid );
				
				// queue it for decompression
				BamThreadPoolMergeDecompressPackage<order_type> * decpack =
					contextbase.decompressFreeList.getPackage();
				*decpack = BamThreadPoolMergeDecompressPackage<order_type>(
					RP.contextbase,RP.blockid,baseid,contextbase.readSeqIds[blockid]++
				);
				
				tpi.enque(decpack);					
			}
			else
			{
				finishedLoop = true;
			}
		}

		// return package to free list
		contextbase.readFreeList.returnPackage(dynamic_cast<BamThreadPoolMergeReadPackage<order_type> *>(P));
	}
};

/*
 * package decompression
 */
template<typename _order_type>
struct BamThreadPoolMergeDecompressPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;

	virtual ~BamThreadPoolMergeDecompressPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolMergeDecompressPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolMergeDecompressPackage<order_type> & RP = *dynamic_cast<BamThreadPoolMergeDecompressPackage<order_type> *>(P);
		BamThreadPoolMergeContextBase<order_type> & contextbase = *(RP.contextbase);
		
		// block id
		uint64_t const blockid = RP.blockid;
		// base id
		uint64_t const baseid = RP.baseid;
		
		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "decompress " << blockid << "," << baseid << "," << RP.seqid << std::endl;
		contextbase.cerrlock.unlock();
		#endif

		// uncompress the block
		contextbase.readBlocks[baseid]->uncompressBlock();
		
		// generate processing package
		BamThreadPoolMergeProcessPackage<order_type> * procpack = contextbase.processFreeList.getPackage();
		*procpack = BamThreadPoolMergeProcessPackage<order_type>(
			RP.contextbase,RP.blockid,RP.baseid,RP.seqid
		);
		
		{
			// get locks
			libmaus::parallel::ScopePosixSpinLock llprocessPendingLock(contextbase.processPendingLock);
			libmaus::parallel::ScopePosixSpinLock llprocessPendingNextLock(contextbase.processPendingNextLock);
			// put package in pending list
			contextbase.processPending[blockid].push(procpack);
			
			// queue package for dispatcher if top of pending list is next package to be processed
			if ( contextbase.processPending[blockid].top()->seqid == contextbase.processPendingNext[blockid] )
			{
				BamThreadPoolMergeProcessPackage<order_type> * nextprocpack =
					contextbase.processPending[blockid].top();
				contextbase.processPending[blockid].pop();
				
				tpi.enque(nextprocpack);
			}
		}
		
		// return package to free list
		contextbase.decompressFreeList.returnPackage(dynamic_cast<BamThreadPoolMergeDecompressPackage<order_type> *>(P));
	}
};

/*
 * package processing (decompose byte stream into alignments)
 */
template<typename _order_type>
struct BamThreadPoolMergeProcessPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;
	
	uint32_t decodeLEInteger4(uint8_t const * pc)
	{
		return
			(static_cast<uint32_t>(pc[0]) << 0) |
			(static_cast<uint32_t>(pc[1]) << 8) |
			(static_cast<uint32_t>(pc[2]) << 16) |
			(static_cast<uint32_t>(pc[3]) << 24);
	}

	virtual ~BamThreadPoolMergeProcessPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolMergeProcessPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolMergeProcessPackage<order_type> & RP = *dynamic_cast<BamThreadPoolMergeProcessPackage<order_type> *>(P);
		BamThreadPoolMergeContextBase<order_type> & contextbase = *(RP.contextbase);
	
		// block id	
		uint64_t const blockid = RP.blockid;
		// base id
		uint64_t const baseid = RP.baseid;
		
		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "process blockid=" << blockid << ",baseid=" << baseid << ",seqid=" << RP.seqid << std::endl;
		contextbase.cerrlock.unlock();
		#endif

		// current state of the alignment parser
		BamThreadPoolMergeParserState & parserState = contextbase.parserStates[blockid];
		
		// get uncompressed data
		typename BamThreadPoolMergeContextBase<order_type>::read_block_type & datablock = *(contextbase.readBlocks[baseid]);
		uint8_t * const pa = datablock.O.begin();
		uint8_t * const pe = pa + datablock.uncompsize;
		uint8_t * pc = pa;
		
		// try to get process buffer
		uint64_t processBufferId = 0;
		bool stalled = ! (contextbase.processBuffersFreeLists[blockid]->tryDequeFront(processBufferId));
		// we are running if there is a free buffer and the input is not empty
		bool running = (!stalled) && (pc != pe);
		// get pointer to process buffer
		BamProcessBuffer * processBuffer = (!stalled) ? contextbase.processBuffers[processBufferId].get() : 0;
		// ids of finished parse buffers
		std::deque<uint64_t> finishedBuffers;
		
		#if 0
		if ( processBuffer )
		{
			contextbase.cerrlock.lock();
			std::cerr << "process blockid=" << blockid << ",baseid=" << baseid << ",seqid=" << RP.seqid << " buffer fill " 
				<< processBuffer->getFill() 
				<< " running=" << running
				<< " stalled=" << stalled
				<< " pe-pc=" << (pe-pc)
				<< " blocklen=" << parserState.blocklen
				<< " blocklenred=" << parserState.blocklenred
				<< std::endl;
			contextbase.cerrlock.unlock();
		}
		#endif
		
		while ( running )
		{
			switch ( parserState.state )
			{
				// block size not completely read
				case BamThreadPoolMergeParserState::bam_thread_pool_merge_parser_state_reading_blocklen:
				{
					// nothing of block len read an complete alignment block is in this buffer
					if ( 
						(!parserState.blocklenred)
						&&
						(pe-pc) >= 4
						&&
						(pe-pc) >= (4+(parserState.blocklen=decodeLEInteger4(pc)))
					)
					{
						// while complete alignment blocks in this buffer
						while ( 
							running 
							&&
							(pe-pc) >= 4
							&&
							(pe-pc) >= (4+(parserState.blocklen=decodeLEInteger4(pc)))
						)
						{
							// try to put alignment in process buffer
							if ( processBuffer->put(pc+4,parserState.blocklen) )
							{
								// increase pointer if succesfull
								pc += 4+parserState.blocklen;
							}
							else
							{
								#if 0
								contextbase.cerrlock.lock();
								std::cerr << "block " << blockid << " in block push for id " << processBufferId << std::endl;
								contextbase.cerrlock.unlock();
								#endif
								
								// process buffer cannot take another alignment block,
								// mark it as finished
								finishedBuffers.push_back(processBufferId);
								// try to get next process buffer, stall if there is none
								stalled = ! (contextbase.processBuffersFreeLists[blockid]->tryDequeFront(processBufferId));
								// continue if there was a free buffer
								running = running && (!stalled);
								// set buffer pointer
								processBuffer = (!stalled) ? contextbase.processBuffers[processBufferId].get() : 0;
							}
						}
						
						// block len read should be zero
						assert ( ! parserState.blocklenred );
						parserState.blocklen = 0;
					}
					// alignment block is split over two input buffers
					else
					{
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "rest bl blockid=" << blockid << " seqid=" << RP.seqid 
							<< " blocklen=" << parserState.blocklen 
							<< " pe-pc=" << (pe-pc) 
							<< " blocklenred=" << parserState.blocklenred
							<< std::endl;
						contextbase.cerrlock.unlock();
						#endif
						
						// add more block len bytes until we reach the end of the buffer
						while ( pc != pe && parserState.blocklenred < 4 )
						{
							parserState.blocklen |= (static_cast<uint32_t>(*(pc++)) << (8*(parserState.blocklenred++)));
						}
						// block length complete, start to read the data
						if ( parserState.blocklenred == 4 )
						{
							// change state
							parserState.state = BamThreadPoolMergeParserState::bam_thread_pool_merge_parser_state_reading_block;
							// set up for next block
							parserState.blockred = 0;

							#if 0
							contextbase.cerrlock.lock();
							std::cerr << "bl blockid=" << blockid << " seqid=" << RP.seqid << " blocklen=" << parserState.blocklen << std::endl;
							contextbase.cerrlock.unlock();
							#endif
							
							// increase size of gather buffer if necessary
							if ( parserState.gatherbuffer.size() < parserState.blocklen )
								parserState.gatherbuffer = libmaus::autoarray::AutoArray<uint8_t>(parserState.blocklen,false);
						}
					}
					break;
				}
				// read block data
				case BamThreadPoolMergeParserState::bam_thread_pool_merge_parser_state_reading_block:
				{
					assert ( pe >= pc );
					// number of bytes left
					uint32_t const av = pe-pc;
					// number of bytes to be extracted
					uint32_t const tocopy = std::min(
						av,
						static_cast<uint32_t>(parserState.blocklen - parserState.blockred)
					);
					// copy bytes to gather buffer
					std::copy(pc,pc+tocopy,parserState.gatherbuffer.begin()+parserState.blockred);
					
					// increase counters
					parserState.blockred += tocopy;
					pc += tocopy;
					
					// if we have the complete block now
					if ( parserState.blockred == parserState.blocklen )
					{
						// try to put alignment in process buffer
						if ( processBuffer->put(parserState.gatherbuffer.begin(),parserState.blocklen) )
						{
							// if successful then move on to next alignment
							parserState.blocklen = 0;
							parserState.blocklenred = 0;
							parserState.state = BamThreadPoolMergeParserState::bam_thread_pool_merge_parser_state_reading_blocklen;
						}
						// alignment does not fit
						else
						{
							// take back pointers by amount we have just advandec them
							parserState.blockred -= tocopy;
							pc -= tocopy;
						
							#if 0
							contextbase.cerrlock.lock();
							std::cerr << "block " << blockid << " ext block push for id " << processBufferId << std::endl;
							contextbase.cerrlock.unlock();
							#endif
							
							// mark process buffer as finished
							finishedBuffers.push_back(processBufferId);
							// stall if there is no free buffer
							stalled = ! (contextbase.processBuffersFreeLists[blockid]->tryDequeFront(processBufferId));
							// running if not stalled
							running = running && (!stalled);
							// get process buffer pointer
							processBuffer = (!stalled) ? contextbase.processBuffers[processBufferId].get() : 0;
						}
					}
					break;
				}
			}
			
			// running if we have not reached the end of the input
			running = running && (pc != pe);
		}

		// check for eof flag on input
		bool const eof = datablock.eof;

		// if we have reached the end of the file and processed the complete input buffer
		if ( eof && !stalled )
		{
			// insert last finished buffer
			finishedBuffers.push_back(processBufferId);

			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "parsing finished for block " << blockid << std::endl;
			contextbase.cerrlock.unlock();
			#endif
			
			// mark block parsing as finished for this block
			if ( ++contextbase.blockParsersFinished == contextbase.numblocks )
			{
				contextbase.parsingFinished.set(true);
				
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "parsing finished for all blocks." << std::endl;
				contextbase.cerrlock.unlock();
				#endif
			}

		}

		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "process blockid=" << blockid << ",baseid=" << baseid << ",seqid=" << RP.seqid << " finished " << finishedBuffers.size() << std::endl;
		contextbase.cerrlock.unlock();
		#endif
		
		// queue finished blocks
		for ( uint64_t i = 0; i < finishedBuffers.size(); ++i )
		{
			// construct merge process buffer info object
			libmaus::parallel::ScopePosixSpinLock lbufferInfoLock(contextbase.bufferInfoLock);
			BamThreadPoolMergeProcessBufferInfo const info(
				blockid,
				contextbase.bufferInfoSeq[blockid]++,
				finishedBuffers[i],
				datablock.eof && // input buffer is marked as eof
				(!stalled) && // no stall
				(i==finishedBuffers.size()-1), // this is the last queue buffer
				contextbase.processBuffers[finishedBuffers[i]].get()
			);
			// enque in pending list
			contextbase.bufferInfoPending[blockid].push(info);

			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "checking blockid=" 
				<< blockid << " seqid " << info.seqid 
				<< " bufferid " << finishedBuffers[i]
				<< std::endl;
			contextbase.cerrlock.unlock();
			
			info.checkOrder<order_type>(*(contextbase.mergeinfo.sheader));
			#endif
			
			// if this is the first buffer
			if ( info.seqid == 0 )
				// if there is some data for each block then set up the block merging
				if ( (++contextbase.bufferInfoInitComplete) == contextbase.numblocks )
				{
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "ready to run " << std::endl;
					contextbase.cerrlock.unlock();		
					#endif

					// get lock
					libmaus::parallel::ScopePosixSpinLock lalPerBlockLock(contextbase.alPerBlockLock);
					
					for ( uint64_t j = 0; j < contextbase.numblocks; ++j )
					{
						// if block has any alignments
						if ( contextbase.alPerBlockFinished[j] != contextbase.mergeinfo.tmpfileblockcntsums[j] )
						{
							// set active merge input buffer for block j
							contextbase.activeMergeBuffers[j] = contextbase.bufferInfoPending[j].top();

							// get info
							BamThreadPoolMergeProcessBufferInfo & info = contextbase.activeMergeBuffers[j];
							
							// buffer should not be empty
							assert ( info.pc != info.pe );
							
							// put entry in merge heap
							contextbase.mergeQueuePush(j,info.buffer->ca + (*(info.pc++)));
							
							// remove buffer now active from pending list
							contextbase.bufferInfoPending[j].pop();
							
							#if 0
							contextbase.cerrlock.lock();
							std::cerr << "ready to run for " << j 
								<< " " << contextbase.mergeinfo.tmpfileblockcntsums[j] 
								<< " " << info.blockid
								<< " " << info.seqid
								<< " " << info.bufferid
								<< " " << libmaus::bambam::BamAlignmentDecoderBase::getReadName(info.buffer->ca + info.pc[-1] + 4)
								<< std::endl;
							contextbase.cerrlock.unlock();
							#endif
						}
						else
						{
							#if 0
							contextbase.cerrlock.lock();
							std::cerr << "empty to run for " << j << " " << contextbase.mergeinfo.tmpfileblockcntsums[j] << std::endl;
							contextbase.cerrlock.unlock();			
							#endif
						}
					}

					// enque a merge package
					BamThreadPoolMergeMergePackage<order_type> * mergepack =
						contextbase.mergeFreeList.getPackage();
					*mergepack = BamThreadPoolMergeMergePackage<order_type>(RP.contextbase);
					tpi.enque(mergepack);					
				}
		}

		// if we did not stall (input buffer was processed completely)
		if ( ! stalled )
		{
			// reinsert process buffer if it is unfinished
			if ( ! eof )
				contextbase.processBuffersFreeLists[blockid]->push_front(processBufferId);
			
			// return read buffer
			{
				libmaus::parallel::ScopePosixSpinLock llock(contextbase.readQueueLocks[blockid]);
				// check if list is empty
				bool const listempty = contextbase.readBlocksFreeLists[blockid]->empty();

				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "parsed " << baseid << " for " << blockid << " without stall." << std::endl;
				contextbase.cerrlock.unlock();
				#endif
				
				// return buffer
				contextbase.readBlocksFreeLists[blockid]->enque(baseid);

				// queue read
				if ( listempty )
				{
					BamThreadPoolMergeReadPackage<order_type> * ptr = contextbase.readFreeList.getPackage();
					*ptr = BamThreadPoolMergeReadPackage<order_type>(RP.contextbase,blockid);
					tpi.enque(ptr);
				}
			}
			
			libmaus::parallel::ScopePosixSpinLock llprocessPendingLock(contextbase.processPendingLock);
			libmaus::parallel::ScopePosixSpinLock llprocessPendingNextLock(contextbase.processPendingNextLock);
			
			// ready for next block
			contextbase.processPendingNext[blockid] += 1;
			
			if ( 
				(!contextbase.processPending[blockid].empty())
				&&
				contextbase.processPending[blockid].top()->seqid == contextbase.processPendingNext[blockid] 
			)
			{
				BamThreadPoolMergeProcessPackage<order_type> * nextprocpack =
					contextbase.processPending[blockid].top();
				contextbase.processPending[blockid].pop();
				
				tpi.enque(nextprocpack);
			}
		}
		// input buffer was not processed completely
		else
		{
			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "stalling blockid " << blockid << std::endl;
			contextbase.cerrlock.unlock();
			#endif
			
			// move rest of block to front of buffer
			std::memmove(pa,pc,pe-pc);
			datablock.uncompsize = (pe-pc);
			
			// insert buffer in stall list
			contextbase.processStallList[blockid]->enque(dynamic_cast<BamThreadPoolMergeProcessPackage<order_type> *>(P));
			assert ( 
				blockid == dynamic_cast<BamThreadPoolMergeProcessPackage<order_type> *>(P)->blockid
			);
		}

		// check whether this block was missing to continue merging
		{
			// get locks
			libmaus::parallel::ScopePosixSpinLock lbufferInfoLock(contextbase.bufferInfoLock);
			libmaus::parallel::ScopePosixSpinLock lmissingLock(contextbase.missingLock);
			if ( contextbase.missingSet.get() )
			{
				// get id of missing block
				uint64_t const missing = contextbase.missing;

				// check if we have any finished buffers for the missing block
				// and if so whether we have the next required one
				if ( 
					contextbase.bufferInfoPending[missing].size()
					&&
					(
						contextbase.bufferInfoPending[missing].top().seqid ==
						contextbase.bufferInfoNext[missing]
					)
				)
				{
					// set it in active vector
					contextbase.activeMergeBuffers[missing] = contextbase.bufferInfoPending[missing].top();
					// get info object
					BamThreadPoolMergeProcessBufferInfo & info = contextbase.activeMergeBuffers[missing];
					// no longer missing
					contextbase.missingSet.set(false);

					// info should be non empty
					assert ( info.pc != info.pe );
					
					// insert alignment in merge heap
					contextbase.mergeQueuePush(missing,info.buffer->ca + (*(info.pc++)));
					
					// remove info object from pending list
					contextbase.bufferInfoPending[missing].pop();

					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "process() adding missing for " << missing << std::endl;
					contextbase.cerrlock.unlock();
					#endif

					// enque merge package
					BamThreadPoolMergeMergePackage<order_type> * mergepack =
						contextbase.mergeFreeList.getPackage();
					*mergepack = BamThreadPoolMergeMergePackage<order_type>(RP.contextbase);
					tpi.enque(mergepack);				
				}
			}
		}

		if ( ! stalled )
			// return package to free list
			contextbase.processFreeList.returnPackage(dynamic_cast<BamThreadPoolMergeProcessPackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolMergeMergePackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;
	
	// decode a little endian 4 byte number
	static uint32_t decodeLE4(uint8_t const * v)
	{
		return
			(static_cast<uint32_t>(v[0]) << 0) |
			(static_cast<uint32_t>(v[1]) << 8) |
			(static_cast<uint32_t>(v[2]) << 16) |
			(static_cast<uint32_t>(v[3]) << 24) ;
	}

	// enque a compress package
	static void enqueCompress(
		BamThreadPoolMergeMergePackage<order_type> & RP,
		uint64_t const deflatebufferid,
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		BamThreadPoolMergeCompressPackage<order_type> * ccpack =
			RP.contextbase->compressFreeList.getPackage();
		*ccpack = BamThreadPoolMergeCompressPackage<order_type>(
			RP.contextbase,
			deflatebufferid,
			RP.contextbase->deflateBufferSeq++
		);
				
		tpi.enque(ccpack);			
	}

	virtual ~BamThreadPoolMergeMergePackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolMergeMergePackage<order_type> *>(P) != 0 );
		
		BamThreadPoolMergeMergePackage<order_type> & RP = *dynamic_cast<BamThreadPoolMergeMergePackage<order_type> *>(P);
		BamThreadPoolMergeContextBase<order_type> & contextbase = *(RP.contextbase);

		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "MERGE, mergeQ.size()=" << contextbase.mergeQ.size() << std::endl;
		contextbase.cerrlock.unlock();
		#endif

		// get the alignments per block lock
		libmaus::parallel::ScopePosixSpinLock lalPerBlockLock(contextbase.alPerBlockLock);

		// try to get a deflate buffer id
		uint64_t deflatebufferid = 0;
		// stall if there is no free buffer
		bool stalled = !contextbase.deflateBuffersFreeList.tryDequeFront(deflatebufferid);
		// get pointer to buffer
		libmaus::lz::BgzfDeflateBase * deflateBuffer =
			stalled ? 0 : contextbase.deflateBuffers[deflatebufferid].get();
		// either there is no buffer or the one we have is not yet full
		assert ( deflateBuffer == 0 || deflateBuffer->pc != deflateBuffer->pe );
		// list of buffers to be compressed
		// std::deque<uint64_t> compressPendingBuffer;

		// handle pending write data from previous run
		while ( (!stalled) && (contextbase.mergeWritePendingFill) )
		{
			// number of free bytes in output deflate buffer
			uint64_t const spaceinbuffer = (deflateBuffer->pe - deflateBuffer->pc);	
			// number of bytes to be written
			uint64_t const bufferwrite = std::min(spaceinbuffer,contextbase.mergeWritePendingFill);
			
			// copy the bytes
			std::copy(
				contextbase.mergeWritePending.begin(),
				contextbase.mergeWritePending.begin()+bufferwrite,
				deflateBuffer->pc
			);
			
			// increase buffer pointer in output deflate buffer
			deflateBuffer->pc += bufferwrite;
			
			if ( bufferwrite != contextbase.mergeWritePendingFill )
				// move rest to front of pending buffer
				std::memmove(
					contextbase.mergeWritePending.begin(),
					contextbase.mergeWritePending.begin() + bufferwrite,
					contextbase.mergeWritePendingFill - bufferwrite
				);
				
			contextbase.mergeWritePendingFill -= bufferwrite;
			
			// check whether buffer is full now
			if ( deflateBuffer->pc == deflateBuffer->pe )
			{
				// pass to compression
				contextbase.mergeWriteIn += 1;
				enqueCompress(RP,deflatebufferid,tpi);
				// compressPendingBuffer.push_back(deflatebufferid);
				// try to get next buffer
				stalled = !contextbase.deflateBuffersFreeList.tryDequeFront(deflatebufferid);
				// set new buffer pointer
				deflateBuffer = stalled ? 0 : contextbase.deflateBuffers[deflatebufferid].get();
			}
		}

		bool running = (contextbase.mergeQ.size() != 0) && (!stalled);
		bool finishflag = false;
		
		while ( running )
		{
			// get merge queue top entry
			std::pair<uint64_t,uint8_t const *> const P = contextbase.mergeQ.top();
			// pop entry
			contextbase.mergeQ.pop();
			
			// get block length
			uint32_t const blocklen = decodeLE4(P.second);
			// block start pointer
			uint8_t const * pa = P.second;
			// current output pointer
			uint8_t const * pc = pa;
			// block end pointer
			uint8_t const * pe = pa + blocklen + 4;
			
			#if 0
			{
				uint8_t const * alp = pa + sizeof(uint32_t);
				::libmaus::bambam::BamFormatAuxiliary aux;
				libmaus::bambam::BamAlignmentDecoderBase::formatAlignment(
					alp,blocklen,
					*(contextbase.mergeinfo.sheader),
					aux
				);
				
				contextbase.cerrlock.lock();
				std::cerr << "[D]\t" << libmaus::bambam::BamAlignmentDecoderBase::formatAlignment(
						alp,blocklen,
						*(contextbase.mergeinfo.sheader),
						aux
					)
					<< std::endl;
				contextbase.cerrlock.unlock();		
			}
			#endif
			
			#if 0
			if ( contextbase.thisAlgn.D.size() < blocklen )
				contextbase.thisAlgn.D = libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
			contextbase.thisAlgn.blocksize = blocklen;
			contextbase.thisAlgnBlock = P.first;
			std::copy(pa+4,pe,contextbase.thisAlgn.D.begin());
			
			if ( contextbase.prevAlgn.blocksize )
			{
				if ( 
					order_type::compareInt(contextbase.prevAlgn.D.begin(),contextbase.thisAlgn.D.begin()) > 0 
				)
				{
					contextbase.cerrlock.lock();
					std::cerr << "order wrong " << contextbase.prevAlgnBlock << " " << contextbase.thisAlgnBlock << std::endl;
					std::cerr << contextbase.prevAlgn.formatAlignment(*(contextbase.mergeinfo.sheader)) << std::endl;
					std::cerr << contextbase.thisAlgn.formatAlignment(*(contextbase.mergeinfo.sheader)) << std::endl;
					contextbase.cerrlock.unlock();
				}
			}
			
			contextbase.prevAlgn.copyFrom(contextbase.thisAlgn);
			contextbase.prevAlgnBlock = contextbase.thisAlgnBlock;
			#endif
			
			// while alignment block data remaining
			while ( pc != pe )
			{
				// number of bytes still pending
				uint64_t const writepending = pe-pc;
				// amount of space in buffer
				uint64_t const bufferspace = deflateBuffer->pe - deflateBuffer->pc;

				if ( bufferspace )
				{
					// copy as much as possible to output buffer
					uint64_t const tocopy = std::min(writepending,bufferspace);
					// copy
					std::copy(pc,pc+tocopy,deflateBuffer->pc);

					// advance pointers
					pc += tocopy;
					deflateBuffer->pc += tocopy;
				}
				
				// if buffer is full
				if ( deflateBuffer->pc == deflateBuffer->pe )
				{
					// pass to compression
					contextbase.mergeWriteIn += 1;
					enqueCompress(RP,deflatebufferid,tpi);
					// compressPendingBuffer.push_back(deflatebufferid);
					// try to get next buffer
					stalled = !contextbase.deflateBuffersFreeList.tryDequeFront(deflatebufferid);
					// set new buffer pointer
					deflateBuffer = stalled ? 0 : contextbase.deflateBuffers[deflatebufferid].get();
					
					// no more buffers?
					if ( stalled )
					{
						// copy pending data to pending buffer
						uint64_t const writepending = pe-pc;
						
						// check buffer size
						if ( contextbase.mergeWritePending.size() < writepending )
							contextbase.mergeWritePending = libmaus::autoarray::AutoArray<uint8_t>(writepending,false);
						
						// copy data
						std::copy(pc,pe,contextbase.mergeWritePending.begin());
						// set number of bytes in buffer
						contextbase.mergeWritePendingFill = writepending;

						// quit loop, there is no more buffer space available					
						running = false;
						// quit inner loop
						pc += writepending;
						// sanity check
						assert ( pc == pe );
					}
				}
			}
			
			if ( ((++contextbase.outputMerged) & (1024*1024-1)) == 0 )
			{
				if ( contextbase.verbose )
				{
					contextbase.cerrlock.lock();
					std::cerr << "[V] " << contextbase.outputMerged << std::endl;
					contextbase.cerrlock.unlock();
				}
			}
			
			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "NAME " << libmaus::bambam::BamAlignmentDecoderBase::getReadName(P.second+4) 
				<< " " << P.first
				<< " " << libmaus::bambam::BamAlignmentDecoderBase::getRefID(P.second+4) 
				<< " " << libmaus::bambam::BamAlignmentDecoderBase::getPos(P.second+4) 
				<< std::endl;
			contextbase.cerrlock.unlock();
			#endif
			
			// increase alignments processed for this block
			contextbase.alPerBlockFinished[P.first]++;
			
			// if there are more alignments for this block
			if ( 
				contextbase.alPerBlockFinished[P.first] != 
				contextbase.mergeinfo.tmpfileblockcntsums[P.first] 
			)
			{
				// get merge buffer info
				BamThreadPoolMergeProcessBufferInfo & info = contextbase.activeMergeBuffers[P.first];
				
				// if there are more alignments in the process buffer
				if ( info.pc != info.pe )
				{
					// put alignment in heap
					contextbase.mergeQueuePush(P.first,info.buffer->ca + (*(info.pc++)));
				}
				// otherwise
				else
				{
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "out of data for block " << P.first << std::endl;
					contextbase.cerrlock.unlock();
					#endif
				
					// reset buffer
					info.buffer->reset();
					//
					uint64_t const freebufferid = info.bufferid;

					libmaus::parallel::ScopePosixSpinLock lbufferInfoLock(contextbase.bufferInfoLock);
					libmaus::parallel::ScopePosixSpinLock lmissingLock(contextbase.missingLock);

					// switch to next
					contextbase.bufferInfoNext[P.first]++;
					
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "pending size " << contextbase.bufferInfoPending[P.first].size() << std::endl;
					if ( contextbase.bufferInfoPending[P.first].size() )
						std::cerr << "top " << contextbase.bufferInfoPending[P.first].top().seqid << std::endl;
					std::cerr << "expecting " << contextbase.bufferInfoNext[P.first] << std::endl;
					contextbase.cerrlock.unlock();
					#endif
					
					#if 0
					assert ( 
						contextbase.bufferInfoPending[P.first].size()
						&&
						contextbase.bufferInfoNext[P.first] == contextbase.bufferInfoPending[P.first].top().seqid
					);
					
					contextbase.bufferInfoPending[P.first].pop();
					#endif
					
					// check if next block is already in queue
					if ( 
						contextbase.bufferInfoPending[P.first].size()
						&&
						(
							contextbase.bufferInfoPending[P.first].top().seqid ==
							contextbase.bufferInfoNext[P.first]
						)
					)
					{
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "next block already waiting." << std::endl;
						contextbase.cerrlock.unlock();
						#endif
					
						// mark next buffer as active
						contextbase.activeMergeBuffers[P.first] =
							contextbase.bufferInfoPending[P.first].top();
						// remove if from pending queue
						contextbase.bufferInfoPending[P.first].pop();

						// get info
						BamThreadPoolMergeProcessBufferInfo & info = contextbase.activeMergeBuffers[P.first];

						// put next alignment in merge heap
						contextbase.mergeQueuePush(P.first,info.buffer->ca + (*(info.pc++)));

						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "new queue size is " << contextbase.mergeQ.size() << std::endl;
						contextbase.cerrlock.unlock();
						#endif
					}
					// otherwise note that we are waiting for it
					else
					{
						#if 0
						std::cerr << "marking " << P.first << " as missing." << std::endl;
						std::cerr << "contextbase.bufferInfoPending[P.first].size()=" << contextbase.bufferInfoPending[P.first].size() <<std::endl;
						std::cerr << "expecting " << contextbase.bufferInfoNext[P.first] << std::endl;
						if ( contextbase.bufferInfoPending[P.first].size() )
							std::cerr << "top is " << contextbase.bufferInfoPending[P.first].top().seqid << std::endl;
						#endif
						
						contextbase.missing = P.first;
						contextbase.missingSet.set(true);

						// we are out of data for this block
						running = false;
					}	
			
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "adding " << freebufferid << " to freelist of " << P.first << std::endl;
					contextbase.cerrlock.unlock();					
					#endif

					// add it to the free list
					contextbase.processBuffersFreeLists[P.first]->push_back(freebufferid);

					// reinsert stalled package if there is any
					BamThreadPoolMergeProcessPackage<order_type> * stalled = 0;
					bool ok = contextbase.processStallList[P.first]->trydeque(stalled);
					if ( ok )
					{
						assert ( stalled->blockid == P.first );
						
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "reinserting stalled block " << stalled->blockid << std::endl;
						contextbase.cerrlock.unlock();
						#endif
						
						tpi.enque(stalled);
					}
					else
					{
						#if 0
						contextbase.cerrlock.lock();
						std::cerr << "stall list is empty." << std::endl;
						contextbase.cerrlock.unlock();
						#endif

					}
				}
			}
			else
			{
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "all data for " << P.first << " seen." << std::endl;
				contextbase.cerrlock.unlock();
				#endif
				
				if ( contextbase.mergeQ.size() == 0 )
				{
					#if 0
					contextbase.cerrlock.lock();
					std::cerr << "queue is empty, terminating loop." << std::endl;
					contextbase.cerrlock.unlock();
					#endif
					
					running = false;
					
					finishflag = true;
					
					// tpi.terminate();
				}
			}
		}

		if ( finishflag )
			contextbase.mergeFinished.set(true);

		if ( ! stalled )
		{
			// if we are finished, enque the final buffer for compression
			if ( finishflag )
			{
				contextbase.mergeWriteIn += 1;
				enqueCompress(RP,deflatebufferid,tpi);
				// compressPendingBuffer.push_back(deflatebufferid);
				
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "should write last buffer" << std::endl;
				contextbase.cerrlock.unlock();					
				#endif
			}
			// otherwise
			else
			{
				#if 0
				contextbase.cerrlock.lock();
				std::cerr << "no stall, reinsering buffer" << std::endl;
				contextbase.cerrlock.unlock();	
				#endif
			
				// put unfinished deflate buffer back into the free list
				contextbase.deflateBuffersFreeList.push_front(deflatebufferid);
			}
		}
		
		#if 0
		for ( uint64_t i = 0; i < compressPendingBuffer.size(); ++i )
			enqueCompress(RP,compressPendingBuffer[i],tpi);
		#endif

		if ( stalled )
			contextbase.mergeStallList.enque(dynamic_cast<BamThreadPoolMergeMergePackage<order_type> *>(P));
		else
			contextbase.mergeFreeList.returnPackage(dynamic_cast<BamThreadPoolMergeMergePackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolMergeCompressPackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;
	
	virtual ~BamThreadPoolMergeCompressPackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolMergeCompressPackage<order_type> *>(P) != 0 );
		
		BamThreadPoolMergeCompressPackage<order_type> & RP = *dynamic_cast<BamThreadPoolMergeCompressPackage<order_type> *>(P);
		BamThreadPoolMergeContextBase<order_type> & contextbase = *(RP.contextbase);

		// compress
		libmaus::lz::BgzfDeflateZStreamBaseFlushInfo const FI = 
			contextbase.deflateBuffers[RP.baseid]->flush(true);

		// pending info
		BamThreadPoolMergeWritePendingInfo const pendinf(RP.baseid,RP.seqid,FI);
		
		// queue write
		{
		libmaus::parallel::ScopePosixSpinLock swritePendingLock(contextbase.writePendingLock);
		contextbase.writePending.push(pendinf);
		}
		
		// check if next write in line is available
		{
			libmaus::parallel::ScopePosixSpinLock ldeflateBufferNextLock(contextbase.deflateBufferNextLock);
			libmaus::parallel::ScopePosixSpinLock swritePendingLock(contextbase.writePendingLock);
			if ( 
				contextbase.writePending.size()
				&&
				contextbase.deflateBufferNext == contextbase.writePending.top().seqid 
			)
			{
				BamThreadPoolMergeWritePackage<order_type> * wpack = contextbase.writeFreeList.getPackage();
				*wpack = BamThreadPoolMergeWritePackage<order_type>(
					RP.contextbase,
					contextbase.writePending.top()
				);
				
				contextbase.writePending.pop();
				
				tpi.enque(wpack);
			}
		}
				
		contextbase.compressFreeList.returnPackage(dynamic_cast<BamThreadPoolMergeCompressPackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolMergeWritePackageDispatcher : public libmaus::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef _order_type order_type;
	
	virtual ~BamThreadPoolMergeWritePackageDispatcher() {}
	virtual void dispatch(
		libmaus::parallel::SimpleThreadWorkPackage * P, 
		libmaus::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & tpi
	)
	{
		assert ( dynamic_cast<BamThreadPoolMergeWritePackage<order_type> *>(P) != 0 );
		
		BamThreadPoolMergeWritePackage<order_type> & RP = *dynamic_cast<BamThreadPoolMergeWritePackage<order_type> *>(P);
		BamThreadPoolMergeContextBase<order_type> & contextbase = *(RP.contextbase);

		BamThreadPoolMergeWritePendingInfo const & pendinf = RP.info;

		#if 0
		contextbase.cerrlock.lock();
		std::cerr << "write " << RP.info.seqid << std::endl;
		contextbase.cerrlock.unlock();
		#endif

		uint64_t const compsize = pendinf.flushinfo.getCompressedSize();
		contextbase.out.write(
			reinterpret_cast<char const *>(contextbase.deflateBuffers[pendinf.baseid]->outbuf.begin()),
			compsize
		);

		contextbase.mergeWriteOutLock.lock();
		contextbase.mergeWriteOut += 1;
		contextbase.mergeWriteOutLock.unlock();
		
		{
			libmaus::parallel::ScopePosixSpinLock ldeflateBufferNextLock(contextbase.deflateBufferNextLock);
			libmaus::parallel::ScopePosixSpinLock swritePendingLock(contextbase.writePendingLock);
			
			contextbase.deflateBufferNext += 1;

			#if 0
			contextbase.cerrlock.lock();
			std::cerr << "pending in write " << contextbase.writePending.size() << std::endl;
			if ( contextbase.writePending.size() )
				std::cerr << "top " << contextbase.writePending.top().seqid << " next " << contextbase.deflateBufferNext << std::endl;
			contextbase.cerrlock.unlock();
			#endif
			
			if ( 
				contextbase.writePending.size()
				&&
				contextbase.deflateBufferNext == contextbase.writePending.top().seqid 
			)
			{
				BamThreadPoolMergeWritePackage<order_type> * wpack = contextbase.writeFreeList.getPackage();
				*wpack = BamThreadPoolMergeWritePackage<order_type>(
					RP.contextbase,
					contextbase.writePending.top()
				);
				
				contextbase.writePending.pop();
				
				tpi.enque(wpack);
			}
		}

		// mark buffer as free		
		contextbase.deflateBuffersFreeList.push_back(pendinf.baseid);
		
		// unstall merging if stalled
		BamThreadPoolMergeMergePackage<order_type> * mergepack = 0;
		if ( contextbase.mergeStallList.trydeque(mergepack) )
			tpi.enque(mergepack);
		
		{
			libmaus::parallel::ScopePosixSpinLock lmergeWriteOutLock(contextbase.mergeWriteOutLock);
			
			if (
				contextbase.mergeFinished.get() &&
				(contextbase.mergeWriteIn.get() == contextbase.mergeWriteOut)
			)
			{
				if ( contextbase.verbose )
				{
					contextbase.cerrlock.lock();
					std::cerr << "[V] writing complete" << std::endl;
					contextbase.cerrlock.unlock();
				}
				
				tpi.terminate();
			}
		}
		
		contextbase.writeFreeList.returnPackage(dynamic_cast<BamThreadPoolMergeWritePackage<order_type> *>(P));
	}
};

template<typename _order_type>
struct BamThreadPoolMergeContext : public BamThreadPoolMergeContextBase<_order_type>
{
	typedef _order_type order_type;

	BamThreadPoolMergeContext(
		libmaus::parallel::SimpleThreadPool & TP,
		MergeInfo const & mergeinfo,
		libmaus::lz::DecompressorObjectFactory	& decfact,
		uint64_t const inputBlocksPerBlock,
		uint64_t const processBuffersPerBlock,
		uint64_t const processBuffersSize,
		uint64_t const numDeflateBuffers,
		int const level,
		std::ostream & out,
		bool const verbose
	)
	: BamThreadPoolMergeContextBase<order_type>(
		TP,mergeinfo,decfact,inputBlocksPerBlock,
		processBuffersPerBlock,processBuffersSize,numDeflateBuffers,level,out,verbose)
	{
	
	}
	
	void startup()
	{
		for ( uint64_t i = 0; i < BamThreadPoolMergeContextBase<_order_type>::numblocks; ++i )
		{
			BamThreadPoolMergeReadPackage<order_type> * ptr = BamThreadPoolMergeContextBase<_order_type>::readFreeList.getPackage();
			*ptr = BamThreadPoolMergeReadPackage<order_type>(this,i);
			BamThreadPoolMergeContextBase<_order_type>::TP.enque(ptr);
		}
	}
};

template<typename _order_type>
void checkSortedBlocks(MergeInfo const & mergeinfo, libmaus::lz::DecompressorObjectFactory & decfact)
{
	libmaus::autoarray::AutoArray<libmaus::aio::CheckedInputStream::unique_ptr_type> inputfiles(mergeinfo.tmpfilenames.size());
	for ( uint64_t i = 0; i < inputfiles.size(); ++i )
	{
		libmaus::aio::CheckedInputStream::unique_ptr_type tptr(new libmaus::aio::CheckedInputStream(mergeinfo.tmpfilenames[i]));
		inputfiles[i] = UNIQUE_PTR_MOVE(tptr);
	}
	
	libmaus::autoarray::AutoArray<libmaus::lz::SimpleCompressedConcatInputStream<std::istream>::unique_ptr_type> concfiles(mergeinfo.tmpfileblocks.size());
	for ( uint64_t i = 0; i < mergeinfo.tmpfileblocks.size(); ++i )
	{
		std::vector<libmaus::lz::SimpleCompressedConcatInputStreamFragment<std::istream> > fragments;
		
		for ( uint64_t j = 0; j < mergeinfo.tmpfileblocks[i].size(); ++j )
			fragments.push_back(
				libmaus::lz::SimpleCompressedConcatInputStreamFragment<std::istream>(
					mergeinfo.tmpfileblocks[i][j],
					inputfiles[j].get()
				)
			);
			
		libmaus::lz::SimpleCompressedConcatInputStream<std::istream>::unique_ptr_type tptr(
			new libmaus::lz::SimpleCompressedConcatInputStream<std::istream>(fragments,decfact)
		);
		concfiles[i] = UNIQUE_PTR_MOVE(tptr);
	}

	for ( uint64_t i = 0; i < mergeinfo.tmpfileblocks.size(); ++i )
	{
		libmaus::bambam::BamAlignment prev, cur;

		for ( uint64_t j = 0; j < mergeinfo.tmpfileblockcntsums[i]; ++j )
		{
			::libmaus::bambam::BamDecoder::readAlignmentGz(*concfiles[i],cur,0 /* no header for validation */,false /* no validation */);

			if ( prev.blocksize )
			{
				int const r = _order_type::compareInt(
					prev.D.begin(),
					cur.D.begin()
				);
				
				if ( r > 0 )
				{
					std::cerr << "NO" << std::endl;
				}
			}

			prev.copyFrom(cur);
		}
	}
}

template<typename _order_type>
void checkSortedBlocksByBlocks(MergeInfo const & mergeinfo, libmaus::lz::DecompressorObjectFactory & decfact)
{
	for ( uint64_t i = 0; i < mergeinfo.tmpfilenamedblocks.size(); ++i )
	{
		std::cerr << "block " << (i) << "/" << mergeinfo.tmpfilenamedblocks.size() << std::endl;
		
		libmaus::lz::SimpleCompressedInputBlockConcat conc(mergeinfo.tmpfilenamedblocks[i]);
		libmaus::lz::SimpleCompressedInputBlockConcatBlock block(decfact);
		std::ostringstream ostr;
		
		while ( conc.readBlock(block) )
		{
			block.uncompressBlock();
			ostr.write(
				reinterpret_cast<char const *>(block.O.begin()),
				block.uncompsize
			);
		
			if ( block.eof )
				break;
		}
		
		std::string const sdata = ostr.str();
		char const * cdata = sdata.c_str();
		uint8_t const * udata = reinterpret_cast<uint8_t const *>(cdata);
		uint64_t datalen = sdata.size();
		uint8_t const * udataend = udata + datalen;
		
		uint8_t const * prev = 0;
		uint32_t prevlen = 0;
		uint64_t alcnt = 0;
		uint64_t failcnt = 0;
		int64_t firstfail = -1;
		
		while ( udata != udataend )
		{
			uint32_t const len =
				(static_cast<uint32_t>(udata[0]) << 0) |
				(static_cast<uint32_t>(udata[1]) << 8) |
				(static_cast<uint32_t>(udata[2]) << 16) |
				(static_cast<uint32_t>(udata[3]) << 24);

			if ( prev )
			{
				int const r = _order_type::compareInt(
					prev+4,
					udata+4
				);

				if ( r > 0 )
				{
					if ( firstfail < 0 )
						firstfail = alcnt;
					// std::cerr << "NO" << std::endl;
					++failcnt;						
				}
			}
			alcnt += 1;
				
			prevlen = len;
			prev = udata;
			udata += len+4;
		}

		std::cerr << "block " << (i) << "/" << mergeinfo.tmpfilenamedblocks[i].size() << " done, alcnt=" << alcnt << " failcnt=" << failcnt << " first fail=" << firstfail << std::endl;
	}	
}

static bool startsWith(std::string const & a, std::string const & b)
{
	return	a.size() >= b.size() && a.substr(0,b.size()) == b;
}

static std::string getDefaultTempComp()
{
	return "zlib:-1";
}

template<typename _order_type>
int bamparsortTemplate(libmaus::util::ArgInfo const & arginfo, std::string const & neworder)
{
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	
	libmaus::lz::CompressorObjectFactory::unique_ptr_type PcompressorFactory;
	libmaus::lz::DecompressorObjectFactory::unique_ptr_type PdecompressorFactory;
	
	std::string const tempcomp = arginfo.getValue<std::string>("tempcomp",getDefaultTempComp());

	if ( startsWith(tempcomp,"snappy") )
	{
		libmaus::lz::CompressorObjectFactory::unique_ptr_type TcompressorFactory(new libmaus::lz::SnappyCompressorObjectFactory());
		PcompressorFactory = UNIQUE_PTR_MOVE(TcompressorFactory);
		libmaus::lz::DecompressorObjectFactory::unique_ptr_type TdecompressorFactory(new libmaus::lz::SnappyDecompressorObjectFactory);
		PdecompressorFactory = UNIQUE_PTR_MOVE(TdecompressorFactory);
	}
	else if ( startsWith(tempcomp,"zlib:") )
	{
		std::string const levelstr = tempcomp.substr(strlen("zlib:"));
		std::istringstream levelistr(levelstr);
		int64_t templevel = -1;
		levelistr >> templevel;

		if ( ! levelistr )
		{
			libmaus::exception::LibMausException lme;
			lme.getStream() << "Cannot parse compression setting in " << tempcomp << std::endl;
			lme.finish();
			throw lme;
		}
		
		templevel = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(templevel);
	
		libmaus::lz::CompressorObjectFactory::unique_ptr_type TcompressorFactory(new libmaus::lz::ZlibCompressorObjectFactory(templevel));
		PcompressorFactory = UNIQUE_PTR_MOVE(TcompressorFactory);
		libmaus::lz::DecompressorObjectFactory::unique_ptr_type TdecompressorFactory(new libmaus::lz::ZlibDecompressorObjectFactory);
		PdecompressorFactory = UNIQUE_PTR_MOVE(TdecompressorFactory);
	}
	else
	{
		libmaus::exception::LibMausException lme;
		lme.getStream() << "Unsupported temp file compression scheme in " << tempcomp << std::endl;
		lme.finish();
		throw lme;	
	}
		
	libmaus::lz::CompressorObjectFactory & compressorFactory = *PcompressorFactory;
	libmaus::lz::DecompressorObjectFactory & decompressorFactory = *PdecompressorFactory;

	std::cerr << "[V] using " << compressorFactory.getDescription() << " to compress temporary files." << std::endl;
	
	typedef _order_type order_type;
	MergeInfo const mergeinfo = produceSortedBlocks<order_type>(arginfo,compressorFactory);
	

	#if 0
	mergeSortedBlocksNoThreadPool<order_type>(arginfo,mergeinfo,neworder,decompressorFactory);
	#else
	
	// write bam header
	{
		std::ostringstream headerostr;
		libmaus::lz::BgzfOutputStream writer(headerostr);
		mergeinfo.sheader->changeSortOrder(neworder);

		std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
			mergeinfo.sheader->text,
			"bamparsort", // ID
			"bamparsort", // PN
			arginfo.commandline, // CL
			::libmaus::bambam::ProgramHeaderLineSet(mergeinfo.sheader->text).getLastIdInChain(), // PP
			std::string(PACKAGE_VERSION) // VN			
		);
		libmaus::bambam::BamHeader const uphead(upheadtext);

		uphead.serialise(writer);

		writer.flush();
		std::string const & sheader = headerostr.str();
		std::cout.write(sheader.c_str(),sheader.size());
	}
	
	// the merging step cannot handle an "empty" file, so skip merging if there are no alignments in the file
	if ( ! mergeinfo.empty() )
	{
		std::string const tmpfilenamebase = 
			arginfo.getUnparsedValue("tmpfileprefix",arginfo.getDefaultTmpFileName());

		uint64_t const numthreads = arginfo.getValue<unsigned int>("numthreads", libmaus::parallel::NumCpus::getNumLogicalProcessors());

		// thread pool	
		libmaus::parallel::SimpleThreadPool TP(numthreads);

		// package dispatchers
		BamThreadPoolMergeReadPackageDispatcher<order_type> readdispatcher;
		TP.registerDispatcher(BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_read,&readdispatcher);
		BamThreadPoolMergeDecompressPackageDispatcher<order_type> decompressdispatcher;
		TP.registerDispatcher(BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_decompress,&decompressdispatcher);
		BamThreadPoolMergeProcessPackageDispatcher<order_type> processdispatcher;
		TP.registerDispatcher(BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_process,&processdispatcher);
		BamThreadPoolMergeMergePackageDispatcher<order_type> mergedispatcher;
		TP.registerDispatcher(BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_merge,&mergedispatcher);
		BamThreadPoolMergeCompressPackageDispatcher<order_type> compressdispatcher;
		TP.registerDispatcher(BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_compress,&compressdispatcher);
		BamThreadPoolMergeWritePackageDispatcher<order_type> writedispatcher;
		TP.registerDispatcher(BamThreadPoolMergeContextBaseConstantsBase::bamthreadpooldecodecontextbase_dispatcher_id_write,&writedispatcher);

		uint64_t const numblocks = mergeinfo.tmpfilenamedblocks.size();
		uint64_t const mem = arginfo.getValueUnsignedNumeric<uint64_t>("mem",getDefaultMem());
		int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));;
		// input blocks, output blocks, process buffers
		uint64_t const inputBlockMemory = 0.1 * mem;
		uint64_t const inputMemoryPerBlock = std::max ( static_cast<uint64_t>(1), (inputBlockMemory + numblocks-1)/numblocks);
		uint64_t const inputBlocksPerBlock = std::max ( static_cast<uint64_t>(1), (inputMemoryPerBlock + (2*64*1024-1))/(2*64*1024));
		
		uint64_t const processBlockMemory = 0.6 * mem;
		uint64_t const processBuffersPerBlock = 16;
		uint64_t const processBuffersSize = ( processBlockMemory + (processBuffersPerBlock*numblocks) - 1 ) / (processBuffersPerBlock*numblocks);
		
		uint64_t const deflateMemory = 0.1 * mem;
		uint64_t const numDeflateBlocks = std::max ( static_cast<uint64_t>(1), (deflateMemory + (2*64*1024-1)) / (2*64*1024) );

		if ( verbose )
		{
			std::cerr << "[V] input blocks per block: " << inputBlocksPerBlock << std::endl;
			std::cerr << "[V] processBuffersPerBlock: " << processBuffersPerBlock << std::endl;
			std::cerr << "[V] processBuffersSize:     " << processBuffersSize << std::endl;
			std::cerr << "[V] numDeflateBlocks:       " << numDeflateBlocks << std::endl;
		}

		BamThreadPoolMergeContext<order_type> mergecontext(TP,mergeinfo,decompressorFactory,
			inputBlocksPerBlock,
			processBuffersPerBlock,processBuffersSize,
			numDeflateBlocks,
			level,
			std::cout,
			verbose
		);
		
		mergecontext.startup();
		TP.join();
	}

	// write bam footer
	{
		std::ostringstream footerostr;
		libmaus::lz::BgzfOutputStream writer(footerostr);
		writer.flush();
		writer.addEOFBlock();
		std::string const & sfooter = footerostr.str();
		std::cout.write(sfooter.c_str(),sfooter.size());
	}

	
	std::cout.flush();
	#endif
	
	return EXIT_SUCCESS;
}

int bamparsort(libmaus::util::ArgInfo const & arginfo)
{
	libmaus::timing::RealTimeClock rtc;
	rtc.start();
	
	std::string const sortorder = arginfo.getValue<std::string>("SO",getDefaultSortOrder());
	int r = EXIT_FAILURE;
		
	if ( sortorder == "coordinate" )
		r = bamparsortTemplate<libmaus::bambam::BamAlignmentPosComparator>(arginfo,"coordinate");
	else if ( sortorder == "queryname" )
		r = bamparsortTemplate<libmaus::bambam::BamAlignmentNameComparator>(arginfo,"queryname");
	else
	{
		std::cerr << "[E] unknown sort order " << sortorder << std::endl;
		r = EXIT_FAILURE;
	}
	
	if ( r != EXIT_FAILURE )
		std::cerr << "[V] processing finished in time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

	return r;
}


int main(int argc, char * argv[])
{
	try
	{
		libmaus::util::ArgInfo arginfo(argc,argv);
		
		std::cerr << "\n[V] - THIS PROGRAM IS HIGHLY EXPERIMENTAL -\n" << std::endl;

		uint64_t const numthreads = arginfo.getValue<unsigned int>("numthreads", libmaus::parallel::NumCpus::getNumLogicalProcessors());
		
		if ( ! arginfo.hasArg("numthreads") )
		{
			std::ostringstream snumthreads;
			snumthreads << numthreads;
			
			arginfo.argmap["numthreads"] = snumthreads.str();
			arginfo.argmultimap.insert(
				std::pair<std::string,std::string>("numthreads",snumthreads.str())
			);
		}
		
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
				V.push_back ( std::pair<std::string,std::string> ( "mem=<["+libmaus::util::ArgInfo::numToUnitNum(getDefaultMem())+"]>", "memory size target" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfileprefix=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "numthreads=<["+::biobambam::Licensing::formatNumber(arginfo.getValue<unsigned int>("numthreads",1))+"]>", "number of threads" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("tempcomp=<[")+getDefaultTempComp()+"]>", "compression setting for temporary files (zlib:{-1,0,...,9,11},snappy)" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

	
		return bamparsort(arginfo);	
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
