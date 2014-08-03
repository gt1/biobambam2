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

// local
#include <config.h>

// std
#include <algorithm>
#include <iostream>
#include <vector>

// libmaus
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/DuplicationMetrics.hpp>
#include <libmaus/bambam/ReadEnds.hpp>
#include <libmaus/lru/SparseLRUFileBunch.hpp>
#include <libmaus/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus/types/types.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/FreeList.hpp>
#include <libmaus/util/GrowingFreeList.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/util/SimpleHashMapInsDel.hpp>
#include <libmaus/util/SimpleHashSetInsDel.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/util/unordered_map.hpp>
#include <libmaus/util/unordered_set.hpp>

struct SignCoding
{
	static int64_t const signshift = std::numeric_limits<int32_t>::min();

	static uint32_t signEncode(int32_t const coord)
	{
		return static_cast<int64_t>(coord)-signshift;
	}
	
	static int32_t signDecode(uint32_t const coord)
	{
		return static_cast<int64_t>(coord)+signshift;
	}
};


struct PairHashKeyType : public SignCoding
{
	typedef PairHashKeyType this_type;

	typedef libmaus::uint::UInt<5> key_type;

	key_type key;

	enum pair_orientation_type { pair_orientaton_FF=0, pair_orientaton_FR=1, pair_orientaton_RF=2, pair_orientaton_RR=3 };
	
	PairHashKeyType() : key() {}
	
	PairHashKeyType(uint64_t const rkey[key_type::words])
	{
		std::copy(&rkey[0],&rkey[key_type::words],&key.A[0]);
	}
	
	PairHashKeyType(
		libmaus::bambam::BamAlignment const & algn,
		libmaus::bambam::BamHeader const & header
	) : key()
	{
		int64_t const thisref = algn.getRefID();
		int64_t const thiscoord = algn.getCoordinate();
		int64_t const otherref = algn.getNextRefID();
		int64_t const othercoord = algn.getAuxAsNumber<int32_t>("MC");
		
		// is this the left mapping end?
		bool const isleft =
			(thisref < otherref) ||
			(thisref == otherref && thiscoord < othercoord ) ||
			(thisref == otherref && thiscoord == othercoord && algn.isRead1());
		
		// as number for hash key
		uint64_t leftflag = isleft ? 0 : 1;
		
		pair_orientation_type orientation;
		
		// orientation of end pair
		if ( isleft )
		{
			if ( ! algn.isReverse() )
			{
				if ( ! algn.isMateReverse() )
					orientation = pair_orientaton_FF;
				else
					orientation = pair_orientaton_FR;
			}
			else
			{
				if ( ! algn.isMateReverse() )
					orientation = pair_orientaton_RF;
				else
					orientation = pair_orientaton_RR;
			}
		}
		else
		{
			if ( ! algn.isMateReverse() )
			{
				if ( ! algn.isReverse() )
					orientation = pair_orientaton_FF;
				else
					orientation = pair_orientaton_FR;
			}
			else
			{
				if ( ! algn.isReverse() )
					orientation = pair_orientaton_RF;
				else
					orientation = pair_orientaton_RR;
			}
		}
		
		// orientation as number		
		uint64_t uorientation = static_cast<uint64_t>(orientation);

		key.A[0] = 
			(static_cast<uint64_t>(signEncode(thisref)) << 32)
			|
			(static_cast<uint64_t>(signEncode(thiscoord)) << 0)
			;
		key.A[1] = 
			(static_cast<uint64_t>(signEncode(otherref)) << 32)
			|
			(static_cast<uint64_t>(signEncode(othercoord)) << 0)
			;
		key.A[2] = 
			(static_cast<uint64_t>(signEncode(algn.getLibraryId(header))) << 32)
			|
			(leftflag) | (uorientation << 1)
			;
		key.A[3] = 0; // tag
		key.A[4] = 0; // tag
	}
	
	bool operator==(this_type const & o) const
	{
		return key == o.key;
	}

	int32_t getRefId() const
	{
		return signDecode(key.A[0] >> 32);
	}

	int32_t getCoord() const
	{
		return signDecode(key.A[0] & 0xFFFFFFFFUL);
	}
	
	int32_t getMateRefId() const
	{
		return signDecode(key.A[1] >> 32);
	}

	int32_t getMateCoord() const
	{
		return signDecode(key.A[1] & 0xFFFFFFFFUL);
	}

	int32_t getLibrary() const
	{
		return signDecode(key.A[2] >> 32);
	}

	pair_orientation_type getOrientation() const
	{
		return static_cast<pair_orientation_type>((key.A[1] >> 1) & 0x3);
	}

	int32_t getLeft() const
	{
		return (((key.A[2] & 0xFFFFFFFFUL) & 1) == 0);
	}

	uint64_t getTag1() const
	{
		return key.A[3];
	}

	uint64_t getTag2() const
	{
		return key.A[4];
	}
};

std::ostream & operator<<(
	std::ostream & out, 
	PairHashKeyType::pair_orientation_type const & ori
)
{
	switch ( ori )
	{
		case PairHashKeyType::pair_orientaton_FF: return out << "FF";
		case PairHashKeyType::pair_orientaton_FR: return out << "FR";
		case PairHashKeyType::pair_orientaton_RF: return out << "RF";
		case PairHashKeyType::pair_orientaton_RR: return out << "RR";
		default: return out << "unknown_orientation";
	}
}

std::ostream & operator<<(std::ostream & ostr, PairHashKeyType const & H)
{
	ostr << "PairHashKeyType(";
	ostr << "refid=" << H.getRefId();
	ostr << ",";
	ostr << "coord=" << H.getCoord();
	ostr << ",";
	ostr << "materefid=" << H.getMateRefId();
	ostr << ",";
	ostr << "matecoord=" << H.getMateCoord();
	ostr << ",";
	ostr << "lib=" << H.getLibrary();
	ostr << ",";
	ostr << "left=" << H.getLeft();
	ostr << ",";
	ostr << "orientation=" << H.getOrientation();
	ostr << ")";
	
	return ostr;
}

struct PairHashKeyHeapComparator
{
	bool operator()(PairHashKeyType const & A, PairHashKeyType const & B)
	{
		for ( unsigned int i = 0; i < A.key.words; ++i )
			if ( A.key.A[i] != B.key.A[i] )
				return A.key.A[i] > B.key.A[i];

		return false;
	}
};

struct PairHashKeyTypeHashFunction
{
	size_t operator()(PairHashKeyType const & H) const
	{
		return libmaus::hashing::EvaHash::hash642(H.key.A,PairHashKeyType::key_type::words);
	}
};

struct FragmentHashKeyType : public SignCoding
{
	typedef FragmentHashKeyType this_type;

	typedef libmaus::uint::UInt<3> key_type;

	key_type key;

	enum fragment_orientation_type { fragment_orientation_F=0, fragment_orientation_R = 1 };
	
	FragmentHashKeyType() : key() {}	

	FragmentHashKeyType(
		libmaus::bambam::BamAlignment const & algn,
		libmaus::bambam::BamHeader const & header,
		bool const usemate = false
	) : key()
	{
		if ( usemate )
		{
			int64_t const thisref = algn.getNextRefID();
			int64_t const thiscoord = algn.getAuxAsNumber<int32_t>("MC");
			fragment_orientation_type fragor = algn.isMateReverse() ? fragment_orientation_R : fragment_orientation_F;
			uint64_t const ufragor = static_cast<uint64_t>(fragor);
			key.A[0] = 
				(static_cast<uint64_t>(signEncode(thisref)) << 32)
				|
				(static_cast<uint64_t>(signEncode(thiscoord)) << 0)
				;
			key.A[1] = (static_cast<uint64_t>(signEncode(algn.getLibraryId(header))) << 32) | ufragor;
			key.A[2] = 0; // tag		
		}
		else
		{
			int64_t const thisref = algn.getRefID();
			int64_t const thiscoord = algn.getCoordinate();
			fragment_orientation_type fragor = algn.isReverse() ? fragment_orientation_R : fragment_orientation_F;
			uint64_t const ufragor = static_cast<uint64_t>(fragor);
			key.A[0] = 
				(static_cast<uint64_t>(signEncode(thisref)) << 32)
				|
				(static_cast<uint64_t>(signEncode(thiscoord)) << 0)
				;
			key.A[1] = (static_cast<uint64_t>(signEncode(algn.getLibraryId(header))) << 32) | ufragor;
			key.A[2] = 0; // tag
		}
	}

	bool operator==(this_type const & o) const
	{
		return key == o.key;
	}

	int32_t getRefId() const
	{
		return signDecode(key.A[0] >> 32);
	}

	int32_t getCoord() const
	{
		return signDecode(key.A[0] & 0xFFFFFFFFUL);
	}

	int32_t getLibrary() const
	{
		return signDecode(key.A[1] >> 32);
	}

	fragment_orientation_type getOrientation() const
	{
		return static_cast<fragment_orientation_type>(key.A[1] & 0x1);
	}
};

std::ostream & operator<<(
	std::ostream & out, 
	FragmentHashKeyType::fragment_orientation_type const & ori
)
{
	switch ( ori )
	{
		case FragmentHashKeyType::fragment_orientation_F: return out << "F";
		case FragmentHashKeyType::fragment_orientation_R: return out << "R";
		default: return out << "unknown_orientation";
	}
}

std::ostream & operator<<(std::ostream & ostr, FragmentHashKeyType const & H)
{
	ostr << "fragment_hash_key(";
	ostr << "refid=" << H.getRefId();
	ostr << ",";
	ostr << "coord=" << H.getCoord();
	ostr << ",";
	ostr << "lib=" << H.getLibrary();
	ostr << ",";
	ostr << "orientation=" << H.getOrientation();
	ostr << " v=" << H.key.A[0] << "," << H.key.A[1] << "," << H.key.A[2];
	ostr << ")";
	
	return ostr;
}

struct FragmentHashKeyHeapComparator
{
	bool operator()(FragmentHashKeyType const & A, FragmentHashKeyType const & B)
	{
		for ( unsigned int i = 0; i < A.key.words; ++i )
			if ( A.key.A[i] != B.key.A[i] )
				return A.key.A[i] > B.key.A[i];

		return false;
	}

};

struct FragmentHashKeyTypeHashFunction
{
	size_t operator()(FragmentHashKeyType const & H) const
	{
		return libmaus::hashing::EvaHash::hash642(H.key.A,FragmentHashKeyType::key_type::words);
	}
};

struct OutputQueueOrder
{
	bool operator()(
		std::pair<uint64_t,libmaus::bambam::BamAlignment *> const & A, 
		std::pair<uint64_t,libmaus::bambam::BamAlignment *> const & B 
	)
	{
		return A.first > B.first;
	}
};


template<typename _writer_type>
struct OutputQueue
{
	typedef _writer_type writer_type;
	typedef OutputQueue<writer_type> this_type;

	writer_type & wr;

	libmaus::util::GrowingFreeList<libmaus::bambam::BamAlignment> & BAFL;	
	
	OutputQueueOrder const order;

	std::priority_queue<
		std::pair<uint64_t,libmaus::bambam::BamAlignment *>,
		std::vector< std::pair<uint64_t,libmaus::bambam::BamAlignment *> >,
		OutputQueueOrder
	> OQ;
		
	int64_t nextout;

	uint64_t const entriespertmpfile;

	libmaus::bambam::BamAuxFilterVector zrtag;

	libmaus::autoarray::AutoArray<libmaus::bambam::BamAlignment *> OL;
	uint64_t olsizefill;
	
	std::string const tmpfileprefix;
	
	std::map<uint64_t,uint64_t> tmpFileFill;
	
	libmaus::lru::SparseLRUFileBunch reorderfiles;

	OutputQueue(
		writer_type & rwr,
		libmaus::util::GrowingFreeList<libmaus::bambam::BamAlignment> & rBAFL,	
		std::string const & rtmpfileprefix,
		uint64_t const rentriespertmpfile = 16*1024
	) 
	: wr(rwr), BAFL(rBAFL), order(), OQ(), nextout(0), 
	  entriespertmpfile(rentriespertmpfile), 
	  OL(entriespertmpfile,false), olsizefill(0), tmpfileprefix(rtmpfileprefix),
	  reorderfiles(tmpfileprefix,16)
	{
		zrtag.set('Z','R');
	}
	
	// flush output list
	void flushOutputList()
	{
		for ( uint64_t i = 0; i < olsizefill; ++i )
		{
			libmaus::bambam::BamAlignment * palgn = OL[i];
			palgn->filterOutAux(zrtag);
			palgn->serialise(wr.getStream());
			BAFL.put(palgn);
		}
		
		olsizefill = 0;
	}
	
	struct BamAlignmentRankComparator
	{
		bool operator()(
			std::pair<uint64_t,libmaus::bambam::BamAlignment *> const & A, 
			std::pair<uint64_t,libmaus::bambam::BamAlignment *> const & B)
		{
			return A.first < B.first;
		}
	};
	
	void flushTempFile(uint64_t const tmpfileindex)
	{
		libmaus::aio::CheckedInputOutputStream & CIOS = reorderfiles[tmpfileindex];
		CIOS.flush();
		CIOS.clear();
		CIOS.seekg(0,std::ios::beg);
		CIOS.clear();
		
		uint64_t const n = tmpFileFill.find(tmpfileindex)->second;
		libmaus::autoarray::AutoArray< std::pair<uint64_t,libmaus::bambam::BamAlignment *> > algns(n);
		
		for ( uint64_t i = 0; i < n; ++i )
		{
			algns[i].second = BAFL.get();
			libmaus::bambam::BamAlignmentDecoder::readAlignmentGz(CIOS,*(algns[i].second),0,false);
			algns[i].first = algns[i].second->getRank();
		}
		
		std::sort(
			algns.begin(),
			algns.begin()+n,
			BamAlignmentRankComparator()
		);

		for ( uint64_t i = 0; i < n; ++i )
		{
			algns[i].second->filterOutAux(zrtag);
			algns[i].second->serialise(wr.getStream());
			BAFL.put(algns[i].second);
		}
		
		nextout += n;
		
		reorderfiles.remove(tmpfileindex);
		tmpFileFill.erase(tmpFileFill.find(tmpfileindex));
		
		// std::cerr << "[V] flushed block " << tmpfileindex << std::endl;
	}
	
	void addTmpFileEntry(libmaus::bambam::BamAlignment * algn)
	{
		addTmpFileEntry(algn, algn->getRank() / entriespertmpfile);
	}
	
	void addTmpFileEntry(
		libmaus::bambam::BamAlignment * algn, uint64_t const tmpfileindex
	)
	{
		// make sure file does not exist before we start adding entries
		if ( tmpFileFill.find(tmpfileindex) == tmpFileFill.end() )
			reorderfiles.remove(tmpfileindex);
	
		libmaus::aio::CheckedInputOutputStream & CIOS = reorderfiles[tmpfileindex];
		algn->serialise(CIOS);
		BAFL.put(algn);
		tmpFileFill[tmpfileindex]++;
		
		// std::cerr << "[V] added to " << tmpfileindex << std::endl;
		
		while ( 
			tmpFileFill.size() 
			&&
			static_cast<uint64_t>(tmpFileFill.begin()->second) == static_cast<uint64_t>(entriespertmpfile)
			&&
			static_cast<uint64_t>(tmpFileFill.begin()->first * entriespertmpfile) == static_cast<uint64_t>(nextout)
		)
		{
			flushTempFile(tmpFileFill.begin()->first);
		}	
	}
	
	void outputListToTmpFiles()
	{
		for ( uint64_t i = 0; i < olsizefill; ++i )
			addTmpFileEntry(OL[i]);
		nextout -= olsizefill;
		olsizefill = 0;	
	}
	
	// flush reorder heap to output list
	void flushInMemQueueInternal()
	{
		while ( 
			OQ.size() && OQ.top().first == static_cast<uint64_t>(nextout)
		)
		{
			libmaus::bambam::BamAlignment * palgn = OQ.top().second;

			OL[olsizefill++] = palgn;
			if ( olsizefill == OL.size() )
				flushOutputList();
			
			nextout += 1;
			OQ.pop();
		}	
	}

	// add alignment
	void push(libmaus::bambam::BamAlignment * algn)
	{
		// tmp file it is assigned to
		uint64_t tmpfileindex;
		
		if ( tmpFileFill.size() && tmpFileFill.find(tmpfileindex=(algn->getRank() / entriespertmpfile)) != tmpFileFill.end() )
		{
			addTmpFileEntry(algn,tmpfileindex);			
		}
		else
		{
			int64_t const rank = algn->getRank();
			assert ( rank >= 0 );
			OQ.push(std::pair<uint64_t,libmaus::bambam::BamAlignment *>(rank,algn));
			flushInMemQueueInternal();
						
			if ( OQ.size() >= 32*1024 )
			{
				outputListToTmpFiles();
				
				while ( OQ.size() )
				{
					addTmpFileEntry(OQ.top().second);
					OQ.pop();
				}
			}
		}
	}
	
	void flush()
	{
		outputListToTmpFiles();

		std::vector<uint64_t> PQ;
		for ( std::map<uint64_t,uint64_t>::const_iterator ita = tmpFileFill.begin(); ita != tmpFileFill.end(); ++ita )
			PQ.push_back(ita->first);

		for ( uint64_t i = 0; i < PQ.size(); ++i )
			flushTempFile(PQ[i]);
	}
};

namespace libmaus
{
	namespace util
	{
		template<>
		struct SimpleHashMapConstants< PairHashKeyType >
		{
			typedef PairHashKeyType::key_type key_type;

			PairHashKeyType const unusedValue;
			PairHashKeyType const deletedValue;
			
			static PairHashKeyType computeUnusedValue()
			{
				key_type U;
				key_type Ulow(std::numeric_limits<uint64_t>::max());
				
				for ( unsigned int i = 0; i < key_type::words; ++i )
				{
					U <<= 64;
					U |= Ulow;
				}
				
				PairHashKeyType PHKT;
				PHKT.key = U;
				
				return PHKT;
			}

			static PairHashKeyType computeDeletedValue()
			{
				// get full mask
				PairHashKeyType U = computeUnusedValue();
				// erase top bit
				U.key.setBit(key_type::words*64-1, 0);
				return U;
			}
						
			PairHashKeyType const & unused() const
			{
				return unusedValue;
			}

			PairHashKeyType const & deleted() const
			{
				return deletedValue;
			}

			bool isFree(PairHashKeyType const & v) const
			{
				return (v.key & deletedValue.key) == deletedValue.key;
			}
			
			bool isInUse(PairHashKeyType const & v) const
			{
				return !isFree(v);
			}
			
			SimpleHashMapConstants() 
			: unusedValue(computeUnusedValue()), deletedValue(computeDeletedValue()) {}
			virtual ~SimpleHashMapConstants() {}
		};
	}
}


namespace libmaus
{
	namespace util
	{	
		template<>
		struct SimpleHashMapHashCompute<PairHashKeyType>
		{
			typedef PairHashKeyType key_type;
			
			inline static uint64_t hash(key_type const v)
			{
				return libmaus::hashing::EvaHash::hash642(&v.key.A[0],v.key.words);
			}
		};
	}
}


namespace libmaus
{
	namespace util
	{
		template<>
		struct SimpleHashMapNumberCast<PairHashKeyType>
		{
			typedef PairHashKeyType key_type;
		
			static uint64_t cast(key_type const & key)
			{
				return key.key.A[0];
			}
		};
	}
}

namespace libmaus
{
	namespace util
	{
		template<>
		struct SimpleHashMapConstants< FragmentHashKeyType >
		{
			typedef FragmentHashKeyType::key_type key_type;

			FragmentHashKeyType const unusedValue;
			FragmentHashKeyType const deletedValue;
			
			static FragmentHashKeyType computeUnusedValue()
			{
				key_type U;
				key_type Ulow(std::numeric_limits<uint64_t>::max());
				
				for ( unsigned int i = 0; i < key_type::words; ++i )
				{
					U <<= 64;
					U |= Ulow;
				}
				
				FragmentHashKeyType PHKT;
				PHKT.key = U;
				
				return PHKT;
			}

			static FragmentHashKeyType computeDeletedValue()
			{
				// get full mask
				FragmentHashKeyType U = computeUnusedValue();
				// erase top bit
				U.key.setBit(key_type::words*64-1, 0);
				return U;
			}
						
			FragmentHashKeyType const & unused() const
			{
				return unusedValue;
			}

			FragmentHashKeyType const & deleted() const
			{
				return deletedValue;
			}

			bool isFree(FragmentHashKeyType const & v) const
			{
				return (v.key & deletedValue.key) == deletedValue.key;
			}
			
			bool isInUse(FragmentHashKeyType const & v) const
			{
				return !isFree(v);
			}
			
			SimpleHashMapConstants() 
			: unusedValue(computeUnusedValue()), deletedValue(computeDeletedValue()) {}
			virtual ~SimpleHashMapConstants() {}
		};
	}
}

namespace libmaus
{
	namespace util
	{	
		template<>
		struct SimpleHashMapHashCompute<FragmentHashKeyType>
		{
			typedef FragmentHashKeyType key_type;
			
			inline static uint64_t hash(key_type const v)
			{
				return libmaus::hashing::EvaHash::hash642(&v.key.A[0],v.key.words);
			}
		};
	}
}

namespace libmaus
{
	namespace util
	{
		template<>
		struct SimpleHashMapNumberCast<FragmentHashKeyType>
		{
			typedef FragmentHashKeyType key_type;
		
			static uint64_t cast(key_type const & key)
			{
				return key.key.A[0];
			}
		};
	}
}

struct OpticalInfoListNode
{
	OpticalInfoListNode * next;

	uint16_t readgroup;	
	uint16_t tile;
	uint32_t x;
	uint32_t y;
	
	OpticalInfoListNode()
	: next(0), readgroup(0), tile(0), x(0), y(0)
	{
	
	}
	
	OpticalInfoListNode(OpticalInfoListNode * rnext, uint16_t rreadgroup, uint16_t rtile, uint32_t rx, uint32_t ry)
	: next(rnext), readgroup(rreadgroup), tile(rtile), x(rx), y(ry)
	{
	
	}
	
	bool operator<(OpticalInfoListNode const & o) const
	{
		if ( readgroup != o.readgroup )
			return readgroup < o.readgroup;
		else if ( tile != o.tile )
			return tile < o.tile;
		else if ( x != o.x )
			return x < o.x;
		else
			return y < o.y;
	}
};

std::ostream & operator<<(std::ostream & out, OpticalInfoListNode const & O)
{
	out << "OpticalInfoListNode(";
	
	out << O.readgroup << "," << O.tile << "," << O.x << "," << O.y;
	
	out << ")";
	return out;
}

struct OpticalInfoListNodeComparator
{
	bool operator()(OpticalInfoListNode const * A, OpticalInfoListNode const * B) const
	{
		return A->operator<(*B);
	}	
};

struct OpticalExternalInfoElement
{
	uint64_t key[PairHashKeyType::key_type::words];
	uint16_t readgroup;
	uint16_t tile;
	uint32_t x;
	uint32_t y;
	
	OpticalExternalInfoElement()
	: readgroup(0), tile(0), x(0), y(0)
	{
		std::fill(&key[0],&key[PairHashKeyType::key_type::words],0);
	}
	
	OpticalExternalInfoElement(PairHashKeyType const & HK, uint16_t const rreadgroup, uint16_t const rtile, uint32_t const rx, uint32_t const ry)
	: readgroup(rreadgroup), tile(rtile), x(rx), y(ry)
	{
		for ( uint64_t i = 0; i < PairHashKeyType::key_type::words; ++i )
			key[i] = HK.key.A[i];
	}

	OpticalExternalInfoElement(PairHashKeyType const & HK, OpticalInfoListNode const & O)
	: readgroup(O.readgroup), tile(O.tile), x(O.x), y(O.y)
	{
		for ( uint64_t i = 0; i < PairHashKeyType::key_type::words; ++i )
			key[i] = HK.key.A[i];
	}
	
	bool operator<(OpticalExternalInfoElement const & o) const
	{
		for ( uint64_t i = 0; i < PairHashKeyType::key_type::words; ++i )
			if ( key[i] != o.key[i] )
				return key[i] < o.key[i];
		
		if ( readgroup != o.readgroup )
			return readgroup < o.readgroup;
		else if ( tile != o.tile )
			return tile < o.tile;
		else if ( x != o.x )
			return x < o.x;
		else
			return y < o.y;
	}
	
	bool operator==(OpticalExternalInfoElement const & o) const
	{
		for ( uint64_t i = 0; i < PairHashKeyType::key_type::words; ++i )
			if ( key[i] != o.key[i] )
				return false;
		
		if ( readgroup != o.readgroup )
			return false;
		else if ( tile != o.tile )
			return false;
		else if ( x != o.x )
			return false;
		else
			return y == o.y;
	}
	
	bool operator!=(OpticalExternalInfoElement const & o) const
	{
		return ! operator==(o);
	}

	bool sameTile(OpticalExternalInfoElement const & o) const
	{
		for ( uint64_t i = 0; i < PairHashKeyType::key_type::words; ++i )
			if ( key[i] != o.key[i] )
				return false;
		
		if ( readgroup != o.readgroup )
			return false;
		else 
			return tile == o.tile;
	}
};

std::ostream & operator<<(std::ostream & out, OpticalExternalInfoElement const & O)
{
	PairHashKeyType HK(O.key);
	out << "OpticalExternalInfoElement(";
	out << HK << ",";
	out << O.readgroup << ",";
	out << O.tile << ",";
	out << O.x << ",";
	out << O.y;
	out << ")";
	
	return out;
}

#include <libmaus/math/iabs.hpp>

struct OpticalInfoList
{
	OpticalInfoListNode * optlist;
	bool expunged;
	
	OpticalInfoList() : optlist(0), expunged(false) {}
	
	static void markYLevel(
		OpticalInfoListNode ** const b,
		OpticalInfoListNode ** l,
		OpticalInfoListNode ** const t,
		bool * const B,
		unsigned int const optminpixeldif
	)
	{
		for ( OpticalInfoListNode ** c = l+1; c != t; ++c )
			if (
				libmaus::math::iabs(
					static_cast<int64_t>((*l)->y)
					-
					static_cast<int64_t>((*c)->y)
				)
				<= optminpixeldif
			)
			{
				B[c-b] = true;			
			}
	}
	
	static void markXLevel(
		OpticalInfoListNode ** const b,
		OpticalInfoListNode ** l,
		OpticalInfoListNode ** const t,
		bool * const B,
		unsigned int const optminpixeldif
	)
	{
		for ( ; l != t; ++l )
		{
			OpticalInfoListNode ** h = l+1;
			
			while ( (h != t) && ((*h)->x-(*l)->x <= optminpixeldif) )
				++h;
				
			markYLevel(b,l,h,B,optminpixeldif);
		}
	}

	static void markTileLevel(
		OpticalInfoListNode ** const b,
		OpticalInfoListNode ** l,
		OpticalInfoListNode ** const t,
		bool * const B,
		unsigned int const optminpixeldif
	)
	{
		while ( l != t )
		{
			OpticalInfoListNode ** h = l+1;
			
			while ( h != t && (*h)->tile == (*l)->tile )
				++h;
			
			markXLevel(b,l,h,B,optminpixeldif);
				
			l = h;
		}
	}

	static void markReadGroupLevel(
		OpticalInfoListNode ** l,
		OpticalInfoListNode ** const t,
		bool * const B,
		unsigned int const optminpixeldif
	)
	{
		OpticalInfoListNode ** const b = l;
	
		while ( l != t )
		{
			OpticalInfoListNode ** h = l+1;
			
			while ( h != t && (*h)->readgroup == (*l)->readgroup )
				++h;

			markTileLevel(b,l,h,B,optminpixeldif);
				
			l = h;
		}	
	}
	
	uint64_t countOpticalDuplicates(
		libmaus::autoarray::AutoArray<OpticalInfoListNode *> & A, 
		libmaus::autoarray::AutoArray<bool> & B, 
		unsigned int const optminpixeldif = 100
	)
	{
		uint64_t opt = 0;

		if ( ! expunged )
		{
			uint64_t n = 0;
			for ( OpticalInfoListNode * cur = optlist; cur; cur = cur->next )
				++n;
			if ( n > A.size() )
				A = libmaus::autoarray::AutoArray<OpticalInfoListNode *>(n);
			if ( n > B.size() )
				B = libmaus::autoarray::AutoArray<bool>(n);

			n = 0;
			for ( OpticalInfoListNode * cur = optlist; cur; cur = cur->next )
				A[n++] = cur;
				
			std::sort ( A.begin(), A.begin()+n, OpticalInfoListNodeComparator() );
			std::fill ( B.begin(), B.begin()+n, false);
			
			markReadGroupLevel(A.begin(),A.begin()+n,B.begin(),optminpixeldif);
			
			for ( uint64_t i = 0; i < n; ++i )
				opt += B[i];
		}
			
		return opt;
	}

	void addOpticalInfo(
		libmaus::bambam::BamAlignment const & algn,
		libmaus::util::FreeList<OpticalInfoListNode> & OILNFL,
		PairHashKeyType const & HK,
		libmaus::sorting::SortingBufferedOutputFile<OpticalExternalInfoElement> & optSBOF,
		libmaus::bambam::BamHeader const & header
	)
	{
		if ( HK.getLeft() )
		{
			uint16_t tile = 0;
			uint32_t x = 0, y = 0;
			if ( libmaus::bambam::ReadEndsBase::parseOptical(reinterpret_cast<uint8_t const *>(algn.getName()),tile,x,y) )
			{
				int64_t const rg = algn.getReadGroupId(header);
				OpticalInfoListNode const newnode(optlist,rg+1,tile,x,y);
				
				if ( expunged || OILNFL.empty() )
				{
					if ( optlist )
					{
						for ( OpticalInfoListNode * cur = optlist; cur; cur = cur->next )
						{
							OpticalExternalInfoElement E(HK,*cur);
							optSBOF.put(E);
						}
						deleteOpticalInfo(OILNFL);
					}
					
					expunged = true;

					OpticalExternalInfoElement E(HK,newnode);
					optSBOF.put(E);
				}
				else
				{
					OpticalInfoListNode * node = OILNFL.get();
					*node = newnode;
					optlist = node;
				}
			} 
		}
	}

	uint64_t deleteOpticalInfo(libmaus::util::FreeList<OpticalInfoListNode> & OILNFL)
	{
		OpticalInfoListNode * node = optlist;
		uint64_t n = 0;
		
		while ( node )
		{
			OpticalInfoListNode * next = node->next;
			OILNFL.put(node);
			node = next;
			++n;
		}
		
		optlist = 0;
				
		return n;			
	}
};

std::ostream & printOpt(std::ostream & out, PairHashKeyType const & HK, uint64_t const opt)
{
	bool const isleft = 
		(HK.getRefId()  < HK.getMateRefId()) ||
		(HK.getRefId() == HK.getMateRefId() && HK.getCoord() < HK.getMateCoord())
	;
	
	if ( isleft )
		std::cerr << "\nopt " 
			<< HK.getLibrary()+1 << " "
			<< HK.getRefId()+1 << " "
			<< HK.getCoord()+1 << " "
			<< HK.getMateRefId()+1 << " "
			<< HK.getMateCoord()+1 << " "
			<< opt
			<< std::endl;
	else
		std::cerr << "\nopt " 
			<< HK.getLibrary()+1 << " "
			<< HK.getMateRefId()+1 << " "
			<< HK.getMateCoord()+1 << " "
			<< HK.getRefId()+1 << " "
			<< HK.getCoord()+1 << " "
			<< opt
			<< std::endl;

	return out;
}

void processOpticalList(
	std::vector<OpticalExternalInfoElement> & optlist,
	std::vector<bool> & optb,
	libmaus::bambam::BamHeader const & header,
	std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> & metrics,
	unsigned int const optminpixeldif
)
{
	// erase boolean vector
	for ( uint64_t i = 0; i < std::min(optb.size(),optlist.size()); ++i )
		optb[i] = false;

	for ( uint64_t i = 0; i+1 < optlist.size(); ++i )
		for ( uint64_t j = i+1; j < optlist.size(); ++j )
			if (
				optlist[j].x - optlist[i].x <= optminpixeldif
				&&
				libmaus::math::iabs(static_cast<int64_t>(optlist[j].y)-static_cast<int64_t>(optlist[i].y)) <= optminpixeldif
			)
			{
				while ( j >= optb.size() )
					optb.push_back(false);

				optb[j] = true;
			}
			
	uint64_t opt = 0;
	for ( uint64_t i = 1; i < std::min(optb.size(),optlist.size()); ++i )
		if ( optb[i] )
			opt += 1;
			
	if ( opt )
	{
		int64_t const thislib = header.getLibraryId(static_cast<int64_t>(optlist[0].readgroup)-1);
		::libmaus::bambam::DuplicationMetrics & met = metrics[thislib];		
		met.opticalduplicates += opt;

		#if 0
		PairHashKeyType const HK(optlist[0].key);
		printOpt(std::cerr,HK,opt);
		#endif
	}
	
	optlist.resize(0);
}


int main(int argc, char *argv[])
{
	try
	{
		libmaus::util::ArgInfo const arginfo(argc,argv);
		int const level = arginfo.getValue<int>("level",-1);
	
		libmaus::bambam::BamDecoder dec(std::cin,true /* put rank */);
		libmaus::bambam::BamHeader const & header = dec.getHeader();
		::libmaus::bambam::BamHeader::unique_ptr_type genuphead(
			libmaus::bambam::BamHeaderUpdate::updateHeader(arginfo,dec.getHeader(),"bamstreamingmarkduplicates",std::string(PACKAGE_VERSION))
		);
		libmaus::bambam::BamWriter wr(std::cout,*genuphead,level);
		libmaus::bambam::BamAlignment & algn = dec.getAlignment();
		
		int64_t const maxreadlen = arginfo.getValue<uint64_t>("maxreadlen",300);
	
		libmaus::util::GrowingFreeList<libmaus::bambam::BamAlignment> BAFL;
		libmaus::util::FreeList<OpticalInfoListNode> OILNFL(32*1024);	

		std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpfile",arginfo.getDefaultTmpFileName());	
		OutputQueue<libmaus::bambam::BamWriter> OQ(wr,BAFL,tmpfilenamebase);

		std::string const optfn = tmpfilenamebase+"_opt";
		libmaus::sorting::SortingBufferedOutputFile<OpticalExternalInfoElement> optSBOF(optfn);
		libmaus::util::TempFileRemovalContainer::addTempFile(optfn);

		std::priority_queue< PairHashKeyType, std::vector<PairHashKeyType>, PairHashKeyHeapComparator > Qpair;
		libmaus::util::SimpleHashMapInsDel<PairHashKeyType,std::pair<libmaus::bambam::BamAlignment *, OpticalInfoList> > SHpair(0);

		std::priority_queue< FragmentHashKeyType, std::vector<FragmentHashKeyType>, FragmentHashKeyHeapComparator > Qfragment;
		libmaus::util::SimpleHashMapInsDel<FragmentHashKeyType,libmaus::bambam::BamAlignment *> SHfragment(0);
		
		std::priority_queue< FragmentHashKeyType, std::vector<FragmentHashKeyType>, FragmentHashKeyHeapComparator > Qpairfragments;
		libmaus::util::SimpleHashSetInsDel<FragmentHashKeyType> SHpairfragments(0);

		libmaus::autoarray::AutoArray<OpticalInfoListNode *> optA;
		libmaus::autoarray::AutoArray<bool> optB;
		unsigned int const optminpixeldif = arginfo.getValue<unsigned int>("optminpixeldif",100);

		uint64_t cnt = 0;

		std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> metrics;
		
		libmaus::timing::RealTimeClock globalrtc;
		globalrtc.start();
		libmaus::timing::RealTimeClock batchrtc;
		batchrtc.start();
		
		double const hashloadfactor = .8;
					
		while ( dec.readAlignment() )
		{
			int64_t const thisref = algn.getRefID();
			int64_t const thispos = algn.getPos();

			uint64_t const thislib = algn.getLibraryId(header);
			::libmaus::bambam::DuplicationMetrics & met = metrics[thislib];
			
			if ( ! algn.isMapped() )
				++met.unmapped;
                        else if ( (!algn.isPaired()) || algn.isMateUnmap() )
                        	++met.unpaired;

			if ( 
				// supplementary alignment
				algn.isSupplementary() 
				|| 
				// secondary alignment
				algn.isSecondary() 
				||
				// single, unmapped
				((!algn.isPaired()) && (!algn.isMapped()))
				||
				// paired, end unmapped
				(algn.isPaired() && (!algn.isMapped()))
			)
			{
				// pass through
				libmaus::bambam::BamAlignment * palgn = BAFL.get();
				palgn->swap(algn);
				OQ.push(palgn);
			}
			// paired end, both mapped
			else if ( 
				algn.isPaired()
				&&
				algn.isMapped()
				&&
				algn.isMateMapped()
			)
			{
				/* 
				 * build the hash key containing
				 *
				 * - ref id
				 * - coordinate (position plus softclipping)
				 * - mate ref id
				 * - mate coordinate (position plus softclipping)
				 * - library id
				 * - flag whether right side read of pair (the one with the higher coordinate)
				 * - orientation
				 * - tag
				 * - mate tag
				 */
				PairHashKeyType HK(algn,header);

				if ( algn.isRead1() )
					met.readpairsexamined++;

				// pair fragment for marking single mappings on same coordinate as duplicates
				FragmentHashKeyType HK1(algn,header);
				if ( !SHpairfragments.contains(HK1) )
				{
					SHpairfragments.insertExtend(HK1,hashloadfactor);
					Qpairfragments.push(HK1);					
					assert ( SHpairfragments.contains(HK1) );
				}
				
				uint64_t keyindex;

				// type has been seen before
				if ( SHpair.containsKey(HK,keyindex) )
				{
					// pair reference
					std::pair<libmaus::bambam::BamAlignment *, OpticalInfoList> & SHpairinfo = SHpair.getValue(keyindex);
					// add optical info
					SHpairinfo.second.addOpticalInfo(algn,OILNFL,HK,optSBOF,header);
					// stored previous alignment
					libmaus::bambam::BamAlignment * oalgn = SHpairinfo.first;
					// score for other alignment
					int64_t const oscore = oalgn->getScore() + oalgn->getAuxAsNumber<int32_t>("MS");
					// score for this alignment
					int64_t const tscore =   algn.getScore() +   algn.getAuxAsNumber<int32_t>("MS");
					// decide on read name if score is the same
					int64_t const tscoreadd = (tscore == oscore) ? ((strcmp(algn.getName(),oalgn->getName()) < 0)?1:-1) : 0;

					// update metrics
					met.readpairduplicates++;

					// this score is greater
					if ( (tscore+tscoreadd) > oscore )
					{
						// mark previous alignment as duplicate
						oalgn->putFlags(oalgn->getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);

						// get alignment from free list
						libmaus::bambam::BamAlignment * palgn = BAFL.get();
						// swap 
						palgn->swap(*oalgn);
						// 
						OQ.push(palgn);
						
						// swap alignment data
						oalgn->swap(algn);
					}
					// other score is greater, mark this alignment as duplicate
					else
					{
						algn.putFlags(algn.getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);
						
						libmaus::bambam::BamAlignment * palgn = BAFL.get();
						palgn->swap(algn);
						OQ.push(palgn);
					}
				}
				// type is new
				else
				{
					libmaus::bambam::BamAlignment * palgn = BAFL.get();
					palgn->swap(algn);
					
					Qpair.push(HK);
					
					keyindex = SHpair.insertExtend(
						HK,
						std::pair<libmaus::bambam::BamAlignment *, OpticalInfoList>(palgn,OpticalInfoList()),
						hashloadfactor
					);
					// pair reference
					std::pair<libmaus::bambam::BamAlignment *, OpticalInfoList> & SHpairinfo = SHpair.getValue(keyindex);
					// add optical info
					SHpairinfo.second.addOpticalInfo(*(SHpairinfo.first),OILNFL,HK,optSBOF,header);
				}
			}
			// single end or pair with one end mapping only
			else
			{
				FragmentHashKeyType HK(algn,header);

				libmaus::bambam::BamAlignment * oalgn;
		
				if ( SHfragment.contains(HK,oalgn) )
				{
					// score for other alignment
					int64_t const oscore = oalgn->getScore();
					// score for this alignment
					int64_t const tscore = algn.getScore();

					// update metrics
					met.unpairedreadduplicates += 1;
					
					// this score is greater
					if ( tscore > oscore )
					{
						// mark previous alignment as duplicate
						oalgn->putFlags(oalgn->getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);

						// put oalgn value in output queue
						libmaus::bambam::BamAlignment * palgn = BAFL.get();
						palgn->swap(*oalgn);
						OQ.push(palgn);

						// swap alignment data
						oalgn->swap(algn);
					}
					// other score is greater, mark this alignment as duplicate
					else
					{
						algn.putFlags(algn.getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);

						libmaus::bambam::BamAlignment * palgn = BAFL.get();
						palgn->swap(algn);
						OQ.push(palgn);
					}
				}
				else
				{
					libmaus::bambam::BamAlignment * palgn = BAFL.get();
					palgn->swap(algn);
					
					SHfragment.insertExtend(HK,palgn,hashloadfactor);
					Qfragment.push(HK);
				}
			}
			
			while ( 
				Qpair.size ()
				&&
				(
					(Qpair.top().getRefId() != thisref)
					||
					(Qpair.top().getRefId() == thisref && Qpair.top().getCoord()+maxreadlen < thispos)
				)
			)
			{
				PairHashKeyType const HK = Qpair.top(); Qpair.pop();

				uint64_t const SHpairindex = SHpair.getIndexUnchecked(HK);
				std::pair<libmaus::bambam::BamAlignment *, OpticalInfoList> SHpairinfo = SHpair.getValue(SHpairindex);
				libmaus::bambam::BamAlignment * palgn = SHpairinfo.first;

				uint64_t const opt = SHpairinfo.second.countOpticalDuplicates(optA,optB,optminpixeldif);
				
				if ( opt )
				{
					uint64_t const thislib = palgn->getLibraryId(header);
					::libmaus::bambam::DuplicationMetrics & met = metrics[thislib];	
					met.opticalduplicates += opt;
					
					#if 0		
					printOpt(std::cerr,HK,opt);
					#endif
				}

				OQ.push(palgn);

				SHpairinfo.second.deleteOpticalInfo(OILNFL);
				SHpair.eraseIndex(SHpairindex);
			}

			while ( 
				Qfragment.size ()
				&&
				(
					(Qfragment.top().getRefId() != thisref)
					||
					(Qfragment.top().getRefId() == thisref && Qfragment.top().getCoord()+maxreadlen < thispos)
				)
			)
			{
				FragmentHashKeyType const HK = Qfragment.top(); Qfragment.pop();

				uint64_t const SHfragmentindex = SHfragment.getIndexUnchecked(HK);
				libmaus::bambam::BamAlignment * palgn = SHfragment.getValue(SHfragmentindex);
				
				if ( SHpairfragments.contains(HK) )
				{
					// update metrics
					met.unpairedreadduplicates += 1;
					
					palgn->putFlags(palgn->getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);				
				}
			
				OQ.push(palgn);
				
				SHfragment.eraseIndex(SHfragmentindex);
			}

			while ( 
				Qpairfragments.size ()
				&&
				(
					(Qpairfragments.top().getRefId() != thisref)
					||
					(Qpairfragments.top().getRefId() == thisref && Qpairfragments.top().getCoord()+maxreadlen < thispos)
				)
			)
			{
				FragmentHashKeyType const & HK = Qpairfragments.top();
				
				SHpairfragments.eraseIndex(SHpairfragments.getIndexUnchecked(HK));
				
				Qpairfragments.pop();
			}
			
			if ( (++cnt % (1024*1024)) == 0 )
			{
				std::cerr 
					<< "[V] " 
					<< cnt << " " 
					<< OQ.nextout << " " 
					<< Qpair.size() << " " 
					<< Qfragment.size() << " " 
					<< libmaus::util::MemUsage() << " "
					<< globalrtc.formatTime(globalrtc.getElapsedSeconds()) << " "
					<< batchrtc.formatTime(batchrtc.getElapsedSeconds())
					<< " " << SHpair.getTableSize()
					<< std::endl;
				
				batchrtc.start();
			}
		}

		while ( Qpair.size () )
		{
			PairHashKeyType const HK = Qpair.top(); Qpair.pop();
			
			uint64_t const SHpairindex = SHpair.getIndexUnchecked(HK);
			std::pair<libmaus::bambam::BamAlignment *, OpticalInfoList> SHpairinfo = SHpair.getValue(SHpairindex);
			libmaus::bambam::BamAlignment * palgn = SHpairinfo.first;
				
			OQ.push(palgn);

			uint64_t const opt = SHpairinfo.second.countOpticalDuplicates(optA,optB,optminpixeldif);
				
			if ( opt )
			{
				uint64_t const thislib = palgn->getLibraryId(header);
				::libmaus::bambam::DuplicationMetrics & met = metrics[thislib];		
				met.opticalduplicates += opt;

				#if 0
				printOpt(std::cerr,HK,opt);
				#endif
			}

			SHpairinfo.second.deleteOpticalInfo(OILNFL);
			SHpair.eraseIndex(SHpairindex);
		}

		while ( Qfragment.size () )
		{
			FragmentHashKeyType const HK = Qfragment.top(); Qfragment.pop();

			uint64_t const SHfragmentindex = SHfragment.getIndexUnchecked(HK);
			libmaus::bambam::BamAlignment * palgn = SHfragment.getValue(SHfragmentindex);
			
			if ( SHpairfragments.contains(HK) )
			{
				// update metrics
				metrics[HK.getLibrary()].unpairedreadduplicates += 1;
				
				palgn->putFlags(palgn->getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);				
			}
		
			OQ.push(palgn);

			SHfragment.eraseIndex(SHfragmentindex);
		}

		while ( Qpairfragments.size () )
		{
			FragmentHashKeyType const & HK = Qpairfragments.top();
			
			SHpairfragments.eraseIndex(SHpairfragments.getIndexUnchecked(HK));

			Qpairfragments.pop();
		}

		OQ.flush();
		
		// process expunged optical data
		libmaus::sorting::SortingBufferedOutputFile<OpticalExternalInfoElement>::merger_ptr_type Poptmerger =
			optSBOF.getMerger();
		libmaus::sorting::SortingBufferedOutputFile<OpticalExternalInfoElement>::merger_type & optmerger =
			*Poptmerger;

		OpticalExternalInfoElement optin;
		std::vector<OpticalExternalInfoElement> optlist;
		std::vector<bool> optb;
		
		while ( optmerger.getNext(optin) )
		{
			if ( optlist.size() && !optlist.back().sameTile(optin) )
				processOpticalList(optlist,optb,header,metrics,optminpixeldif);
			
			optlist.push_back(optin);
		}
		
		if ( optlist.size() )
			processOpticalList(optlist,optb,header,metrics,optminpixeldif);
		
		// we have counted duplicated pairs for both ends, so divide by two
		for ( std::map<uint64_t,::libmaus::bambam::DuplicationMetrics>::iterator ita = metrics.begin(); ita != metrics.end();
			++ita )
			ita->second.readpairduplicates /= 2;

		/**
		 * write metrics
		 **/
		::libmaus::aio::CheckedOutputStream::unique_ptr_type pM;
		std::ostream * pmetricstr = 0;
		
		if ( arginfo.hasArg("M") && (arginfo.getValue<std::string>("M","") != "") )
		{
			::libmaus::aio::CheckedOutputStream::unique_ptr_type tpM(
					new ::libmaus::aio::CheckedOutputStream(arginfo.getValue<std::string>("M",std::string("M")))
				);
			pM = UNIQUE_PTR_MOVE(tpM);
			pmetricstr = pM.get();
		}
		else
		{
			pmetricstr = & std::cerr;
		}

		std::ostream & metricsstr = *pmetricstr;

		::libmaus::bambam::DuplicationMetrics::printFormatHeader(arginfo.commandline,metricsstr);
		for ( std::map<uint64_t,::libmaus::bambam::DuplicationMetrics>::const_iterator ita = metrics.begin(); ita != metrics.end();
			++ita )
			ita->second.format(metricsstr, header.getLibraryName(ita->first));
		
		if ( metrics.size() == 1 )
		{
			metricsstr << std::endl;
			metricsstr << "## HISTOGRAM\nBIN\tVALUE" << std::endl;
			metrics.begin()->second.printHistogram(metricsstr);
		}
		
		metricsstr.flush();
		pM.reset();
		/*
		 * end of metrics file writing
		 */
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}
}
