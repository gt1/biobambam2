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
#include <libmaus/types/types.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GrowingFreeList.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/util/SimpleHashMap.hpp>
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

		std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpfile",arginfo.getDefaultTmpFileName());	
		OutputQueue<libmaus::bambam::BamWriter> OQ(wr,BAFL,tmpfilenamebase);

		typedef libmaus::util::unordered_map<PairHashKeyType, libmaus::bambam::BamAlignment *, PairHashKeyTypeHashFunction>::type pair_hash_type;
		typedef libmaus::util::unordered_map<FragmentHashKeyType, libmaus::bambam::BamAlignment *,FragmentHashKeyTypeHashFunction>::type fragment_hash_type;
		typedef libmaus::util::unordered_set<FragmentHashKeyType, FragmentHashKeyTypeHashFunction>::type pair_fragment_hash_type;

		std::priority_queue< PairHashKeyType, std::vector<PairHashKeyType>, PairHashKeyHeapComparator > Qpair;
		pair_hash_type Hpair;

		std::priority_queue< FragmentHashKeyType, std::vector<FragmentHashKeyType>, FragmentHashKeyHeapComparator > Qfragment;
		fragment_hash_type Hfragment;
		
		std::priority_queue< FragmentHashKeyType, std::vector<FragmentHashKeyType>, FragmentHashKeyHeapComparator > Qpairfragments;
		pair_fragment_hash_type Hpairfragments;
		
		uint64_t cnt = 0;

		std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> metrics;
		
		libmaus::timing::RealTimeClock globalrtc;
		globalrtc.start();
		libmaus::timing::RealTimeClock batchrtc;
		batchrtc.start();
					
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
				if ( algn.isRead1() )
					met.readpairsexamined++;
			
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
				
				pair_hash_type::iterator it = Hpair.find(HK);		

				FragmentHashKeyType HK1(algn,header,false);
				if ( Hpairfragments.find(HK1) == Hpairfragments.end() )
				{
					Hpairfragments.insert(HK1);
					Qpairfragments.push(HK1);

					// std::cerr << "Adding " << HK1 << std::endl;
					
					assert ( Hpairfragments.find(HK1) != Hpairfragments.end() );
				}
				FragmentHashKeyType HK2(algn,header,true);
				if ( Hpairfragments.find(HK2) == Hpairfragments.end() )
				{
					Hpairfragments.insert(HK2);
					Qpairfragments.push(HK2);

					// std::cerr << "Adding " << HK2 << std::endl;

					assert ( Hpairfragments.find(HK2) != Hpairfragments.end() );
				}
				
				// type has been seen before
				if ( it != Hpair.end() )
				{
					// other alignment
					libmaus::bambam::BamAlignment * oalgn = it->second;

					// score for other alignment
					int64_t const oscore = oalgn->getScore() + oalgn->getAuxAsNumber<int32_t>("MS");
					// score for this alignment
					int64_t const tscore =   algn.getScore() +   algn.getAuxAsNumber<int32_t>("MS");
					// decide on read name if score is the same
					int64_t const tscoreadd = (tscore == oscore) ? ((strcmp(algn.getName(),oalgn->getName()) < 0)?1:-1) : 0;

					#if 0
					std::cerr << "[V] Comparing for hash key " << HK << "\n" << algn.formatAlignment(dec.getHeader()) << "\n" << 
						oalgn->formatAlignment(dec.getHeader()) << std::endl;
					#endif
					
					#if 0
					if ( tscore == oscore )
					{
						std::cerr << "tscore=" << tscore << " oscore=" << oscore << " tscoreadd=" << tscoreadd << std::endl;
					}
					#endif

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
					
					Hpair[HK] = palgn;
					Qpair.push(HK);
				}
			}
			// single end or pair with one end mapping only
			else
			{
				FragmentHashKeyType HK(algn,header);
								
				fragment_hash_type::iterator it = Hfragment.find(HK);		
		
				if ( it != Hfragment.end() )
				{
					// other alignment
					libmaus::bambam::BamAlignment * oalgn = it->second;
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
					
					Hfragment[HK] = palgn;
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

				pair_hash_type::iterator it = Hpair.find(HK);
				assert ( it != Hpair.end() );
				
				libmaus::bambam::BamAlignment * palgn = it->second;
				OQ.push(palgn);

				Hpair.erase(it);
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

				fragment_hash_type::iterator it = Hfragment.find(HK);
				assert ( it != Hfragment.end() );
				
				libmaus::bambam::BamAlignment * palgn = it->second;
				
				if ( Hpairfragments.find(HK) != Hpairfragments.end() )
				{
					// update metrics
					met.unpairedreadduplicates += 1;
					
					palgn->putFlags(palgn->getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);				
				}
			
				OQ.push(palgn);
				
				Hfragment.erase(it);
			}

			while ( 
				Qpairfragments.size ()
				&&
				(
					(Qpairfragments.top().getRefId() != thisref)
					||
					(Qpairfragments.top().getRefId() == thisref && Qpairfragments.top().getCoord()+3*maxreadlen < thispos)
				)
			)
			{
				FragmentHashKeyType const & HK = Qpairfragments.top();
								
				pair_fragment_hash_type::iterator it = Hpairfragments.find(HK);
				
				assert ( it != Hpairfragments.end() );
				
				Hpairfragments.erase(it);
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
					<< std::endl;
				
				batchrtc.start();
			}
		}

		while ( Qpair.size () )
		{
			PairHashKeyType const HK = Qpair.top(); Qpair.pop();

			pair_hash_type::iterator it = Hpair.find(HK);
			assert ( it != Hpair.end() );
			
			libmaus::bambam::BamAlignment * palgn = it->second;
			OQ.push(palgn);

			Hpair.erase(it);
		}

		while ( Qfragment.size () )
		{
			FragmentHashKeyType const HK = Qfragment.top(); Qfragment.pop();

			fragment_hash_type::iterator it = Hfragment.find(HK);
			assert ( it != Hfragment.end() );
			
			libmaus::bambam::BamAlignment * palgn = it->second;
			
			if ( Hpairfragments.find(HK) != Hpairfragments.end() )
			{
				// update metrics
				metrics[HK.getLibrary()].unpairedreadduplicates += 1;
				
				palgn->putFlags(palgn->getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);				
			}
		
			OQ.push(palgn);
			
			Hfragment.erase(it);
		}

		while ( Qpairfragments.size () )
		{
			FragmentHashKeyType const & HK = Qpairfragments.top();
							
			pair_fragment_hash_type::iterator it = Hpairfragments.find(HK);
			
			assert ( it != Hpairfragments.end() );
			
			Hpairfragments.erase(it);
			Qpairfragments.pop();
		}

		OQ.flush();

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
