/**
    bambam
    Copyright (C) 2009-2013 German Tischler
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
#include <config.h>

#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/aio/CheckedInputStream.hpp>
#include <libmaus/aio/CheckedOutputStream.hpp>
#include <libmaus/aio/SynchronousGenericInput.hpp>
#include <libmaus/bambam/BamAlignmentFreeList.hpp>
#include <libmaus/bambam/BamAlignmentPairFreeList.hpp>
#include <libmaus/bambam/ReadEndsFreeList.hpp>
#include <libmaus/bambam/BamAlignmentSnappyInput.hpp>
#include <libmaus/bambam/BamAlignmentInputCallbackBam.hpp>
#include <libmaus/bambam/BamAlignmentInputCallbackSnappy.hpp>
#include <libmaus/bambam/BamAlignmentInputPositionCallbackNull.hpp>
#include <libmaus/bambam/BamAlignmentInputCallbackSnappy.hpp>
#include <libmaus/bambam/BamAlignmentInputCallbackBam.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamParallelRewrite.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamAlignmentInputPositionUpdateCallback.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/DuplicationMetrics.hpp>
#include <libmaus/bambam/DupMarkBase.hpp>
#include <libmaus/bambam/DupSetCallbackStream.hpp>
#include <libmaus/bambam/DupSetCallbackSet.hpp>
#include <libmaus/bambam/DupSetCallbackVector.hpp>
#include <libmaus/bambam/OpticalComparator.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/bambam/ReadEndsContainer.hpp>
#include <libmaus/bambam/SortedFragDecoder.hpp>
#include <libmaus/fastx/FastATwoBitTable.hpp>
#include <libmaus/bitio/BitVector.hpp>
#include <libmaus/lz/BgzfInflateDeflateParallel.hpp>
#include <libmaus/lz/BgzfRecode.hpp>
#include <libmaus/lz/BgzfRecodeParallel.hpp>
#include <libmaus/math/iabs.hpp>
#include <libmaus/math/numbits.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/trie/SimpleTrie.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/ContainerGetObject.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1048576;  }
static bool getDefaultRewriteBam() { return 0; }
static int getDefaultRewriteBamLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultColHashBits() { return 20; }
static uint64_t getDefaultColListSize() { return 32*1024*1024; }
static uint64_t getDefaultFragBufSize() { return 48*1024*1024; }
static uint64_t getDefaultMarkThreads() { return 1; }
static uint64_t getDefaultMaxReadLength() { return 500; }
static bool getDefaultRmDup() { return 0; }
static std::string getProgId() { return "bammarkduplicates2"; }
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }
static uint64_t getDefaultInputBufferSize() { return 64*1024; }

#include <libmaus/aio/PosixFdInput.hpp>

/**
 * class storing elements for a given coordinate
 **/
template<typename _free_list_type>
struct PairActiveCountTemplate
{
	typedef _free_list_type free_list_type;
	typedef typename free_list_type::element_type list_node_type;
	typedef PairActiveCountTemplate<free_list_type> this_type;

	int32_t refid;
	int32_t coordinate;
	uint64_t incnt;
	uint64_t outcnt;	
	list_node_type * root;
	bool expunge;
	
	PairActiveCountTemplate() : refid(-1), coordinate(-1), incnt(0), outcnt(0), root(0), expunge(false) {}
	PairActiveCountTemplate(
		int32_t const rrefid,
		int32_t const rcoordinate,
		uint64_t const rincnt,
		uint64_t const routcnt
	) : refid(rrefid), coordinate(rcoordinate), incnt(rincnt), outcnt(routcnt), root(0), expunge(false)
	{	
	}
	~PairActiveCountTemplate()
	{
	}
	
	void setExpunge(bool const rexpunge)
	{
		expunge = rexpunge;
	}
	
	bool getExpunge() const
	{
		return expunge;
	}
	
	bool operator==(std::pair<int32_t,int32_t> const & o)
	{
		return o.first == refid && o.second == coordinate;
	}
	
	bool operator<(this_type const & o) const
	{
		if ( refid != o.refid )
			return refid < o.refid;
		else if ( coordinate != o.coordinate )
			return coordinate < o.coordinate;
		else
			return false;
	}
	
	void incIn()
	{
		++incnt;
	}

	void incOut()
	{
		++outcnt;
	}
	
	void freeAlignments(free_list_type & list)
	{
		list_node_type * p = root;
		
		while ( p )
		{
			list_node_type * q = p->next;
			p->next = 0;
			list.put(p);
			p = q;
		}
		
		root = 0;
	}

	void addAlignmentPair(list_node_type * ptr)
	{
		ptr->next = root;
		root = ptr;
	}
};

template<typename free_list_type>
std::ostream & operator<<(std::ostream & out, PairActiveCountTemplate<free_list_type> const & A)
{
	out << "PairActiveCount(" << A.refid << "," << A.coordinate << "," << A.incnt << "," << A.outcnt << ")";
	return out;
}

typedef PairActiveCountTemplate<libmaus::bambam::BamAlignmentPairFreeList> BamPairActiveCount;
typedef PairActiveCountTemplate<libmaus::bambam::ReadEndsFreeList> ReadEndsActiveCount;

struct BamAlignmentInputPositionCallbackDupMark : public libmaus::bambam::BamAlignmentInputPositionUpdateCallback
{
	static unsigned int const defaultfreelistsize = 16*1024;

	// static int const default_maxreadlength = 250;
	static int getDefaultMaxReadLength() { return 250; }
	static uint64_t getDefaultFreeListSize() { return defaultfreelistsize; }
	
	#define POS_READ_ENDS

	#if defined(POS_READ_ENDS)
	typedef ReadEndsActiveCount active_count_type;
	#else
	typedef BamPairActiveCount active_count_type;
	#endif


	libmaus::bambam::BamHeader const & bamheader;

	libmaus::bambam::ReadEndsBasePointerComparator const REcomp;

	std::pair<int32_t,int32_t> position;

	std::pair<int32_t,int32_t> expungepositionpairs;
	std::deque<active_count_type> activepairs;
	int64_t totalactivepairs;
	active_count_type::free_list_type APFLpairs;
	uint64_t excntpairs;
	uint64_t fincntpairs;
	uint64_t strcntpairs;
	#if defined(POS_READ_ENDS)
	std::vector<libmaus::bambam::ReadEnds *> REpairs;
	#else
	std::vector<libmaus::bambam::ReadEnds> REpairs;
	#endif

	std::pair<int32_t,int32_t> expungepositionfrags;
	std::deque<active_count_type> activefrags;
	int64_t totalactivefrags;
	active_count_type::free_list_type APFLfrags;
	uint64_t excntfrags;
	uint64_t fincntfrags;
	uint64_t strcntfrags;
	std::vector<libmaus::bambam::ReadEnds *> REfrags;

	::libmaus::bambam::DupSetCallback * DSC;
	
	int maxreadlength;

	BamAlignmentInputPositionCallbackDupMark(libmaus::bambam::BamHeader const & rbamheader, uint64_t const freelistsize = defaultfreelistsize)
	: 
		bamheader(rbamheader), REcomp(), position(-1,-1), 
		expungepositionpairs(-1,-1), activepairs(), totalactivepairs(0), APFLpairs(freelistsize), excntpairs(0), fincntpairs(0), strcntpairs(0), REpairs(),
		expungepositionfrags(-1,-1), activefrags(), totalactivefrags(0), APFLfrags(2*freelistsize), excntfrags(0), fincntfrags(0), strcntfrags(0), REfrags(),
		DSC(0),
		maxreadlength(getDefaultMaxReadLength())
	{
	}
	
	void setDupSetCallback(::libmaus::bambam::DupSetCallback * rDSC)
	{
		DSC = rDSC;
	}
	
	void setMaxReadLength(int const rmaxreadlength)
	{
		maxreadlength = rmaxreadlength;
	}

	virtual ~BamAlignmentInputPositionCallbackDupMark() {}

	/**
	 * check whether this alignment is part of an innie pair
	 **/
	static bool isSimplePair(::libmaus::bambam::BamAlignment const & A)
	{
		// both ends need to be mapped
		if ( ! (A.isMapped() && A.isMateMapped()) )
			return false;

		// mapped to the same reference sequence
		if ( A.getRefID() != A.getNextRefID() )
			return false;
		
		int const rev1 = A.isReverse() ? 1 : 0;
		int const rev2 = A.isMateReverse() ? 1 : 0;
		
		// one forward, one reverse
		if ( rev1 + rev2 != 1 )
			return false;

		// reverse read needs to map behind the forward read
		if ( rev2 )
			return A.getPos() < A.getNextPos();
		else
			return A.getNextPos() < A.getPos();
	}
	
	/**
	 * update input position
	 **/	
	void updatePosition(::libmaus::bambam::BamAlignment const & A)
	{
		int32_t const refid = A.getRefID();
		int32_t const pos = A.getPos();
		
		// std::cerr << "refid=" << refid << ",pos=" << pos << std::endl;
	
		position.first = refid;
		position.second = pos;
		
		if ( 
			static_cast<uint32_t>(refid) < static_cast<uint32_t>(position.first) 
			||
			(
				static_cast<uint32_t>(refid) == static_cast<uint32_t>(position.first) 
				&&
				static_cast<uint32_t>(pos) < static_cast<uint32_t>(position.second) 
			)
		)
			{
				libmaus::exception::LibMausException se;
				se.getStream() << "Input file is not sorted by coordinate." << std::endl;
				se.finish();
				throw se;
			}

		int32_t const coord = A.getCoordinate();
		std::pair<int32_t,int32_t> pcoord(refid,coord);
		active_count_type acomp(refid,coord,0,0);
		
		/*
		 * update fragment list
		 */
		if ( ! activefrags.size() || activefrags.back() < acomp )
		{
			// insert at back of list
			activefrags.push_back(active_count_type(refid,coord,1,0));		
		}
		else if ( activefrags.back().refid == refid && activefrags.back().coordinate == coord )
		{
			// increment at back of list
			activefrags.back().incIn();
		}
		else
		{
			// find position
			std::deque<active_count_type>::iterator it = 
				std::lower_bound(activefrags.begin(),activefrags.end(),acomp);

			if ( it != activefrags.end() && *it == pcoord )
			{
				// increment existing
				it->incIn();
			}
			else
			{			
				// insert inside list	
				uint64_t const tomove = activefrags.end()-it;
					
				activefrags.push_back(active_count_type());
				
				for ( uint64_t j = 0; j < tomove; ++j )
					activefrags [ activefrags.size()-j-1 ] = activefrags[activefrags.size()-j-2];
					
				activefrags [ activefrags.size() - tomove - 1 ] = active_count_type(refid,coord,1,0);
			}
			
		}

		totalactivefrags++;

		/*
		 * update pair list
		 */
		if ( isSimplePair(A) && A.isReverse() )
		{	
			// we have not seen the coordinate before
			if ( ! activepairs.size() || activepairs.back() < acomp )
			{
				activepairs.push_back(active_count_type(refid,coord,1,0));
			}
			// increment at end
			else if ( activepairs.back().refid == refid && activepairs.back().coordinate == coord )
			{
				activepairs.back().incIn();
			}
			else
			{
				std::deque<active_count_type>::iterator it = 
					std::lower_bound(activepairs.begin(),activepairs.end(),acomp);

				if ( it != activepairs.end() && *it == pcoord )
				{
					it->incIn();
				}
				else
				{				
					uint64_t const tomove = activepairs.end()-it;
					
					activepairs.push_back(active_count_type());
					
					for ( uint64_t j = 0; j < tomove; ++j )
						activepairs [ activepairs.size()-j-1 ] = activepairs[activepairs.size()-j-2];
						
					activepairs [ activepairs.size() - tomove - 1 ] = active_count_type(refid,coord,1,0);
				}
			}

			totalactivepairs++;
		}
	}

	void finishActiveFrontFrags()
	{	
		assert ( activefrags.size() );
	
		active_count_type & AC = activefrags.front();
		
		uint64_t lfincntfrags = 0;
		for ( active_count_type::free_list_type::element_type * ptr = AC.root; ptr; ptr = ptr->next )
		{
			if ( lfincntfrags < REfrags.size() )
				REfrags[lfincntfrags] = &(ptr->A);
			else
				REfrags.push_back(&(ptr->A));

			lfincntfrags += 1;
		}

		std::sort(REfrags.begin(),REfrags.begin()+lfincntfrags,REcomp);
		
		fincntfrags += lfincntfrags;
		
		uint64_t l = 0;
		while ( l != lfincntfrags )
		{
			uint64_t h = l+1;
			while ( h != lfincntfrags && libmaus::bambam::DupMarkBase::isDupFrag(*REfrags[l],*REfrags[h]) )
				++h;
				
			if ( h-l > 1 )
				libmaus::bambam::DupMarkBase::markDuplicateFragsPointers(REfrags.begin()+l,REfrags.begin()+h,*DSC);
			
			l = h;
		}
		
		AC.freeAlignments(APFLfrags);
		totalactivefrags -= AC.incnt;
		assert ( totalactivefrags >= 0 );

		// set new frag expunge position
		if ( 
			AC.refid > expungepositionfrags.first
			||
			(
				AC.refid == expungepositionfrags.first
				&&
				AC.coordinate > expungepositionfrags.second
			)
		)
		{
			expungepositionfrags.first = AC.refid;
			expungepositionfrags.second = AC.coordinate;
		}

		activefrags.pop_front();	
	}

	void finishActiveFrontPairs()
	{	
		assert ( activepairs.size() );
	
		active_count_type & AC = activepairs.front();

		uint64_t lfincntpairs = 0;
		for ( active_count_type::free_list_type::element_type * ptr = AC.root; ptr; ptr = ptr->next )
		{
			if ( lfincntpairs < REpairs.size() )
				#if defined(POS_READ_ENDS)
				REpairs[lfincntpairs] = &(ptr->A);
				#else
				libmaus::bambam::ReadEndsBase::fillFragPair(
					ptr->A[0],
					ptr->A[1],
					bamheader,
					REpairs[lfincntpairs]);
				#endif
			else
			{
				#if defined(POS_READ_ENDS)
				REpairs.push_back(&(ptr->A));
				#else
				REpairs.push_back(libmaus::bambam::ReadEnds(ptr->A[0],ptr->A[1],bamheader));
				#endif
			}

			lfincntpairs += 1;
		}

		#if defined(POS_READ_ENDS)
		std::sort(REpairs.begin(),REpairs.begin()+lfincntpairs,REcomp);
		#else
		std::sort(REpairs.begin(),REpairs.begin()+lfincntpairs);		
		#endif
		
		fincntpairs += lfincntpairs;
		
		uint64_t l = 0;
		while ( l != lfincntpairs )
		{
			uint64_t h = l+1;
			while ( h != lfincntpairs && libmaus::bambam::DupMarkBase::isDupPair(*REpairs[l],*REpairs[h]) )
				++h;
				
			if ( h-l > 1 )
			{
				#if defined(POS_READ_ENDS)
				libmaus::bambam::DupMarkBase::markDuplicatePairsPointers(REpairs.begin()+l,REpairs.begin()+h,*DSC);
				#else
				libmaus::bambam::DupMarkBase::markDuplicatePairs(REpairs.begin()+l,REpairs.begin()+h,*DSC);
				#endif			
			}
			
			l = h;
		}
		
		AC.freeAlignments(APFLpairs);
		totalactivepairs -= AC.incnt;
		assert ( totalactivepairs >= 0 );

		// set new pair expunge position
		if ( 
			AC.refid > expungepositionpairs.first 
			||
			(
				AC.refid == expungepositionpairs.first
				&&
				AC.coordinate > expungepositionpairs.second
			)
		)
		{
			expungepositionpairs.first = AC.refid;
			expungepositionpairs.second = AC.coordinate;
		}
		
		activepairs.pop_front();		
	}

	void expungeActiveFrontPairs(
		::libmaus::bambam::ReadEndsContainer * pairREC, ::libmaus::bambam::BamHeader const & /* header */
	)
	{
		assert ( activepairs.size() );
	
		active_count_type & AC = activepairs.front();
		
		uint64_t lexcntpairs = 0;
		for ( active_count_type::free_list_type::element_type * ptr = AC.root; ptr; ptr = ptr->next )
		{
			#if defined(POS_READ_ENDS)
			pairREC->put(ptr->A);
			#else
			pairREC->putPair(ptr->A[0],ptr->A[1],header);
			#endif
			lexcntpairs += 1;
		}
		
		excntpairs += lexcntpairs;
		
		if ( 
			AC.refid > expungepositionpairs.first 
			||
			(
				AC.refid == expungepositionpairs.first
				&&
				AC.coordinate > expungepositionpairs.second
			)
		)
		{
			expungepositionpairs.first = AC.refid;
			expungepositionpairs.second = AC.coordinate;
		}
		
		AC.freeAlignments(APFLpairs);

		activepairs.pop_front();	
	}

	void expungeActiveFrontFrags(::libmaus::bambam::ReadEndsContainer * fragREC, ::libmaus::bambam::BamHeader const & /* header */)
	{
		assert ( activefrags.size() );
	
		active_count_type & AC = activefrags.front();

		uint64_t lexcntfrags = 0;
		for ( active_count_type::free_list_type::element_type * ptr = AC.root; ptr; ptr = ptr->next )
		{
			fragREC->put(ptr->A);
			lexcntfrags += 1;
		}
		
		excntfrags += lexcntfrags;

		if ( 
			AC.refid > expungepositionfrags.first
			||
			(
				AC.refid == expungepositionfrags.first
				&&
				AC.coordinate > expungepositionfrags.second
			)
		)
		{
			expungepositionfrags.first = AC.refid;
			expungepositionfrags.second = AC.coordinate;
		}
		AC.freeAlignments(APFLfrags);

		activefrags.pop_front();	
	}

	bool isActivePair(::libmaus::bambam::BamAlignment const & B)
	{
		bool const isactivepairs = 
			B.getRefID() > expungepositionpairs.first
			||
			(
				B.getRefID() == expungepositionpairs.first
				&&
				B.getCoordinate() > expungepositionpairs.second
			)
		;
	
		return isactivepairs;
	}

	bool isActiveFrag(::libmaus::bambam::BamAlignment const & B)
	{
		bool const isactivefrag = 
			B.getRefID() > expungepositionfrags.first
			||
			(
				B.getRefID() == expungepositionfrags.first
				&&
				B.getCoordinate() > expungepositionfrags.second
			)
		;
	
		return isactivefrag;
	}
	
	void setExpungePairs(libmaus::bambam::BamAlignment const & B)
	{
		if ( isActivePair(B) )
		{
			active_count_type const bkey(B.getRefID(),B.getCoordinate(),0,0);

			// key is not there or we can insert at the end, insert it
			if ( ! activepairs.size() || activepairs.back() < bkey )
			{
				activepairs.push_back(bkey);
				activepairs.back().setExpunge(true);
			}
			else 
			{
				std::deque<active_count_type>::iterator const it = std::lower_bound(activepairs.begin(),activepairs.end(),bkey);
			
				// key is already there	
				if ( 
					it != activepairs.end() &&
					it->refid == bkey.refid &&
					it->coordinate == bkey.coordinate )
				{
					it->setExpunge(true);
				}
				// key is not there and needs to be inserted ahead of the end
				else
				{
					uint64_t const tomove = activepairs.end()-it;
						
					activepairs.push_back(active_count_type());
								
					for ( uint64_t j = 0; j < tomove; ++j )
					activepairs [ activepairs.size()-j-1 ] = activepairs[activepairs.size()-j-2];
								
					activepairs [ activepairs.size() - tomove - 1 ] = bkey;
					activepairs [ activepairs.size() - tomove - 1 ].setExpunge(true);
				}
			}
		}
	}

	void setExpungeFrag(libmaus::bambam::BamAlignment const & B)
	{
		if ( isActiveFrag(B) )
		{
			active_count_type const bkey(B.getRefID(),B.getCoordinate(),0,0);

			// key is not there or we can insert at the end, insert it
			if ( ! activefrags.size() || activefrags.back() < bkey )
			{
				activefrags.push_back(bkey);
				activefrags.back().setExpunge(true);
			}
			else 
			{
				std::deque<active_count_type>::iterator const it = std::lower_bound(activefrags.begin(),activefrags.end(),bkey);
			
				// key is already there	
				if ( 
					it != activefrags.end() &&
					it->refid == bkey.refid &&
					it->coordinate == bkey.coordinate )
				{
					it->setExpunge(true);
				}
				// key is not there and needs to be inserted ahead of the end
				else
				{
					uint64_t const tomove = activefrags.end()-it;
						
					activefrags.push_back(active_count_type());
								
					for ( uint64_t j = 0; j < tomove; ++j )
					activefrags [ activefrags.size()-j-1 ] = activefrags[activefrags.size()-j-2];
								
					activefrags [ activefrags.size() - tomove - 1 ] = bkey;
					activefrags [ activefrags.size() - tomove - 1 ].setExpunge(true);
				}
			}
		}
	}
	
	void expungeUntilPairs(
		libmaus::bambam::BamAlignment const & A,
		::libmaus::bambam::ReadEndsContainer * pairREC, 
		::libmaus::bambam::BamHeader const & header
	)
	{
		if ( isActivePair(A) )
		{
			while ( 
				activepairs.size()
				&&
				(
					(
						activepairs.front().refid < A.getRefID()
					)
					||
					(
						activepairs.front().refid == A.getRefID()
						&&
						activepairs.front().coordinate <= A.getCoordinate()
					)
				)
			)
			{
				expungeActiveFrontPairs(pairREC,header);
			}

			expungepositionpairs.first = A.getRefID();
			expungepositionpairs.second = A.getCoordinate();
		}
	}

	void expungeUntilFrag(
		libmaus::bambam::BamAlignment const & A,
		::libmaus::bambam::ReadEndsContainer * fragREC, 
		::libmaus::bambam::BamHeader const & header
	)
	{
		if ( isActiveFrag(A) )
		{
			while ( 
				activefrags.size()
				&&
				(
					(
						activefrags.front().refid < A.getRefID()
					)
					||
					(
						activefrags.front().refid == A.getRefID()
						&&
						activefrags.front().coordinate <= A.getCoordinate()
					)
				)
			)
			{
				expungeActiveFrontFrags(fragREC,header);
			}

			expungepositionfrags.first = A.getRefID();
			expungepositionfrags.second = A.getCoordinate();
		}
	}
	
	void flushPairs(::libmaus::bambam::ReadEndsContainer * pairREC, ::libmaus::bambam::BamHeader const & header)
	{
		while ( activepairs.size() )
		{
			if ( activepairs.front().incnt == activepairs.front().outcnt )
			{
				if ( activepairs.front().getExpunge() )
				{
					expungeActiveFrontPairs(pairREC,header);			
				}
				else
				{
					finishActiveFrontPairs();
				}
			}
			else
			{
				// std::cerr << "WARNING: expunge on flush (this should not happen)" << std::endl;
				expungeActiveFrontPairs(pairREC,header);
			}
		}
	}

	void flushFrags(::libmaus::bambam::ReadEndsContainer * fragREC, ::libmaus::bambam::BamHeader const & header)
	{
		while ( activefrags.size() )
		{
			if ( activefrags.front().incnt == activefrags.front().outcnt )
			{
				if ( activefrags.front().getExpunge() )
				{
					expungeActiveFrontFrags(fragREC,header);			
				}
				else
				{
					finishActiveFrontFrags();
				}
			}
			else
			{
				// std::cerr << "WARNING: expunge on flush (this should not happen)" << std::endl;
				expungeActiveFrontFrags(fragREC,header);
			}
		}
	}
	
	/**
	 * flush lists
	 **/
	void flush(
		::libmaus::bambam::ReadEndsContainer * pairREC,
		::libmaus::bambam::ReadEndsContainer * fragREC,
		::libmaus::bambam::BamHeader const & header
	)
	{
		flushPairs(pairREC,header);
		flushFrags(fragREC,header);
	}
	
	static bool isMatchingPair(
		::libmaus::bambam::BamAlignment const & A,
		::libmaus::bambam::BamAlignment const & B
	)
	{
		bool ok = true;
		
		ok = ok && A.isPaired();
		ok = ok && B.isPaired();

		int const a1 = A.isRead1();
		int const a2 = A.isRead2();
		int const b1 = B.isRead1();
		int const b2 = B.isRead2();
		
		ok = ok && (a1+b1 == 1);
		ok = ok && (a2+b2 == 1);
		
		ok = ok && (A.getRefID() == B.getNextRefID());
		ok = ok && (B.getRefID() == A.getNextRefID());
		ok = ok && (A.getPos() == B.getNextPos());
		ok = ok && (B.getPos() == A.getNextPos());
		
		return ok;
	}
	
	/**
	 * add a pair
	 **/
	void addAlignmentPair(
		::libmaus::bambam::BamAlignment const & A, ::libmaus::bambam::BamAlignment const & B,
		::libmaus::bambam::ReadEndsContainer * pairREC,
		::libmaus::bambam::BamHeader const & header,
		uint64_t tagid
	)
	{
		active_count_type const bkey(B.getRefID(),B.getCoordinate(),0,0);
		
		bool done = false;
		
		while ( ! done )
		{
			if ( isActivePair(B) )
			{
				// expunge front element
				if ( APFLpairs.empty() )
				{
					assert ( activepairs.size() );
				
					// std::cerr << "Expunging " << activepairs.front() << " because free list is empty." << std::endl;	
					expungeActiveFrontPairs(pairREC,header);
					
					// check if this made any finished elements visible at the
					// front of the queue
					checkFinishedPairs(pairREC,bamheader);
				}
				//
				else
				{
					// find active_count_type object
					active_count_type const bkey(B.getRefID(),B.getCoordinate(),0,0);
					std::deque<active_count_type>::iterator const ita = std::lower_bound(activepairs.begin(),activepairs.end(),bkey);
					
					bool const activeok = 
						(ita != activepairs.end())
						&&
						(ita->refid == bkey.refid && ita->coordinate == bkey.coordinate);
						
					if ( ! activeok )
					{
						if ( isMatchingPair(A,B) )
						{
							libmaus::exception::LibMausException se;
							se.getStream() 
								<< "Unable to find active object for B=\n" <<  B.formatAlignment(header) << '\n' 
								<< "mate A=\n" << A.formatAlignment(header) << '\n'
								<< "refid=" << B.getRefID() << " coordinate=" << B.getCoordinate() << std::endl;
							se.finish();
							throw se;
						}
						else
						{
							#if 0
							std::cerr 
								<< "Non matching pair\n"
								<< A.formatAlignment(header) << "\n"
								<< B.formatAlignment(header) << "\n";
							#endif
							expungeUntilPairs(B,pairREC,header);
							pairREC->putPair(A,B,header,tagid);
							excntpairs += 1;
						}
					}
					else
					{
						assert ( activeok );
						
						// copy alignments
						// BamAlignmentPairListNode * ptr = APFLpairs.get();
						active_count_type::list_node_type * ptr = APFLpairs.get();
						#if defined(POS_READ_ENDS)
						libmaus::bambam::ReadEndsBase::fillFragPair(A,B,header,ptr->A,tagid);
						#else
						ptr->A[0].copyFrom(A);
						ptr->A[1].copyFrom(B);
						#endif
						ita->addAlignmentPair(ptr);
						ita->incOut();
					}
										
					// done inserting this one
					done = true;
				}
			}
			else
			{
				// interval was already handled and this read pair
				// is too late, handle it by writing it out
				pairREC->putPair(A,B,header,tagid);
				excntpairs += 1;
				done = true;
			}
		}		
	}

	/**
	 * add a pair
	 **/
	void addAlignmentFrag(
		::libmaus::bambam::BamAlignment const & B,
		::libmaus::bambam::ReadEndsContainer * fragREC,
		::libmaus::bambam::BamHeader const & header,
		uint64_t const tagid
	)
	{
		if ( B.getLseq() > maxreadlength )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "BamAlignmentInputPositionCallbackDupMark::addAlignmentFrag(): maximum allowed read length is " <<
				maxreadlength << " but input contains read of length " << B.getLseq() << "." << std::endl
				<< "Please set the maxreadlength parameter to a sufficiently high value." << std::endl;
			se.finish();
			throw se;
		}
	
		active_count_type const bkey(B.getRefID(),B.getCoordinate(),0,0);
		
		bool done = false;
		
		while ( ! done )
		{
			if ( isActiveFrag(B) )
			{
				// expunge front element
				if ( APFLfrags.empty() )
				{
					assert ( activefrags.size() );
				
					// std::cerr << "Expunging " << activefrags.front() << " because free list is empty." << std::endl;	
					expungeActiveFrontFrags(fragREC,header);
					
					// check if this made any finished elements visible at the
					// front of the queue
					checkFinishedFrags(fragREC,bamheader);
				}
				//
				else
				{
					// find active_count_type object
					active_count_type const bkey(B.getRefID(),B.getCoordinate(),0,0);
					std::deque<active_count_type>::iterator const ita = std::lower_bound(activefrags.begin(),activefrags.end(),bkey);
					assert ( ita != activefrags.end() );
					assert ( ita->refid == bkey.refid && ita->coordinate == bkey.coordinate );
					
					// copy alignments
					active_count_type::list_node_type * ptr = APFLfrags.get();
					ptr->A.reset();
					libmaus::bambam::ReadEndsBase::fillFrag(B,header,ptr->A,tagid);
					ita->addAlignmentPair(ptr);
					ita->incOut();

					// done inserting this one
					done = true;
				}
			}
			else
			{
				// interval was already handled and this frag
				// is too late, handle it by writing it out
				fragREC->putFrag(B,header,tagid);
				excntfrags += 1;
				done = true;
			}
		}		
	}
	
	void checkFinishedPairs(
		::libmaus::bambam::ReadEndsContainer * pairREC,
		::libmaus::bambam::BamHeader const & header
	)
	{
		// check for finished pair intervals
		while ( 
			activepairs.size()
			&&
			// input position is beyond end of activepairs front interval
			(
				position.first > activepairs.front().refid
				||
				(
					position.first == activepairs.front().refid &&
					position.second > activepairs.front().coordinate 
				)
			)
			&&
			// we have seen all pairs in the interval
			(
				activepairs.front().outcnt == activepairs.front().incnt
			)
		)
		{
			if ( activepairs.front().getExpunge() )
			{
				expungeActiveFrontPairs(pairREC,header);			
			}
			else
			{
				finishActiveFrontPairs();
			}
		}	
	}

	void checkFinishedFrags(
		::libmaus::bambam::ReadEndsContainer * fragREC,
		::libmaus::bambam::BamHeader const & header
	)
	{
	
		// check for finished pair intervals
		while ( 
			activefrags.size()
			&&
			// input position is beyond end of activefrags front interval
			(
				position.first > activefrags.front().refid
				||
				(
					position.first == activefrags.front().refid &&
					position.second > activefrags.front().coordinate + maxreadlength
				)
			)
			&&
			// we have seen all frags in the interval
			(
				activefrags.front().outcnt == activefrags.front().incnt
			)
		)
		{
			if ( activefrags.front().getExpunge() )
			{
				expungeActiveFrontFrags(fragREC,header);			
			}
			else
			{
				finishActiveFrontFrags();
			}
		}	
	}
};

struct PositionTrackCallback : 
	public ::libmaus::bambam::CollatingBamDecoderAlignmentInputCallback,
	public BamAlignmentInputPositionCallbackDupMark
{
	PositionTrackCallback(libmaus::bambam::BamHeader const & bamheader, uint64_t const freelistsize = getDefaultFreeListSize()) 
	: BamAlignmentInputPositionCallbackDupMark(bamheader,freelistsize) 
	{
	
	}
	virtual ~PositionTrackCallback() {}
	
	void operator()(::libmaus::bambam::BamAlignment const & A)
	{
		updatePosition(A);
	}
};

static int markDuplicates(::libmaus::util::ArgInfo const & arginfo)
{
	libmaus::timing::RealTimeClock globrtc; globrtc.start();

	::libmaus::util::TempFileRemovalContainer::setup();

	if ( (!(arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != ""))) && isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "refusing to read compressed data from terminal. please use I=<filename> or redirect standard input to a file" << std::endl;
		se.finish();
		throw se;
	}
	
	if ( (!(arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != ""))) && isatty(STDOUT_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "refusing to write compressed data to terminal. please use O=<filename> or redirect standard output to a file" << std::endl;
		se.finish();
		throw se;		
	}
	
	// logarithm of collation hash table size
	unsigned int const colhashbits = arginfo.getValue<unsigned int>("colhashbits",getDefaultColHashBits());
	// length of collation output list
	uint64_t const collistsize = arginfo.getValueUnsignedNumeric<uint64_t>("collistsize",getDefaultColListSize());
	// buffer size for fragment and pair data
	uint64_t const fragbufsize = arginfo.getValueUnsignedNumeric<uint64_t>("fragbufsize",getDefaultFragBufSize());
	// print verbosity messages
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	// rewritten file should be in bam format, if input is given via stdin
	unsigned int const rewritebam = arginfo.getValue<unsigned int>("rewritebam",getDefaultRewriteBam());

	// buffer size for fragment and pair data
	uint64_t const maxreadlength = arginfo.getValueUnsignedNumeric<uint64_t>("maxreadlength",getDefaultMaxReadLength());

	enum tag_type_enum
	{
		tag_type_none,
		tag_type_string,
		tag_type_nucleotide
	};

	// tag field
	bool const havetag = arginfo.hasArg("tag");
	std::string const tag = arginfo.getUnparsedValue("tag","no tag");
	libmaus::trie::SimpleTrie::unique_ptr_type Ptagtrie;

	if ( havetag && (tag.size() != 2 || (!isalpha(tag[0])) || (!isalnum(tag[1])) ) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "tag " << tag << " is invalid" << std::endl;
		se.finish();
		throw se;			
	}
	
	if ( havetag )
	{
		libmaus::trie::SimpleTrie::unique_ptr_type Ttagtrie(new libmaus::trie::SimpleTrie);
		Ptagtrie = UNIQUE_PTR_MOVE(Ttagtrie);

		// allocate tag id 0 for empty tag
		uint8_t const * p = 0;
		Ptagtrie->insert(p,p);
	}

	// nucl tag field
	bool const havenucltag = arginfo.hasArg("nucltag");
	std::string const nucltag = arginfo.getUnparsedValue("nucltag","no tag");

	if ( havenucltag && (nucltag.size() != 2 || (!isalpha(nucltag[0])) || (!isalnum(nucltag[1])) ) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "nucltag " << tag << " is invalid" << std::endl;
		se.finish();
		throw se;			
	}
	
	if ( havetag && havenucltag )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "tag and nucltag are mutually exclusive" << std::endl;
		se.finish();
		throw se;					
	}
	
	tag_type_enum tag_type;
	
	if ( havetag )
		tag_type = tag_type_string;
	else if ( havenucltag )
		tag_type = tag_type_nucleotide;
	else
		tag_type = tag_type_none;

	char const * ctag = havetag ? tag.c_str() : 0;
	char const * cnucltag = havenucltag ? nucltag.c_str() : 0;
	char const * tag1 = 0;
	char const * tag2 = 0;
	libmaus::autoarray::AutoArray<char> tagbuffer;
	uint64_t taglen = 0;
	uint64_t tagid = 0;
	libmaus::fastx::FastATwoBitTable const FATBT;

	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilename = tmpfilenamebase + "_bamcollate";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const tmpfilenamereadfrags = tmpfilenamebase + "_readfrags";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenamereadfrags);
	std::string const tmpfilenamereadpairs = tmpfilenamebase + "_readpairs";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenamereadpairs);
	std::string const tmpfilesnappyreads = tmpfilenamebase + "_alignments";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilesnappyreads);
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string const tmpfiledupset = tmpfilenamebase + "_dupset";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfiledupset);

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const rewritebamlevel = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("rewritebamlevel",getDefaultRewriteBamLevel()));
		
	if ( verbose )
		std::cerr << "[V] output compression level " << level << std::endl;

	::libmaus::timing::RealTimeClock fragrtc; fragrtc.start();

	libmaus::bambam::BamAlignmentInputCallbackSnappy<BamAlignmentInputPositionCallbackDupMark>::unique_ptr_type SRC;
	libmaus::bambam::BamAlignmentInputCallbackBam<BamAlignmentInputPositionCallbackDupMark>::unique_ptr_type BWR;
	BamAlignmentInputPositionCallbackDupMark * PTI = 0;
	::libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	libmaus::aio::PosixFdInputStream::unique_ptr_type PFIS;
	libmaus::aio::CheckedOutputStream::unique_ptr_type copybamstr;

	typedef ::libmaus::bambam::BamCircularHashCollatingBamDecoder col_type;
	typedef ::libmaus::bambam::BamParallelCircularHashCollatingBamDecoder par_col_type;
	typedef ::libmaus::bambam::CircularHashCollatingBamDecoder col_base_type;
	typedef ::libmaus::bambam::BamMergeCoordinateCircularHashCollatingBamDecoder merge_col_type;
	typedef col_base_type::unique_ptr_type col_base_ptr_type;	
	col_base_ptr_type CBD;

	uint64_t const markthreads = 
		std::max(static_cast<uint64_t>(1),arginfo.getValue<uint64_t>("markthreads",getDefaultMarkThreads()));

	uint64_t const colexcludeflags =
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY;
	
	if ( arginfo.getPairCount("I") > 1 )
	{
		std::vector<std::string> const inputfilenames = arginfo.getPairValues("I");
		col_base_ptr_type tCBD(new merge_col_type(inputfilenames,tmpfilename,colexcludeflags,true,colhashbits,collistsize));
		CBD = UNIQUE_PTR_MOVE(tCBD);
	}
	// if we are reading the input from a file
	else if ( arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != "") )
	{
		std::string const inputfilename = arginfo.getValue<std::string>("I","I");
		uint64_t const inputbuffersize = arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",getDefaultInputBufferSize());
		libmaus::aio::PosixFdInputStream::unique_ptr_type tPFIS(new libmaus::aio::PosixFdInputStream(inputfilename,inputbuffersize,0));
		PFIS = UNIQUE_PTR_MOVE(tPFIS);
		
		if ( markthreads > 1 )
		{
			col_base_ptr_type tCBD(new par_col_type(
                                *PFIS,
                                markthreads,
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));
			CBD = UNIQUE_PTR_MOVE(tCBD);
		}
		else
		{
			#if 0
			col_base_ptr_type tCBD(new col_type(
                                *CIS,
                                // numthreads
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));
                        #else	
			col_base_ptr_type tCBD(new col_type(
                                *PFIS,
                                // numthreads
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));	
			#endif
			CBD = UNIQUE_PTR_MOVE(tCBD);
		}
	}
	// not a file, we are reading from standard input
	else
	{

		// rewrite to bam
		if ( rewritebam )
		{
			if ( rewritebam > 1 )
			{
				libmaus::aio::CheckedOutputStream::unique_ptr_type tcopybamstr(new libmaus::aio::CheckedOutputStream(tmpfilesnappyreads));
				copybamstr = UNIQUE_PTR_MOVE(tcopybamstr);

				if ( markthreads > 1 )
				{
					col_base_ptr_type tCBD(new par_col_type(std::cin,
                                                *copybamstr,
                                                markthreads,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}
				else
				{
					col_base_ptr_type tCBD(new col_type(std::cin,
                                                *copybamstr,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}

				if ( verbose )
					std::cerr << "[V] Copying bam compressed alignments to file " << tmpfilesnappyreads << std::endl;
			}
			else
			{
				if ( markthreads > 1 )
				{
					col_base_ptr_type tCBD(new par_col_type(std::cin,markthreads,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}
				else
				{
					col_base_ptr_type tCBD(new col_type(std::cin,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}

				// rewrite file and mark duplicates
				libmaus::bambam::BamAlignmentInputCallbackBam<BamAlignmentInputPositionCallbackDupMark>::unique_ptr_type tBWR(new libmaus::bambam::BamAlignmentInputCallbackBam<BamAlignmentInputPositionCallbackDupMark>(tmpfilesnappyreads,CBD->getHeader(),rewritebamlevel));
				BWR = UNIQUE_PTR_MOVE(tBWR);
				CBD->setInputCallback(BWR.get());

				if ( verbose )
					std::cerr << "[V] Writing bam compressed alignments to file " << tmpfilesnappyreads << std::endl;
			}
		}
		else
		{
			col_base_ptr_type tCBD(new col_type(std::cin,
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));
			CBD = UNIQUE_PTR_MOVE(tCBD);

			libmaus::bambam::BamAlignmentInputCallbackSnappy<BamAlignmentInputPositionCallbackDupMark>::unique_ptr_type 
				tSRC(new libmaus::bambam::BamAlignmentInputCallbackSnappy<BamAlignmentInputPositionCallbackDupMark>(tmpfilesnappyreads,CBD->getHeader()));
			SRC = UNIQUE_PTR_MOVE(tSRC);
			CBD->setInputCallback(SRC.get());
			if ( verbose )
				std::cerr << "[V] Writing snappy compressed alignments to file " << tmpfilesnappyreads << std::endl;
		}
	}

	::libmaus::bambam::BamHeader const bamheader = CBD->getHeader();
	
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CBD->disableValidation();

	uint64_t const trackfreelistsize = arginfo.getValueUnsignedNumeric<uint64_t>("trackfreelistsize",PositionTrackCallback::getDefaultFreeListSize());
	PositionTrackCallback PTC(bamheader,trackfreelistsize);
	
	if ( SRC )
		PTI = SRC.get();
	else if ( BWR )
		PTI = BWR.get();
	else
	{
		CBD->setInputCallback(&PTC);
		PTI = & PTC;
	}
	
	// CBD->setPrimaryExpungeCallback(&PTC);
	
	typedef col_base_type::alignment_ptr_type alignment_ptr_type;
	std::pair<alignment_ptr_type,alignment_ptr_type> P;
	uint64_t const mod = arginfo.getValue<unsigned int>("mod",getDefaultMod()); // modulus for verbosity
	uint64_t fragcnt = 0; // mapped fragments
	uint64_t paircnt = 0; // mapped pairs
	uint64_t lastproc = 0; // printed at last fragment count

	// bool const copyAlignments = false;
	#if defined(MARKDUPLICATECOPYALIGNMENTS)
	bool const copyAlignments = true;
	#else
	bool const copyAlignments = false;
	#endif
	::libmaus::bambam::ReadEndsContainer::unique_ptr_type fragREC(new ::libmaus::bambam::ReadEndsContainer(fragbufsize,tmpfilenamereadfrags,copyAlignments)); // fragment container
	::libmaus::bambam::ReadEndsContainer::unique_ptr_type pairREC(new ::libmaus::bambam::ReadEndsContainer(fragbufsize,tmpfilenamereadpairs,copyAlignments)); // pair container
	
	int64_t maxrank = -1; // maximal appearing rank
	uint64_t als = 0; // number of processed alignments (= mapped+unmapped fragments)
	std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> metrics;

	::libmaus::timing::RealTimeClock rtc; rtc.start(); // clock

	// #define UNPAIREDDEBUG

	#if defined(UNPAIREDDEBUG)
	::libmaus::bambam::BamFormatAuxiliary bamauxiliary;
	#endif
	
	libmaus::timing::RealTimeClock readinrtc; readinrtc.start();

	::libmaus::bambam::DupSetCallbackStream DSCV(tmpfiledupset,metrics);
	// DupSetCallbackSet DSCV(metrics);
	
	PTI->setDupSetCallback(&DSCV);
	PTI->setMaxReadLength(maxreadlength);

	
	while ( CBD->tryPair(P) )
	{
		#if 0
		std::cerr 
			<< PTI->getPosition().first << ","
			<< PTI->getPosition().second << std::endl;
		#endif
	
		assert ( P.first || P.second );
		uint64_t const lib = 
			P.first 
			? 
			P.first->getLibraryId(bamheader) 
			:
			P.second->getLibraryId(bamheader) 				
			;
		::libmaus::bambam::DuplicationMetrics & met = metrics[lib];
		
		if ( P.first )
		{
			maxrank = std::max(maxrank,P.first->getRank());
			als++;
			
			if ( P.first->isUnmap() ) 
			{
				++met.unmapped;
			}                                    
			else if ( (!P.first->isPaired()) || P.first->isMateUnmap() )
			{
				#if defined(UNPAIREDDEBUG)
				std::cerr << "[D]\t1\t" << P.first->formatAlignment(bamheader,bamauxiliary) << std::endl;
				#endif
				met.unpaired++;
			}
		}
		if ( P.second )
		{
			maxrank = std::max(maxrank,P.second->getRank());
			als++;

			if ( P.second->isUnmap() )
			{
				++met.unmapped;
			}
			else if ( (!P.second->isPaired()) || P.second->isMateUnmap() )
			{
				#if defined(UNPAIREDDEBUG)
				std::cerr << "[D]\t2\t" << P.first->formatAlignment(bamheader,bamauxiliary) << std::endl;
				#endif
				met.unpaired++;
			}
		}
		
		switch ( tag_type )
		{
			case tag_type_string:
			{
				// length of tags for read1 and read2
				uint64_t l1 = 0, l2 = 0;
				
				// aux lookup for read1
				if ( P.first )
				{
					tag1 = P.first->getAuxString(ctag);
					l1 = tag1 ? strlen(tag1) : 0;
				}
				// aux lookup for read2
				if ( P.second )
				{
					tag2 = P.second->getAuxString(ctag);
					l2 = tag2 ? strlen(tag2) : 0;
				}
				
				// length of concatenated tag
				taglen = l1 + l2 + 2;
				// expand buffer if necessary
				if ( taglen > tagbuffer.size() )
					tagbuffer = libmaus::autoarray::AutoArray<char>(taglen,false);

				// concatenate tags
				char * outptr = tagbuffer.begin();

				memcpy(outptr,tag1,l1);
				outptr += l1;
				*(outptr++) = 0;

				memcpy(outptr,tag2,l2);
				outptr += l2;
				*(outptr++) = 0;

				assert ( outptr - tagbuffer.begin() == static_cast<ptrdiff_t>(taglen) );

				// look up tag id			
				tagid = Ptagtrie->insert(
					tagbuffer.begin(),
					outptr
				);

				break;
			}
			case tag_type_nucleotide:
			{
				// aux lookup for read1
				if ( P.first )
					tag1 = P.first->getAuxString(cnucltag);
				// aux lookup for read2
				if ( P.second )
					tag2 = P.second->getAuxString(cnucltag);

				tagid = (FATBT(tag1) << 32) | FATBT(tag2);
				
				break;
			}
			default:
			{
				tagid = 0;
				break;
			}
		}

		// we are not interested in unmapped reads, ignore them
		if ( P.first && P.first->isUnmap() )
		{
			P.first = 0;
		}
		if ( P.second && P.second->isUnmap() )
		{
			P.second = 0;
		}
			
		if ( P.first && P.second )
		{
			#if defined(DEBUG)
			std::cerr << "[V] Got pair for name " << P.first->getName() 
				<< "," << P.first->getCoordinate()
				<< "," << P.first->getFlagsS()
				<< "," << P.second->getCoordinate()
				<< "," << P.second->getFlagsS()
				<< std::endl;
			#endif
		
			met.readpairsexamined++;
		
			assert ( ! P.first->isUnmap() );
			assert ( ! P.second->isUnmap() );
			
			// swap reads if necessary so P.first is left of P.second
			// in terms of coordinates
			if ( 
				// second has higher ref id
				(
					P.second->getRefID() > P.first->getRefID()	
				)
				||
				// same refid, second has higher position
				(
					P.second->getRefID() == P.first->getRefID() 
					&&
					P.second->getCoordinate() > P.first->getCoordinate()
				)
				||
				// same refid, same position, first is read 1
				(
					P.second->getRefID() == P.first->getRefID() 
					&&
					P.second->getCoordinate() == P.first->getCoordinate()
					&&
					P.first->isRead1()
				)
			)
			{
			
			}
			else
			{
				std::swap(P.first,P.second);
			}
			
			// simple one, register it
			if (  BamAlignmentInputPositionCallbackDupMark::isSimplePair(*(P.second)) )
			{
				PTI->addAlignmentPair(*(P.first),*(P.second),pairREC.get(),bamheader,tagid);
				PTI->checkFinishedPairs(pairREC.get(),bamheader);
			}
			// strangely mapped pairs, mark right coordinate as expunged
			else
			{
				PTI->setExpungePairs(*(P.second));
				PTI->strcntpairs += 1;
				pairREC->putPair(*(P.first),*(P.second),bamheader,tagid);
			}

			paircnt++;
		}
		
		if ( P.first )
		{
			PTI->addAlignmentFrag(*(P.first),fragREC.get(),bamheader,tagid);
			fragcnt++;
		}
		if ( P.second )
		{
			PTI->addAlignmentFrag(*(P.second),fragREC.get(),bamheader,tagid);
			fragcnt++;
		}	
		
		if ( verbose && fragcnt/mod != lastproc/mod )
		{
			std::cerr
				<< "[V] " 
				<< als << " als, "
				<< fragcnt << " mapped frags, " 
				<< paircnt << " mapped pairs, "
				<< fragcnt/rtc.getElapsedSeconds() << " frags/s "
				<< ::libmaus::util::MemUsage()
				<< " time "
				<< readinrtc.getElapsedSeconds()
				<< " total "
				<< fragrtc.formatTime(fragrtc.getElapsedSeconds())
				<< std::endl;
			readinrtc.start();
			lastproc = fragcnt;
		}		
	}

	PTI->flush(pairREC.get(),fragREC.get(),bamheader);
	std::cerr << "[D] " << "excntpairs=" << PTI->excntpairs << " fincntpairs=" << PTI->fincntpairs << " strcntpairs=" << PTI->strcntpairs << std::endl;
	uint64_t const diskpairs = PTI->excntpairs + PTI->strcntpairs;
	std::cerr << "[D] " << "excntfrags=" << PTI->excntfrags << " fincntfrags=" << PTI->fincntfrags << " strcntfrags=" << PTI->strcntfrags << std::endl;
	uint64_t const diskfrags = PTI->excntfrags + PTI->strcntfrags;
	
	assert ( PTI->excntpairs + PTI->fincntpairs + PTI->strcntpairs == paircnt );
	assert ( PTI->excntfrags + PTI->fincntfrags + PTI->strcntfrags == fragcnt );
	
	if ( copybamstr )
	{
		copybamstr->flush();
		copybamstr.reset();
	}

	uint64_t const numranks = CBD->getRank(); // maxrank+1;
	
	CBD.reset();
	CIS.reset();
	SRC.reset();
	BWR.reset();
			
	fragREC->flush();
	pairREC->flush();
	fragREC->releaseArray();
	pairREC->releaseArray();
	
	if ( verbose )
		std::cerr << "[V] fragment and pair data computed in time " << fragrtc.getElapsedSeconds() << " (" << fragrtc.formatTime(fragrtc.getElapsedSeconds()) << ")" << std::endl;

	#if 0
	uint64_t const numranks = maxrank+1;
	
	if ( numranks != als )
		std::cerr << "[D] numranks=" << numranks << " != als=" << als << std::endl;
	
	assert ( numranks == als );
	#endif

	if ( verbose )
		std::cerr
			<< "[V] " 
			<< numranks << " lines, "
			<< als << " als, "
			<< fragcnt << " mapped frags, " 
			<< paircnt << " mapped pairs, "
			<< fragcnt/rtc.getElapsedSeconds() << " frags/s "
			<< ::libmaus::util::MemUsage()
			<< std::endl;

	// DupSetCallbackVector DSCV(numranks,metrics);

	/*
	 * process fragment and pair data to determine which reads are to be marked as duplicates
	 */		
	::libmaus::bambam::ReadEnds nextfrag;
	std::vector< ::libmaus::bambam::ReadEnds > lfrags;
	uint64_t dupcnt = 0;

	if ( verbose )
		std::cerr << "[V] Checking pairs...";
	rtc.start();
	::libmaus::bambam::SortedFragDecoder::unique_ptr_type pairDec(pairREC->getDecoder());
	pairREC.reset();
	pairDec->getNext(lfrags);

	while ( pairDec->getNext(nextfrag) )
	{
		if ( ! libmaus::bambam::DupMarkBase::isDupPair(nextfrag,lfrags.front()) )
		{
			// dupcnt += markDuplicatePairs(lfrags.begin(),lfrags.end(),DSCV);
			dupcnt += libmaus::bambam::DupMarkBase::markDuplicatePairsVector(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	// dupcnt += markDuplicatePairs(lfrags.begin(),lfrags.end(),DSCV);
	dupcnt += libmaus::bambam::DupMarkBase::markDuplicatePairsVector(lfrags,DSCV);
	lfrags.resize(0);
	pairDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << (diskpairs)/rtc.getElapsedSeconds() << std::endl;

	if ( verbose )
		std::cerr << "[V] Checking single fragments...";
	rtc.start();
	::libmaus::bambam::SortedFragDecoder::unique_ptr_type fragDec(fragREC->getDecoder());
	fragREC.reset();
	fragDec->getNext(lfrags);
	while ( fragDec->getNext(nextfrag) )
	{
		if ( !libmaus::bambam::DupMarkBase::isDupFrag(nextfrag,lfrags.front()) )
		{
			dupcnt += libmaus::bambam::DupMarkBase::markDuplicateFrags(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	dupcnt += libmaus::bambam::DupMarkBase::markDuplicateFrags(lfrags,DSCV);
	lfrags.resize(0);
	fragDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << (diskfrags)/rtc.getElapsedSeconds() << std::endl;		

	DSCV.flush(numranks);

	if ( verbose )
		std::cerr << "[V] number of alignments marked as duplicates: " << DSCV.getNumDups() << " time " << fragrtc.getElapsedSeconds() << " (" << fragrtc.formatTime(fragrtc.getElapsedSeconds()) << ")" << std::endl;
	/*
	 * end of fragment processing
	 */

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
		ita->second.format(metricsstr, bamheader.getLibraryName(ita->first));
	
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

	/*
	 * mark the duplicates
	 */
	libmaus::bambam::DupMarkBase::markDuplicatesInFile(
		arginfo,verbose,bamheader,maxrank,mod,level,DSCV,tmpfilesnappyreads,rewritebam,tmpfileindex,
		getProgId(),
		std::string(PACKAGE_VERSION),
		getDefaultRmDup(),
		getDefaultMD5(),
		getDefaultIndex(),
		getDefaultMarkThreads()
	);
		
	if ( verbose )
		std::cerr << "[V] " << ::libmaus::util::MemUsage() << " " 
			<< globrtc.getElapsedSeconds() 
			<< " ("
			<< globrtc.formatTime(globrtc.getElapsedSeconds())
			<< ")"
			<< std::endl;
		
	return EXIT_SUCCESS;
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
				
				V.push_back ( std::pair<std::string,std::string> ( "I=<filename>", "input file, stdin if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<filename>", "output file, stdout if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "M=<filename>", "metrics file, stderr if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "markthreads=<["+::biobambam::Licensing::formatNumber(getDefaultMarkThreads())+"]>", "number of helper threads" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report (default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for each mod'th record/alignment" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebam=<["+::biobambam::Licensing::formatNumber(getDefaultRewriteBam())+"]>", "compression of temporary alignment file when input is via stdin (0=snappy,1=gzip/bam,2=copy)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebamlevel=<["+::biobambam::Licensing::formatNumber(getDefaultRewriteBamLevel())+"]>", std::string("compression setting for rewritten input file if rewritebam=1 (") + libmaus::bambam::BamBlockWriterBaseFactory::getLevelHelpText() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rmdup=<["+::biobambam::Licensing::formatNumber(getDefaultRmDup())+"]>", "remove duplicates (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhashbits=<["+::biobambam::Licensing::formatNumber(getDefaultColHashBits())+"]>", "log_2 of size of hash table used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collistsize=<["+::biobambam::Licensing::formatNumber(getDefaultColListSize())+"]>", "output list size for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "fragbufsize=<["+::biobambam::Licensing::formatNumber(getDefaultFragBufSize())+"]>", "size of each fragment/pair file buffer in bytes" ) );
				V.push_back ( std::pair<std::string,std::string> ( "maxreadlength=<["+::biobambam::Licensing::formatNumber(getDefaultMaxReadLength())+"]>", "maximum allowed read length" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "trackfreelistsize=<["+::biobambam::Licensing::formatNumber(PositionTrackCallback::getDefaultFreeListSize())+"]>", "tracking lists free pool size" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputbuffersize=<["+::biobambam::Licensing::formatNumber(getDefaultInputBufferSize())+"]>", "size of input buffer" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tag=<[a-zA-Z][a-zA-Z0-9]>", "aux field id for tag string extraction" ) );
				V.push_back ( std::pair<std::string,std::string> ( "nucltag=<[a-zA-Z][a-zA-Z0-9]>", "aux field id for nucleotide tag extraction" ) );
				V.push_back ( std::pair<std::string,std::string> ( "D=<filename>", "duplicates output file if rmdup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupmd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for duplicates output file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupmd5filename=<filename>", "file name for md5 check sum of dup file (default: extend duplicates output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupindex=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index for duplicates file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupindexfilename=<filename>", "file name for BAM index file for duplicates file (default: extend duplicates output file name)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return markDuplicates(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

