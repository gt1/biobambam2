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
#include <libmaus/bambam/BamParallelRewrite.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/DuplicationMetrics.hpp>
#include <libmaus/bambam/OpticalComparator.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/bambam/ReadEndsContainer.hpp>
#include <libmaus/bambam/SortedFragDecoder.hpp>
#include <libmaus/bitio/BitVector.hpp>
#include <libmaus/lz/BgzfInflateDeflateParallel.hpp>
#include <libmaus/lz/BgzfRecode.hpp>
#include <libmaus/math/iabs.hpp>
#include <libmaus/math/numbits.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/ContainerGetObject.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <biobambam/Licensing.hpp>

// static std::string formatNumber(int64_t const n) { std::ostringstream ostr; ostr << n; return ostr.str(); }
static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1048576;  }
static bool getDefaultRewriteBam() { return 0; }
static int getDefaultRewriteBamLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultColHashBits() { return 20; }
static uint64_t getDefaultColListSize() { return 32*1024*1024; }
static uint64_t getDefaultFragBufSize() { return 48*1024*1024; }
static uint64_t getDefaultMarkThreads() { return 1; }
static bool getDefaultRmDup() { return 0; }

struct DupSetCallback
{
	virtual void operator()(::libmaus::bambam::ReadEnds const & A) = 0;
	virtual uint64_t getNumDups() const = 0;
	virtual void addOpticalDuplicates(uint64_t const libid, uint64_t const count) = 0;
	virtual bool isMarked(uint64_t const i) const = 0;
	virtual void flush(uint64_t const n) = 0;
};

struct DupSetCallbackSet : public DupSetCallback
{
	std::set<uint64_t> S;

	virtual void operator()(::libmaus::bambam::ReadEnds const & A)
	{
		S.insert(A.getRead1IndexInFile());
		
		if ( A.isPaired() )
			S.insert(A.getRead2IndexInFile());
	}
	virtual uint64_t getNumDups() const
	{
		return S.size();
	}
	virtual void addOpticalDuplicates(uint64_t const, uint64_t const)
	{
	
	}
	virtual bool isMarked(uint64_t const i) const
	{
		return S.find(i) != S.end();
	}
	virtual void flush(uint64_t const)
	{
	
	}
};

struct DupSetCallbackVector : public DupSetCallback
{
	::libmaus::bitio::BitVector B;

	std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> & metrics;

	DupSetCallbackVector(
		uint64_t const n,
		std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> & rmetrics
	) : B(n), metrics(rmetrics) /* unpairedreadduplicates(), readpairduplicates(), metrics(rmetrics) */ {}
	
	void operator()(::libmaus::bambam::ReadEnds const & A)
	{
		B.set(A.getRead1IndexInFile(),true);
		
		if ( A.isPaired() )
		{
			B.set(A.getRead2IndexInFile(),true);
			metrics[A.getLibraryId()].readpairduplicates++;
		}
		else
		{
			metrics[A.getLibraryId()].unpairedreadduplicates++;
		}
	}	

	uint64_t getNumDups() const
	{
		uint64_t dups = 0;
		for ( uint64_t i = 0; i < B.size(); ++i )
			if ( B.get(i) )
				dups++;
		
		return dups;
	}
	void addOpticalDuplicates(uint64_t const libid, uint64_t const count)
	{
		metrics[libid].opticalduplicates += count;
	}

	bool isMarked(uint64_t const i) const
	{
		return B[i];
	}
	
	void flush(uint64_t const)
	{
	
	}
};

struct DupSetCallbackStream : public DupSetCallback
{
	std::string const filename;
	libmaus::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type SGO;
	std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> & metrics;
	uint64_t numdup;
	::libmaus::bitio::BitVector::unique_ptr_type B;

	DupSetCallbackStream(
		std::string const & rfilename,
		std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> & rmetrics
	) : filename(rfilename), SGO(new libmaus::aio::SynchronousGenericOutput<uint64_t>(filename,8*1024)), metrics(rmetrics), numdup(0) /* unpairedreadduplicates(), readpairduplicates(), metrics(rmetrics) */ {}
	
	void operator()(::libmaus::bambam::ReadEnds const & A)
	{
		SGO->put(A.getRead1IndexInFile());
		numdup++;
		
		if ( A.isPaired() )
		{
			SGO->put(A.getRead2IndexInFile());
			numdup++;
			metrics[A.getLibraryId()].readpairduplicates++;
		}
		else
		{
			metrics[A.getLibraryId()].unpairedreadduplicates++;
		}
	}	

	uint64_t getNumDups() const
	{
		return numdup;
	}
	void addOpticalDuplicates(uint64_t const libid, uint64_t const count)
	{
		metrics[libid].opticalduplicates += count;
	}

	bool isMarked(uint64_t const i) const
	{
		return (*B)[i];
	}
	
	void flush(uint64_t const n)
	{
		SGO->flush();
		SGO.reset();
		
		::libmaus::bitio::BitVector::unique_ptr_type tB(new ::libmaus::bitio::BitVector(n));
		B = UNIQUE_PTR_MOVE(tB);
		for ( uint64_t i = 0; i < n; ++i )
			B->set(i,false);
			
		libmaus::aio::SynchronousGenericInput<uint64_t> SGI(filename,8*1024);
		uint64_t v;
		numdup = 0;
		while ( SGI.getNext(v) )
		{
			if ( ! B->get(v) )
				numdup++;
			B->set(v,true);			
		}
	}
};

// #define DEBUG

#if defined(DEBUG)
#define MARKDUPLICATEPAIRSDEBUG
#define MARKDUPLICATEFRAGSDEBUG
#define MARKDUPLICATECOPYALIGNMENTS
#endif

static uint64_t markDuplicatePairsRef(
	std::vector< ::libmaus::bambam::ReadEnds > & lfrags, 
	DupSetCallback & DSC,
	unsigned int const optminpixeldif = 100
)
{
	if ( lfrags.size() > 1 )
	{
		#if defined(MARKDUPLICATEPAIRSDEBUG)
		std::cerr << "[V] --- checking for duplicate pairs ---" << std::endl;
		for ( uint64_t i = 0; i < lfrags.size(); ++i )
			std::cerr << "[V] " << lfrags[i] << std::endl;
		#endif
	
		uint64_t maxscore = lfrags[0].getScore();
		uint64_t maxindex = 0;
		
		for ( uint64_t i = 1 ; i < lfrags.size(); ++i )
			if ( lfrags[i].getScore() > maxscore )
			{
				maxscore = lfrags[i].getScore();
				maxindex = i;
			}

		for ( uint64_t i = 0; i < lfrags.size(); ++i )
			if ( i != maxindex )
				DSC(lfrags[i]);
	
		// check for optical duplicates
		std::sort ( lfrags.begin(), lfrags.end(), ::libmaus::bambam::OpticalComparator() );
		
		for ( uint64_t low = 0; low < lfrags.size(); )
		{
			uint64_t high = low+1;
			
			// search top end of tile
			while ( 
				high < lfrags.size() && 
				lfrags[high].getReadGroup() == lfrags[low].getReadGroup() &&
				lfrags[high].getTile() == lfrags[low].getTile() )
			{
				high++;
			}
			
			if ( high-low > 1 && lfrags[low].getTile() )
			{
				#if defined(DEBUG)
				std::cerr << "[D] Range " << high-low << " for " << lfrags[low] << std::endl;
				#endif
			
				std::vector<bool> opt(high-low,false);
				bool haveoptdup = false;
				
				for ( uint64_t i = low; i+1 < high; ++i )
				{
					for ( 
						uint64_t j = i+1; 
						j < high && lfrags[j].getX() - lfrags[low].getX() <= optminpixeldif;
						++j 
					)
						if ( 
							::libmaus::math::iabs(
								static_cast<int64_t>(lfrags[i].getY())
								-
								static_cast<int64_t>(lfrags[j].getY())
							)
							<= optminpixeldif
						)
					{	
						opt [ j - low ] = true;
						haveoptdup = true;
					}
				}
				
				if ( haveoptdup )
				{
					unsigned int const lib = lfrags[low].getLibraryId();
					uint64_t numopt = 0;
					for ( uint64_t i = 0; i < opt.size(); ++i )
						if ( opt[i] )
							numopt++;
							
					DSC.addOpticalDuplicates(lib,numopt);	
				}
			}
			
			low = high;
		}	
	}
	else
	{
		#if defined(MARKDUPLICATEPAIRSDEBUG)
		std::cerr << "[V] --- singleton pair set ---" << std::endl;
		for ( uint64_t i = 0; i < lfrags.size(); ++i )
			std::cerr << "[V] " << lfrags[i] << std::endl;
		#endif
	
	}
	
	// all but one are duplicates
	return lfrags.size() ? 2*(lfrags.size() - 1) : 0;
}


template<typename iterator>
static uint64_t markDuplicatePairs(
	iterator const lfrags_a,
	iterator const lfrags_e,
	DupSetCallback & DSC,
	unsigned int const optminpixeldif = 100
)
{
	if ( lfrags_e-lfrags_a > 1 )
	{
		#if defined(MARKDUPLICATEPAIRSDEBUG)
		std::cerr << "[V] --- checking for duplicate pairs ---" << std::endl;
		for ( iterator lfrags_c = lfrags_a; lfrags_c != lfrags_e; ++lfrags_c )
			std::cerr << "[V] " << (*lfrags_c) << std::endl;
		#endif
	
		uint64_t maxscore = lfrags_a->getScore();
		
		iterator lfrags_m = lfrags_a;
		for ( iterator lfrags_c = lfrags_a+1; lfrags_c != lfrags_e; ++lfrags_c )
			if ( lfrags_c->getScore() > maxscore )
			{
				maxscore = lfrags_c->getScore();
				lfrags_m = lfrags_c;
			}

		for ( iterator lfrags_c = lfrags_a; lfrags_c != lfrags_e; ++lfrags_c )
			if ( lfrags_c != lfrags_m )
				DSC(*lfrags_c);
	
		// check for optical duplicates
		std::sort ( lfrags_a, lfrags_e, ::libmaus::bambam::OpticalComparator() );
		
		for ( iterator low = lfrags_a; low != lfrags_e; )
		{
			iterator high = low+1;
			
			// search top end of tile
			while ( 
				high != lfrags_e && 
				high->getReadGroup() == low->getReadGroup() &&
				high->getTile() == low->getTile() )
			{
				++high;
			}
			
			if ( high-low > 1 && low->getTile() )
			{
				#if defined(DEBUG)
				std::cerr << "[D] Range " << high-low << " for " << lfrags[low] << std::endl;
				#endif
			
				std::vector<bool> opt(high-low,false);
				bool haveoptdup = false;
				
				for ( iterator i = low; i+1 != high; ++i )
				{
					for ( 
						iterator j = i+1; 
						j != high && j->getX() - low->getX() <= optminpixeldif;
						++j 
					)
						if ( 
							::libmaus::math::iabs(
								static_cast<int64_t>(i->getY())
								-
								static_cast<int64_t>(j->getY())
							)
							<= optminpixeldif
						)
					{	
						opt [ j - low ] = true;
						haveoptdup = true;
					}
				}
				
				if ( haveoptdup )
				{
					unsigned int const lib = low->getLibraryId();
					uint64_t numopt = 0;
					for ( uint64_t i = 0; i < opt.size(); ++i )
						if ( opt[i] )
							numopt++;
							
					DSC.addOpticalDuplicates(lib,numopt);	
				}
			}
			
			low = high;
		}	
	}
	else
	{
		#if defined(MARKDUPLICATEPAIRSDEBUG)
		std::cerr << "[V] --- singleton pair set ---" << std::endl;
		for ( iterator i = lfrags_a; i != lfrags_e; ++i )
			std::cerr << "[V] " << (*i) << std::endl;
		#endif
	
	}
	
	uint64_t const lfragssize = lfrags_e-lfrags_a;
	// all but one are duplicates
	return lfragssize ? 2*(lfragssize - 1) : 0;
}


static uint64_t markDuplicateFrags(
	std::vector< ::libmaus::bambam::ReadEnds > const & lfrags, DupSetCallback & DSC
)
{
	if ( lfrags.size() > 1 )
	{
		#if defined(MARKDUPLICATEFRAGSDEBUG)
		std::cerr << "[V] --- frag set --- " << std::endl;
		for ( uint64_t i = 0; i < lfrags.size(); ++i )
			std::cerr << "[V] " << lfrags[i] << std::endl;
		#endif
	
		bool containspairs = false;
		bool containsfrags = false;
		
		for ( uint64_t i = 0; i < lfrags.size(); ++i )
			if ( lfrags[i].isPaired() )
				containspairs = true;
			else
				containsfrags = true;
		
		// if there are any single fragments
		if ( containsfrags )
		{
			// mark single ends as duplicates if there are pairs
			if ( containspairs )
			{
				#if defined(MARKDUPLICATEFRAGSDEBUG)
				std::cerr << "[V] there are pairs, marking single ends as duplicates" << std::endl;
				#endif
			
				uint64_t dupcnt = 0;
				// std::cerr << "Contains pairs." << std::endl;
				for ( uint64_t i = 0; i < lfrags.size(); ++i )
					if ( !lfrags[i].isPaired() )
					{
						DSC(lfrags[i]);
						dupcnt++;
					}
					
				return dupcnt;	
			}
			// if all are single keep highest score only
			else
			{
				#if defined(MARKDUPLICATEFRAGSDEBUG)
				std::cerr << "[V] there are only fragments, keeping best one" << std::endl;
				#endif
				// std::cerr << "Frags only." << std::endl;
			
				uint64_t maxscore = lfrags[0].getScore();
				uint64_t maxindex = 0;
				
				for ( uint64_t i = 1; i < lfrags.size(); ++i )
					if ( lfrags[i].getScore() > maxscore )
					{
						maxscore = lfrags[i].getScore();
						maxindex = i;
					}

				for ( uint64_t i = 0; i < lfrags.size(); ++i )
					if ( i != maxindex )
						DSC(lfrags[i]);

				return  lfrags.size()-1;
			}			
		}
		else
		{
			#if defined(MARKDUPLICATEFRAGSDEBUG)
			std::cerr << "[V] group does not contain unpaired reads." << std::endl;
			#endif
		
			return 0;
		}
	}	
	else
	{
		return 0;
	}
}

static bool isDupPair(::libmaus::bambam::ReadEnds const & A, ::libmaus::bambam::ReadEnds const & B)
{
	bool const notdup = 
		A.getLibraryId()       != B.getLibraryId()       ||
		A.getRead1Sequence()   != B.getRead1Sequence()   ||
		A.getRead1Coordinate() != B.getRead1Coordinate() ||
		A.getOrientation()     != B.getOrientation()     ||
		A.getRead2Sequence()   != B.getRead2Sequence()   ||
		A.getRead2Coordinate() != B.getRead2Coordinate()
	;
	
	return ! notdup;
}

static bool isDupFrag(::libmaus::bambam::ReadEnds const & A, ::libmaus::bambam::ReadEnds const & B)
{
	bool const notdup = 
		A.getLibraryId()      != B.getLibraryId()       ||
		A.getRead1Sequence()   != B.getRead1Sequence()   ||
		A.getRead1Coordinate() != B.getRead1Coordinate() ||
		A.getOrientation()     != B.getOrientation()
	;
	
	return ! notdup;
}

struct AlignmentPairListNode
{
	typedef AlignmentPairListNode this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus::bambam::BamAlignment A[2];
	AlignmentPairListNode * next;
	
	AlignmentPairListNode() : next(0) {}
};

struct AlignmentListNode
{
	typedef AlignmentListNode this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus::bambam::BamAlignment A;
	AlignmentListNode * next;
	
	AlignmentListNode() : next(0) {}
};

template<typename _element_type>
struct FreeList
{
	typedef _element_type element_type;
	typedef FreeList<element_type> this_type;
	typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus::autoarray::AutoArray< element_type * > freelist;
	uint64_t freecnt;
	
	void cleanup()
	{
		for ( uint64_t i = 0; i < freelist.size(); ++i )
		{
			delete freelist[i];
			freelist[i] = 0;
		}	
	}
	
	FreeList(uint64_t const numel) : freelist(numel), freecnt(numel)
	{
		try
		{
			for ( uint64_t i = 0; i < numel; ++i )
				freelist[i] = 0;
			for ( uint64_t i = 0; i < numel; ++i )
				freelist[i] = new element_type;
		}
		catch(...)
		{
			cleanup();
			throw;
		}
	}
	
	~FreeList()
	{
		cleanup();
	}
	
	bool empty() const
	{
		return freecnt == 0;
	}
	
	element_type * get()
	{
		element_type * p = 0;
		assert ( ! empty() );

		p = freelist[--freecnt];
		freelist[freecnt] = 0;

		return p;
	}
	
	void put(element_type * ptr)
	{
		freelist[freecnt++] = ptr;
	}
};

typedef FreeList<AlignmentPairListNode> AlignmentPairFreeList;
typedef FreeList<AlignmentListNode> AlignmentFreeList;

struct PairActiveCount
{
	int32_t refid;
	int32_t coordinate;
	uint64_t incnt;
	uint64_t outcnt;	
	AlignmentPairListNode * root;
	bool expunge;
	
	PairActiveCount() : refid(-1), coordinate(-1), incnt(0), outcnt(0), root(0), expunge(false) {}
	PairActiveCount(
		int32_t const rrefid,
		int32_t const rcoordinate,
		uint64_t const rincnt,
		uint64_t const routcnt
	) : refid(rrefid), coordinate(rcoordinate), incnt(rincnt), outcnt(routcnt), root(0), expunge(false)
	{	
	}
	~PairActiveCount()
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
	
	bool operator<(PairActiveCount const & o) const
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
	
	void freeAlignments(AlignmentPairFreeList & list)
	{
		AlignmentPairListNode * p = root;
		
		while ( p )
		{
			AlignmentPairListNode * q = p->next;
			p->next = 0;
			list.put(p);
			p = q;
		}
		
		root = 0;
	}

	void addAlignmentPair(AlignmentPairListNode * ptr)
	{
		ptr->next = root;
		root = ptr;
	}
};

std::ostream & operator<<(std::ostream & out, PairActiveCount const & A)
{
	out << "PairActiveCount(" << A.refid << "," << A.coordinate << "," << A.incnt << "," << A.outcnt << ")";
	return out;
}

struct PositionTrackInterface
{
	static unsigned int const freelistsize = 16*1024;

	libmaus::bambam::BamHeader const & bamheader;
	std::pair<int32_t,int32_t> position;

	std::pair<int32_t,int32_t> expungepositionpairs;
	std::deque<PairActiveCount> activepairs;
	int64_t totalactivepairs;
	AlignmentPairFreeList APFLpairs;
	uint64_t excntpairs;
	uint64_t fincntpairs;
	uint64_t strcntpairs;
	std::vector<libmaus::bambam::ReadEnds> REpairs;
	
	std::pair<int32_t,int32_t> expungepositionfrags;

	DupSetCallback * DSC;

	PositionTrackInterface(libmaus::bambam::BamHeader const & rbamheader)
	: 
		bamheader(rbamheader), position(-1,-1), 
		expungepositionpairs(-1,-1), totalactivepairs(0), APFLpairs(freelistsize), excntpairs(0), fincntpairs(0), strcntpairs(0),
		expungepositionfrags(-1,-1),
		DSC(0)
	{
	}
	
	void setDupSetCallback(DupSetCallback * rDSC)
	{
		DSC = rDSC;
	}

	virtual ~PositionTrackInterface() {}

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
	
		position.first = refid;
		position.second = pos;
		
		if ( isSimplePair(A) && A.isReverse() )
		{
			int32_t const coord = A.getCoordinate();
			std::pair<int32_t,int32_t> pcoord(refid,coord);
			PairActiveCount acomp(refid,coord,0,0);
			
			// we have not seen the coordinate before
			if ( ! activepairs.size() || activepairs.back() < acomp )
			{
				activepairs.push_back(PairActiveCount(refid,coord,1,0));
			}
			// increment at end
			else if ( activepairs.back().refid == refid && activepairs.back().coordinate == coord )
			{
				activepairs.back().incIn();
			}
			else
			{
				std::deque<PairActiveCount>::iterator it = 
					std::lower_bound(activepairs.begin(),activepairs.end(),acomp);

				if ( it != activepairs.end() && *it == pcoord )
				{
					it->incIn();
				}
				else
				{				
					uint64_t const tomove = activepairs.end()-it;
					
					activepairs.push_back(PairActiveCount());
					
					for ( uint64_t j = 0; j < tomove; ++j )
						activepairs [ activepairs.size()-j-1 ] = activepairs[activepairs.size()-j-2];
						
					activepairs [ activepairs.size() - tomove - 1 ] = PairActiveCount(refid,coord,1,0);
				}
			}

			totalactivepairs++;
		}
	}

	void finishActiveFront()
	{	
		assert ( activepairs.size() );
	
		PairActiveCount & AC = activepairs.front();

		uint64_t lfincntpairs = 0;
		for ( AlignmentPairListNode * ptr = AC.root; ptr; ptr = ptr->next )
		{
			if ( lfincntpairs < REpairs.size() )
				libmaus::bambam::ReadEndsBase::fillFragPair(
					ptr->A[0],
					ptr->A[1],
					bamheader,
					REpairs[lfincntpairs]);
			else
				REpairs.push_back(libmaus::bambam::ReadEnds(ptr->A[0],ptr->A[1],bamheader));

			lfincntpairs += 1;
		}

		std::sort(REpairs.begin(),REpairs.begin()+lfincntpairs);
		
		fincntpairs += lfincntpairs;
		
		uint64_t l = 0;
		while ( l != lfincntpairs )
		{
			uint64_t h = l+1;
			while ( h != lfincntpairs && isDupPair(REpairs[l],REpairs[h]) )
				++h;
				
			if ( h-l > 1 )
			{
				markDuplicatePairs(REpairs.begin()+l,REpairs.begin()+h,*DSC);
			
				#if 0
				DupSetCallbackSet S0, S1;
				markDuplicatePairs(REpairs.begin()+l,REpairs.begin()+h,S0);
				std::vector < libmaus::bambam::ReadEnds > VV(REpairs.begin()+l,REpairs.begin()+h);
				markDuplicatePairsRef(VV,S1);
				
				assert ( S0.S == S1.S );
				#endif
			
				#if 0
				std::cerr << std::string(80,'-') << std::endl;
				for ( uint64_t i = l; i < h; ++i )
				{
					std::cerr << REpairs[i] << std::endl;
				}
				std::cerr << std::string(80,'-') << std::endl;
				#endif
			}
			
			l = h;
		}
		
		AC.freeAlignments(APFLpairs);
		totalactivepairs -= AC.incnt;
		assert ( totalactivepairs >= 0 );
		activepairs.pop_front();	
	}

	void expungeActiveFront(::libmaus::bambam::ReadEndsContainer * pairREC, ::libmaus::bambam::BamHeader const & header)
	{
		assert ( activepairs.size() );
	
		PairActiveCount & AC = activepairs.front();
		
		uint64_t lexcntpairs = 0;
		for ( AlignmentPairListNode * ptr = AC.root; ptr; ptr = ptr->next )
		{
			pairREC->putPair(ptr->A[0],ptr->A[1],header);
			lexcntpairs += 1;
		}
		
		excntpairs += lexcntpairs;
		
		expungepositionpairs.first = AC.refid;
		expungepositionpairs.second = AC.coordinate;
		AC.freeAlignments(APFLpairs);

		activepairs.pop_front();	
	}

	bool isActive(::libmaus::bambam::BamAlignment const & B)
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
	
	void setExpunge(libmaus::bambam::BamAlignment const & B)
	{
		if ( isActive(B) )
		{
			PairActiveCount const bkey(B.getRefID(),B.getCoordinate(),0,0);

			// key is not there or we can insert at the end, insert it
			if ( ! activepairs.size() || activepairs.back() < bkey )
			{
				activepairs.push_back(bkey);
				activepairs.back().setExpunge(true);
			}
			else 
			{
				std::deque<PairActiveCount>::iterator const it = std::lower_bound(activepairs.begin(),activepairs.end(),bkey);
			
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
						
					activepairs.push_back(PairActiveCount());
								
					for ( uint64_t j = 0; j < tomove; ++j )
					activepairs [ activepairs.size()-j-1 ] = activepairs[activepairs.size()-j-2];
								
					activepairs [ activepairs.size() - tomove - 1 ] = bkey;
					activepairs [ activepairs.size() - tomove - 1 ].setExpunge(true);
				}
			}
		}
	}
	
	void expungeUntil(
		libmaus::bambam::BamAlignment const & A,
		::libmaus::bambam::ReadEndsContainer * pairREC, 
		::libmaus::bambam::BamHeader const & header
	)
	{
		if ( isActive(A) )
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
				expungeActiveFront(pairREC,header);
			}

			expungepositionpairs.first = A.getRefID();
			expungepositionpairs.second = A.getCoordinate();
		}
	}
	
	/**
	 * flush lists
	 **/
	void flush(::libmaus::bambam::ReadEndsContainer * pairREC, ::libmaus::bambam::BamHeader const & header)
	{
		while ( activepairs.size() )
		{
			if ( activepairs.front().incnt == activepairs.front().outcnt )
			{
				if ( activepairs.front().getExpunge() )
				{
					expungeActiveFront(pairREC,header);			
				}
				else
				{
					finishActiveFront();
				}
			}
			else
			{
				std::cerr << "WARNING: expunge on flush (this should not happen)" << std::endl;
				expungeActiveFront(pairREC,header);
			}
		}
	}
	
	
	/**
	 * add a pair
	 **/
	void addAlignmentPair(
		::libmaus::bambam::BamAlignment const & A, ::libmaus::bambam::BamAlignment const & B,
		::libmaus::bambam::ReadEndsContainer * pairREC,
		::libmaus::bambam::BamHeader const & header
	)
	{
		PairActiveCount const bkey(B.getRefID(),B.getCoordinate(),0,0);
		
		bool done = false;
		
		while ( ! done )
		{
			bool const isactivepairs = isActive(B);
					
			if ( isactivepairs )
			{
				// expunge front element
				if ( APFLpairs.empty() )
				{
					assert ( activepairs.size() );
				
					// std::cerr << "Expunging " << activepairs.front() << " because free list is empty." << std::endl;	
					expungeActiveFront(pairREC,header);
					
					// check if this made any finished elements visible at the
					// front of the queue
					checkFinished(pairREC,bamheader);
				}
				//
				else
				{
					// find PairActiveCount object
					PairActiveCount const bkey(B.getRefID(),B.getCoordinate(),0,0);
					std::deque<PairActiveCount>::iterator const ita = std::lower_bound(activepairs.begin(),activepairs.end(),bkey);
					assert ( ita != activepairs.end() );
					assert ( ita->refid == bkey.refid && ita->coordinate == bkey.coordinate );
					
					// copy alignments
					AlignmentPairListNode * ptr = APFLpairs.get();
					ptr->A[0].copyFrom(A);
					ptr->A[1].copyFrom(B);
					ita->addAlignmentPair(ptr);
					ita->incOut();
										
					// done inserting this one
					done = true;
				}
			}
			else
			{
				// interval was already handled and this read pair
				// is too late, handle it by writing it out
				pairREC->putPair(A,B,header);
				excntpairs += 1;
				done = true;
			}
		}		
	}
	
	void checkFinished(
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
				expungeActiveFront(pairREC,header);			
			}
			else
			{
				finishActiveFront();
			}
			// finishActiveFront();
		}	
	}
};

struct PositionTrackCallback : 
	public ::libmaus::bambam::CollatingBamDecoderAlignmentInputCallback,
	public PositionTrackInterface
{
	PositionTrackCallback(libmaus::bambam::BamHeader const & bamheader) : PositionTrackInterface(bamheader) {}
	
	void operator()(::libmaus::bambam::BamAlignment const & A)
	{
		updatePosition(A);
	}
};

struct SnappyRewriteCallback : 
	public ::libmaus::bambam::CollatingBamDecoderAlignmentInputCallback,
	public PositionTrackInterface
{
	typedef SnappyRewriteCallback this_type;
	typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef ::libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t als;
	::libmaus::lz::SnappyFileOutputStream::unique_ptr_type SFOS;

	SnappyRewriteCallback(std::string const & filename, libmaus::bambam::BamHeader const & bamheader)
	: PositionTrackInterface(bamheader), als(0), SFOS(new ::libmaus::lz::SnappyFileOutputStream(filename))
	{
		
	}
	
	~SnappyRewriteCallback()
	{
		flush();
	}

	void operator()(::libmaus::bambam::BamAlignment const & A)
	{
		// std::cerr << "Callback." << std::endl;
		als++;
		updatePosition(A);
		A.serialise(*SFOS);
	}
	
	void flush()
	{
		SFOS->flush();
	}
};

struct BamRewriteCallback : 
	public ::libmaus::bambam::CollatingBamDecoderAlignmentInputCallback,
	public PositionTrackInterface
{
	typedef BamRewriteCallback this_type;
	typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef ::libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t als;
	::libmaus::bambam::BamWriter::unique_ptr_type BWR;
	
	BamRewriteCallback(
		std::string const & filename,
		::libmaus::bambam::BamHeader const & bamheader,
		int const rewritebamlevel
	)
	: PositionTrackInterface(bamheader), als(0), BWR(new ::libmaus::bambam::BamWriter(filename,bamheader,rewritebamlevel))
	{
		
	}
	
	~BamRewriteCallback()
	{
	}

	void operator()(::libmaus::bambam::BamAlignment const & A)
	{
		// std::cerr << "Callback." << std::endl;
		als++;
		updatePosition(A);
		A.serialise(BWR->getStream());
	}	
};

struct SnappyRewrittenInput
{
	::libmaus::lz::SnappyFileInputStream GZ;
	::libmaus::bambam::BamAlignment alignment;
	
	SnappyRewrittenInput(std::string const & filename)
	: GZ(filename)
	{
	
	}
	
	::libmaus::bambam::BamAlignment & getAlignment()
	{
		return alignment;
	}

	::libmaus::bambam::BamAlignment const & getAlignment() const
	{
		return alignment;
	}
	
	bool readAlignment()
	{
		/* read alignment block size */
		int64_t const bs0 = GZ.get();
		int64_t const bs1 = GZ.get();
		int64_t const bs2 = GZ.get();
		int64_t const bs3 = GZ.get();
		if ( bs3 < 0 )
			// reached end of file
			return false;
		
		/* assemble block size as LE integer */
		alignment.blocksize = (bs0 << 0) | (bs1 << 8) | (bs2 << 16) | (bs3 << 24) ;

		/* read alignment block */
		if ( alignment.blocksize > alignment.D.size() )
			alignment.D = ::libmaus::bambam::BamAlignment::D_array_type(alignment.blocksize);
		GZ.read(reinterpret_cast<char *>(alignment.D.begin()),alignment.blocksize);
		// assert ( static_cast<int64_t>(GZ.gcount()) == static_cast<int64_t>(alignment.blocksize) );
	
		return true;
	}
};

::libmaus::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus::util::ArgInfo const & arginfo,
	::libmaus::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bammarkduplicates", // ID
		"bammarkduplicates", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus::bambam::BamHeader(upheadtext));
	
	return UNIQUE_PTR_MOVE(uphead);
}

void addBamDuplicateFlag(
	::libmaus::util::ArgInfo const & arginfo,
	bool const verbose,
	::libmaus::bambam::BamHeader const & bamheader,
	uint64_t const maxrank,
	uint64_t const mod,
	int const level,
	DupSetCallback const & DSC,
	std::istream & in
)
{
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,bamheader));

	::libmaus::aio::CheckedOutputStream::unique_ptr_type pO;
	std::ostream * poutputstr = 0;
	
	if ( arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != "") )
	{
		::libmaus::aio::CheckedOutputStream::unique_ptr_type tpO(
                                new ::libmaus::aio::CheckedOutputStream(arginfo.getValue<std::string>("O",std::string("O")))
                        );
		pO = UNIQUE_PTR_MOVE(tpO);
		poutputstr = pO.get();
	}
	else
	{
		poutputstr = & std::cout;
	}

	std::ostream & outputstr = *poutputstr;

	/* write bam header */
	{
		libmaus::lz::BgzfDeflate<std::ostream> headout(outputstr);
		uphead->serialise(headout);
		headout.flush();
	}

	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;

	libmaus::timing::RealTimeClock globrtc; globrtc.start();
	libmaus::timing::RealTimeClock locrtc; locrtc.start();
	libmaus::lz::BgzfRecode rec(in,outputstr,level);
	
	bool haveheader = false;
	uint64_t blockskip = 0;
	std::vector<uint8_t> headerstr;
	uint64_t preblocksizes = 0;
	
	/* read and copy blocks until we have reached the end of the BAM header */
	while ( (!haveheader) && rec.getBlock() )
	{
		std::copy ( rec.deflatebase.pa, rec.deflatebase.pa + rec.P.second,
			std::back_insert_iterator < std::vector<uint8_t> > (headerstr) );
		
		try
		{
			libmaus::util::ContainerGetObject< std::vector<uint8_t> > CGO(headerstr);
			libmaus::bambam::BamHeader header;
			header.init(CGO);
			haveheader = true;
			blockskip = CGO.i - preblocksizes;
		}
		catch(std::exception const & ex)
		{
			#if 0
			std::cerr << "[D] " << ex.what() << std::endl;
			#endif
		}
	
		// need to read another block to get header, remember size of current block
		if ( ! haveheader )
			preblocksizes += rec.P.second;
	}
	
	if ( blockskip )
	{
		uint64_t const bytesused = (rec.deflatebase.pc - rec.deflatebase.pa) - blockskip;
		memmove ( rec.deflatebase.pa, rec.deflatebase.pa + blockskip, bytesused );
		rec.deflatebase.pc = rec.deflatebase.pa + bytesused;
		rec.P.second = bytesused;
		blockskip = 0;
		
		if ( ! bytesused )
			rec.getBlock();
	}

	/* parser state types and variables */
	enum parsestate { state_reading_blocklen, state_pre_skip, state_marking, state_post_skip };
	parsestate state = state_reading_blocklen;
	unsigned int blocklenred = 0;
	uint32_t blocklen = 0;
	uint32_t preskip = 0;
	uint64_t alcnt = 0;
	unsigned int const dupflagskip = 15;
	// uint8_t const dupflagmask = static_cast<uint8_t>(~(4u));
			
	/* while we have alignment data blocks */
	while ( rec.P.second )
	{
		uint8_t * pa       = rec.deflatebase.pa;
		uint8_t * const pc = rec.deflatebase.pc;

		while ( pa != pc )
			switch ( state )
			{
				/* read length of next alignment block */
				case state_reading_blocklen:
					/* if this is a little endian machine allowing unaligned access */
					#if defined(LIBMAUS_HAVE_i386)
					if ( (!blocklenred) && ((pc-pa) >= static_cast<ptrdiff_t>(sizeof(uint32_t))) )
					{
						blocklen = *(reinterpret_cast<uint32_t const *>(pa));
						blocklenred = sizeof(uint32_t);
						pa += sizeof(uint32_t);
						
						state = state_pre_skip;
						preskip = dupflagskip;
					}
					else
					#endif
					{
						while ( pa != pc && blocklenred < sizeof(uint32_t) )
							blocklen |= static_cast<uint32_t>(*(pa++)) << ((blocklenred++)*8);

						if ( blocklenred == sizeof(uint32_t) )
						{
							state = state_pre_skip;
							preskip = dupflagskip;
						}
					}
					break;
				/* skip data before the part we modify */
				case state_pre_skip:
					{
						uint32_t const skip = std::min(pc-pa,static_cast<ptrdiff_t>(preskip));
						pa += skip;
						preskip -= skip;
						blocklen -= skip;
						
						if ( ! skip )
							state = state_marking;
					}
					break;
				/* change data */
				case state_marking:
					assert ( pa != pc );
					if ( DSC.isMarked(alcnt) )
						*pa |= 4;
					state = state_post_skip;
					// intended fall through to post_skip case
				/* skip data after part we modify */
				case state_post_skip:
				{
					uint32_t const skip = std::min(pc-pa,static_cast<ptrdiff_t>(blocklen));
					pa += skip;
					blocklen -= skip;

					if ( ! blocklen )
					{
						state = state_reading_blocklen;
						blocklenred = 0;
						blocklen = 0;
						alcnt++;
						
						if ( verbose && ((alcnt & (bmask)) == 0) )
						{
							std::cerr << "[V] Marked " << (alcnt+1) << " (" << (alcnt+1)/(1024*1024) << "," << static_cast<double>(alcnt+1)/maxrank << ")"
								<< " time " << locrtc.getElapsedSeconds()
								<< " total " << globrtc.formatTime(globrtc.getElapsedSeconds())
								<< " "
								<< libmaus::util::MemUsage()
								<< std::endl;
							locrtc.start();
						}
					}
					break;
				}
			}

		rec.putBlock();			
		rec.getBlock();
	}
			
	rec.addEOFBlock();
	outputstr.flush();
	
	if ( verbose )
		std::cerr << "[V] Marked " << 1.0 << " total for marking time "
			<< globrtc.formatTime(globrtc.getElapsedSeconds()) 
			<< " "
			<< libmaus::util::MemUsage()
			<< std::endl;
}

namespace libmaus
{
	namespace lz
	{
		struct BgzfParallelRecodeDeflateBase : public ::libmaus::lz::BgzfConstants
		{
			libmaus::autoarray::AutoArray<uint8_t> B;
			
			uint8_t * const pa;
			uint8_t * pc;
			uint8_t * const pe;
			
			BgzfParallelRecodeDeflateBase()
			: B(getBgzfMaxBlockSize(),false), 
			  pa(B.begin()), 
			  pc(B.begin()), 
			  pe(B.end())
			{
			
			}
		};
	}
}

namespace libmaus
{
	namespace lz
	{
		struct BgzfRecodeParallel : public ::libmaus::lz::BgzfConstants
		{
			libmaus::lz::BgzfInflateDeflateParallel BIDP; // (std::cin,std::cout,Z_DEFAULT_COMPRESSION,32,128);
			libmaus::lz::BgzfParallelRecodeDeflateBase deflatebase;
			std::pair<uint64_t,uint64_t> P;

			BgzfRecodeParallel(
				std::istream & in, std::ostream & out,
				int const level, // = Z_DEFAULT_COMPRESSION,
				uint64_t const numthreads,
				uint64_t const numbuffers
			) : BIDP(in,out,level,numthreads,numbuffers), deflatebase(), P(0,0)
			{
			
			}
			
			~BgzfRecodeParallel()
			{
				BIDP.flush();
			}
			
			bool getBlock()
			{
				P.second = BIDP.read(reinterpret_cast<char *>(deflatebase.B.begin()),deflatebase.B.size());
				deflatebase.pc = deflatebase.pa + P.second;
				
				return P.second != 0;
			}
			
			void putBlock()
			{
				BIDP.write(reinterpret_cast<char const *>(deflatebase.B.begin()),P.second);
			}
			
			void addEOFBlock()
			{
				BIDP.flush();
			}
		};
	}
}

void addBamDuplicateFlagParallel(
	::libmaus::util::ArgInfo const & arginfo,
	bool const verbose,
	::libmaus::bambam::BamHeader const & bamheader,
	uint64_t const maxrank,
	uint64_t const mod,
	int const level,
	DupSetCallback const & DSC,
	std::istream & in,
	uint64_t const numthreads
)
{
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,bamheader));

	::libmaus::aio::CheckedOutputStream::unique_ptr_type pO;
	std::ostream * poutputstr = 0;
	
	if ( arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != "") )
	{
		::libmaus::aio::CheckedOutputStream::unique_ptr_type tpO(
                                new ::libmaus::aio::CheckedOutputStream(arginfo.getValue<std::string>("O",std::string("O")))
                        );
		pO = UNIQUE_PTR_MOVE(tpO);
		poutputstr = pO.get();
	}
	else
	{
		poutputstr = & std::cout;
	}

	std::ostream & outputstr = *poutputstr;

	/* write bam header */
	{
		libmaus::lz::BgzfDeflate<std::ostream> headout(outputstr);
		uphead->serialise(headout);
		headout.flush();
	}

	libmaus::lz::BgzfRecodeParallel rec(in,outputstr,level,numthreads,numthreads*4);

	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;

	libmaus::timing::RealTimeClock globrtc; globrtc.start();
	libmaus::timing::RealTimeClock locrtc; locrtc.start();
	
	bool haveheader = false;
	uint64_t blockskip = 0;
	std::vector<uint8_t> headerstr;
	uint64_t preblocksizes = 0;
	
	/* read and copy blocks until we have reached the end of the BAM header */
	while ( (!haveheader) && rec.getBlock() )
	{
		std::copy ( rec.deflatebase.pa, rec.deflatebase.pa + rec.P.second,
			std::back_insert_iterator < std::vector<uint8_t> > (headerstr) );
		
		try
		{
			libmaus::util::ContainerGetObject< std::vector<uint8_t> > CGO(headerstr);
			libmaus::bambam::BamHeader header;
			header.init(CGO);
			haveheader = true;
			blockskip = CGO.i - preblocksizes;
		}
		catch(std::exception const & ex)
		{
			#if 0
			std::cerr << "[D] " << ex.what() << std::endl;
			#endif
		}
	
		// need to read another block to get header, remember size of current block
		if ( ! haveheader )
			preblocksizes += rec.P.second;
	}
	
	if ( blockskip )
	{
		uint64_t const bytesused = (rec.deflatebase.pc - rec.deflatebase.pa) - blockskip;
		memmove ( rec.deflatebase.pa, rec.deflatebase.pa + blockskip, bytesused );
		rec.deflatebase.pc = rec.deflatebase.pa + bytesused;
		rec.P.second = bytesused;
		blockskip = 0;
		
		if ( ! bytesused )
			rec.getBlock();
	}

	/* parser state types and variables */
	enum parsestate { state_reading_blocklen, state_pre_skip, state_marking, state_post_skip };
	parsestate state = state_reading_blocklen;
	unsigned int blocklenred = 0;
	uint32_t blocklen = 0;
	uint32_t preskip = 0;
	uint64_t alcnt = 0;
	unsigned int const dupflagskip = 15;
	// uint8_t const dupflagmask = static_cast<uint8_t>(~(4u));
			
	/* while we have alignment data blocks */
	while ( rec.P.second )
	{
		uint8_t * pa       = rec.deflatebase.pa;
		uint8_t * const pc = rec.deflatebase.pc;

		while ( pa != pc )
			switch ( state )
			{
				/* read length of next alignment block */
				case state_reading_blocklen:
					/* if this is a little endian machine allowing unaligned access */
					#if defined(LIBMAUS_HAVE_i386)
					if ( (!blocklenred) && ((pc-pa) >= static_cast<ptrdiff_t>(sizeof(uint32_t))) )
					{
						blocklen = *(reinterpret_cast<uint32_t const *>(pa));
						blocklenred = sizeof(uint32_t);
						pa += sizeof(uint32_t);
						
						state = state_pre_skip;
						preskip = dupflagskip;
					}
					else
					#endif
					{
						while ( pa != pc && blocklenred < sizeof(uint32_t) )
							blocklen |= static_cast<uint32_t>(*(pa++)) << ((blocklenred++)*8);

						if ( blocklenred == sizeof(uint32_t) )
						{
							state = state_pre_skip;
							preskip = dupflagskip;
						}
					}
					break;
				/* skip data before the part we modify */
				case state_pre_skip:
					{
						uint32_t const skip = std::min(pc-pa,static_cast<ptrdiff_t>(preskip));
						pa += skip;
						preskip -= skip;
						blocklen -= skip;
						
						if ( ! skip )
							state = state_marking;
					}
					break;
				/* change data */
				case state_marking:
					assert ( pa != pc );
					if ( DSC.isMarked(alcnt) )
						*pa |= 4;
					state = state_post_skip;
					// intended fall through to post_skip case
				/* skip data after part we modify */
				case state_post_skip:
				{
					uint32_t const skip = std::min(pc-pa,static_cast<ptrdiff_t>(blocklen));
					pa += skip;
					blocklen -= skip;

					if ( ! blocklen )
					{
						state = state_reading_blocklen;
						blocklenred = 0;
						blocklen = 0;
						alcnt++;
						
						if ( verbose && ((alcnt & (bmask)) == 0) )
						{
							std::cerr << "[V] Marked " << (alcnt+1) << " (" << (alcnt+1)/(1024*1024) << "," << static_cast<double>(alcnt+1)/maxrank << ")"
								<< " time " << locrtc.getElapsedSeconds()
								<< " total " << globrtc.formatTime(globrtc.getElapsedSeconds())
								<< " "
								<< libmaus::util::MemUsage()
								<< std::endl;
							locrtc.start();
						}
					}
					break;
				}
			}

		rec.putBlock();			
		rec.getBlock();
	}
			
	rec.addEOFBlock();
	outputstr.flush();
	
	if ( verbose )
		std::cerr << "[V] Marked " << 1.0 << " total for marking time "
			<< globrtc.formatTime(globrtc.getElapsedSeconds()) 
			<< " "
			<< libmaus::util::MemUsage()
			<< std::endl;
}


template<typename decoder_type>
static void markDuplicatesInFileTemplate(
	::libmaus::util::ArgInfo const & arginfo,
	bool const verbose,
	::libmaus::bambam::BamHeader const & bamheader,
	uint64_t const maxrank,
	uint64_t const mod,
	int const level,
	DupSetCallback const & DSC,
	decoder_type & decoder
)
{
	::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,bamheader));

	::libmaus::aio::CheckedOutputStream::unique_ptr_type pO;
	std::ostream * poutputstr = 0;
	
	if ( arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != "") )
	{
		::libmaus::aio::CheckedOutputStream::unique_ptr_type tpO(
                                new ::libmaus::aio::CheckedOutputStream(arginfo.getValue<std::string>("O",std::string("O")))
                        );
		pO = UNIQUE_PTR_MOVE(tpO);
		poutputstr = pO.get();
	}
	else
	{
		poutputstr = & std::cout;
	}

	std::ostream & outputstr = *poutputstr;
	
	libmaus::timing::RealTimeClock globrtc, locrtc;
	globrtc.start(); locrtc.start();
	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;

	// rewrite file and mark duplicates
	::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(outputstr,*uphead,level));
	libmaus::bambam::BamAlignment & alignment = decoder.getAlignment();
	for ( uint64_t r = 0; decoder.readAlignment(); ++r )
	{
		if ( DSC.isMarked(r) )
			alignment.putFlags(alignment.getFlags() | ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP);
		
		alignment.serialise(writer->getStream());
		
		if ( verbose && ((r+1) & bmask) == 0 )
		{
			std::cerr << "[V] Marked " << (r+1) << " (" << (r+1)/(1024*1024) << "," << static_cast<double>(r+1)/maxrank << ")"
				<< " time " << locrtc.getElapsedSeconds()
				<< " total " << globrtc.formatTime(globrtc.getElapsedSeconds())
				<< " "
				<< libmaus::util::MemUsage()
				<< std::endl;
			locrtc.start();
		}
	}
	
	writer.reset();
	outputstr.flush();
	pO.reset();
	
	if ( verbose )
		std::cerr << "[V] Marked " << maxrank << "(" << maxrank/(1024*1024) << "," << 1 << ")" << " total for marking time " << globrtc.formatTime(globrtc.getElapsedSeconds()) 
			<< " "
			<< libmaus::util::MemUsage()
			<< std::endl;

}


template<typename decoder_type, typename writer_type>
static void removeDuplicatesFromFileTemplate(
	bool const verbose,
	uint64_t const maxrank,
	uint64_t const mod,
	DupSetCallback const & DSC,
	decoder_type & decoder,
	writer_type & writer
)
{
	libmaus::timing::RealTimeClock globrtc, locrtc;
	globrtc.start(); locrtc.start();
	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;

	// rewrite file and mark duplicates
	libmaus::bambam::BamAlignment & alignment = decoder.getAlignment();
	for ( uint64_t r = 0; decoder.readAlignment(); ++r )
	{
		if ( ! DSC.isMarked(r) )
			alignment.serialise(writer.getStream());
		
		if ( verbose && ((r+1) & bmask) == 0 )
		{
			std::cerr << "[V] Filtered " << (r+1) << " (" << (r+1)/(1024*1024) << "," << static_cast<double>(r+1)/maxrank << ")"
				<< " time " << locrtc.getElapsedSeconds()
				<< " total " << globrtc.formatTime(globrtc.getElapsedSeconds())
				<< " "
				<< libmaus::util::MemUsage()
				<< std::endl;
			locrtc.start();
		}
	}
		
	if ( verbose )
		std::cerr << "[V] Filtered " << maxrank << "(" << maxrank/(1024*1024) << "," << 1 << ")" << " total for marking time " << globrtc.formatTime(globrtc.getElapsedSeconds()) 
			<< " "
			<< libmaus::util::MemUsage()
			<< std::endl;

}

struct UpdateHeader : public libmaus::bambam::BamHeaderRewriteCallback
{
	libmaus::util::ArgInfo const & arginfo;

	UpdateHeader(libmaus::util::ArgInfo const & rarginfo)
	: arginfo(rarginfo)
	{
	
	}

	::libmaus::bambam::BamHeader::unique_ptr_type operator()(::libmaus::bambam::BamHeader const & header)  const
	{
		::libmaus::bambam::BamHeader::unique_ptr_type ptr(updateHeader(arginfo,header));
		return UNIQUE_PTR_MOVE(ptr);
	}
};

static void markDuplicatesInFile(
	::libmaus::util::ArgInfo const & arginfo,
	bool const verbose,
	::libmaus::bambam::BamHeader const & bamheader,
	uint64_t const maxrank,
	uint64_t const mod,
	int const level,
	DupSetCallback const & DSC,
	std::string const & recompressedalignments,
	bool const rewritebam
)
{
	bool const rmdup = arginfo.getValue<int>("rmdup",getDefaultRmDup());
	uint64_t const markthreads = 
		std::max(static_cast<uint64_t>(1),arginfo.getValue<uint64_t>("markthreads",getDefaultMarkThreads()));

	if ( rmdup )
	{
		::libmaus::bambam::BamHeader::unique_ptr_type uphead(updateHeader(arginfo,bamheader));
		
		bool const inputisbam =
			(arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != ""))
			||
			rewritebam;
	
		::libmaus::aio::CheckedOutputStream::unique_ptr_type pO;
		std::ostream * poutputstr = 0;
		
		if ( arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != "") )
		{
			::libmaus::aio::CheckedOutputStream::unique_ptr_type tpO(
                                        new ::libmaus::aio::CheckedOutputStream(arginfo.getValue<std::string>("O",std::string("O")))
                                );
			pO = UNIQUE_PTR_MOVE(tpO);
			poutputstr = pO.get();
		}
		else
		{
			poutputstr = & std::cout;
		}

		std::ostream & outputstr = *poutputstr;
		
		if ( inputisbam )
		{
			std::string const inputfilename =
				(arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != ""))
				?
				arginfo.getValue<std::string>("I","I")
				:
				recompressedalignments;

			if ( markthreads < 2 )
			{
				::libmaus::bambam::BamDecoder decoder(inputfilename);
				decoder.disableValidation();
				::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(outputstr,*uphead,level));	
				removeDuplicatesFromFileTemplate(verbose,maxrank,mod,DSC,decoder,*writer);
			}
			else
			{
				libmaus::aio::CheckedInputStream CIS(inputfilename);
				UpdateHeader UH(arginfo);
				libmaus::bambam::BamParallelRewrite BPR(CIS,UH,outputstr,level,markthreads,4 /* blocks per thread */);
				libmaus::bambam::BamAlignmentDecoder & dec = BPR.getDecoder();
				libmaus::bambam::BamParallelRewrite::writer_type & writer = BPR.getWriter();
				removeDuplicatesFromFileTemplate(verbose,maxrank,mod,DSC,dec,writer);
			}
		}
		else
		{		
			SnappyRewrittenInput decoder(recompressedalignments);
			if ( verbose )
				std::cerr << "[V] Reading snappy alignments from " << recompressedalignments << std::endl;
			::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(outputstr,*uphead,level));
			removeDuplicatesFromFileTemplate(verbose,maxrank,mod,DSC,decoder,*writer);
		}

		outputstr.flush();
		pO.reset();
	}
	else
	{
		if ( arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != "") )
		{
			std::string const inputfilename = arginfo.getValue<std::string>("I","I");
			libmaus::aio::CheckedInputStream CIS(inputfilename);
		
			if ( markthreads == 1 )
				addBamDuplicateFlag(arginfo,verbose,bamheader,maxrank,mod,level,DSC,CIS);
			else
				addBamDuplicateFlagParallel(arginfo,verbose,bamheader,maxrank,mod,level,DSC,CIS,markthreads);
		}
		else
		{
			if ( rewritebam )
			{
				libmaus::aio::CheckedInputStream CIS(recompressedalignments);
				
				if ( markthreads == 1 )
					addBamDuplicateFlag(arginfo,verbose,bamheader,maxrank,mod,level,DSC,CIS);		
				else
					addBamDuplicateFlagParallel(arginfo,verbose,bamheader,maxrank,mod,level,DSC,CIS,markthreads);
			}
			else
			{
				SnappyRewrittenInput decoder(recompressedalignments);
				if ( verbose )
					std::cerr << "[V] Reading snappy alignments from " << recompressedalignments << std::endl;
				markDuplicatesInFileTemplate(arginfo,verbose,bamheader,maxrank,mod,level,DSC,decoder);
			}
		}
	}
}

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
	int const rewritebamlevel = arginfo.getValue<int>("rewritebamlevel",getDefaultRewriteBamLevel());

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

	std::string const tmpfiledupset = tmpfilenamebase + "_dupset";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfiledupset);

	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	
	switch ( level )
	{
		case Z_NO_COMPRESSION:
		case Z_BEST_SPEED:
		case Z_BEST_COMPRESSION:
		case Z_DEFAULT_COMPRESSION:
			break;
		default:
		{
			::libmaus::exception::LibMausException se;
			se.getStream()
				<< "Unknown compression level, please use"
				<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
				<< " level=" << Z_BEST_SPEED << " (fast) or"
				<< " level=" << Z_BEST_COMPRESSION << " (best) or"
				<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
			se.finish();
			throw se;
		}
			break;
	}
	switch ( rewritebamlevel )
	{
		case Z_NO_COMPRESSION:
		case Z_BEST_SPEED:
		case Z_BEST_COMPRESSION:
		case Z_DEFAULT_COMPRESSION:
			break;
		default:
		{
			::libmaus::exception::LibMausException se;
			se.getStream()
				<< "Unknown value for rewritebamlevel, please use"
				<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
				<< " level=" << Z_BEST_SPEED << " (fast) or"
				<< " level=" << Z_BEST_COMPRESSION << " (best) or"
				<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
			se.finish();
			throw se;
		}
			break;
	}
	
	if ( verbose )
		std::cerr << "[V] output compression level " << level << std::endl;

	::libmaus::timing::RealTimeClock fragrtc; fragrtc.start();

	SnappyRewriteCallback::unique_ptr_type SRC;
	BamRewriteCallback::unique_ptr_type BWR;
	PositionTrackInterface * PTI = 0;
	::libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	libmaus::aio::CheckedOutputStream::unique_ptr_type copybamstr;

	typedef ::libmaus::bambam::BamCircularHashCollatingBamDecoder col_type;
	typedef ::libmaus::bambam::BamParallelCircularHashCollatingBamDecoder par_col_type;
	typedef ::libmaus::bambam::CircularHashCollatingBamDecoder col_base_type;
	typedef col_base_type::unique_ptr_type col_base_ptr_type;	
	col_base_ptr_type CBD;

	uint64_t const markthreads = 
		std::max(static_cast<uint64_t>(1),arginfo.getValue<uint64_t>("markthreads",getDefaultMarkThreads()));
	
	// if we are reading the input from a file
	if ( arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != "") )
	{
		std::string const inputfilename = arginfo.getValue<std::string>("I","I");
		::libmaus::aio::CheckedInputStream::unique_ptr_type tCIS(new ::libmaus::aio::CheckedInputStream(inputfilename));
		CIS = UNIQUE_PTR_MOVE(tCIS);
		
		if ( markthreads > 1 )
		{
			col_base_ptr_type tCBD(new par_col_type(
                                *CIS,
                                markthreads,
                                tmpfilename,
                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL*/ 0,
                                true /* put rank */,
                                colhashbits,collistsize));
			CBD = UNIQUE_PTR_MOVE(tCBD);
		}
		else
		{
			col_base_ptr_type tCBD(new col_type(
                                *CIS,
                                // numthreads
                                tmpfilename,
                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL*/ 0,
                                true /* put rank */,
                                colhashbits,collistsize));	
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
                                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL */ 0,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}
				else
				{
					col_base_ptr_type tCBD(new col_type(std::cin,
                                                *copybamstr,
                                                tmpfilename,
                                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL */ 0,
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
                                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL */ 0,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}
				else
				{
					col_base_ptr_type tCBD(new col_type(std::cin,
                                                tmpfilename,
                                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL */ 0,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}

				// rewrite file and mark duplicates
				BamRewriteCallback::unique_ptr_type tBWR(new BamRewriteCallback(tmpfilesnappyreads,CBD->getHeader(),rewritebamlevel));
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
                                /* libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY      | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL */ 0,
                                true /* put rank */,
                                colhashbits,collistsize));
			CBD = UNIQUE_PTR_MOVE(tCBD);

			SnappyRewriteCallback::unique_ptr_type tSRC(new SnappyRewriteCallback(tmpfilesnappyreads,CBD->getHeader()));
			SRC = UNIQUE_PTR_MOVE(tSRC);
			CBD->setInputCallback(SRC.get());
			if ( verbose )
				std::cerr << "[V] Writing snappy compressed alignments to file " << tmpfilesnappyreads << std::endl;
		}
	}

	::libmaus::bambam::BamHeader const bamheader = CBD->getHeader();

	PositionTrackCallback PTC(bamheader);
	
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

	DupSetCallbackStream DSCV(tmpfiledupset,metrics);
	// DupSetCallbackSet DSCV(metrics);
	
	PTI->setDupSetCallback(&DSCV);
	
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
			if (  PositionTrackInterface::isSimplePair(*(P.second)) )
			{
				PTI->addAlignmentPair(*(P.first),*(P.second),pairREC.get(),bamheader);
				PTI->checkFinished(pairREC.get(),bamheader);
			}
			// strangely mapped pairs, mark right coordinate as expunged
			else
			{
				PTI->setExpunge(*(P.second));
				pairREC->putPair(*(P.first),*(P.second),bamheader);
				PTI->strcntpairs += 1;
			}

			paircnt++;
		}
	
		if ( P.first )
		{
			fragREC->putFrag(*(P.first),bamheader);				
			fragcnt++;
		}
		if ( P.second )
		{
			fragREC->putFrag(*(P.second),bamheader);
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

	PTI->flush(pairREC.get(),bamheader);
	std::cerr << "excntpairs=" << PTI->excntpairs << " fincntpairs=" << PTI->fincntpairs << " strcntpairs=" << PTI->strcntpairs << std::endl;
	
	assert ( PTI->excntpairs + PTI->fincntpairs + PTI->strcntpairs == paircnt );
	
	if ( copybamstr )
	{
		copybamstr->flush();
		copybamstr.reset();
	}
	
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

	uint64_t const numranks = maxrank+1;
	
	if ( numranks != als )
		std::cerr << "[D] numranks=" << numranks << " != als=" << als << std::endl;
	
	assert ( numranks == als );

	if ( verbose )
		std::cerr
			<< "[V] " 
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
		if ( ! isDupPair(nextfrag,lfrags.front()) )
		{
			// dupcnt += markDuplicatePairs(lfrags.begin(),lfrags.end(),DSCV);
			dupcnt += markDuplicatePairsRef(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	// dupcnt += markDuplicatePairs(lfrags.begin(),lfrags.end(),DSCV);
	dupcnt += markDuplicatePairsRef(lfrags,DSCV);
	lfrags.resize(0);
	pairDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << (PTI->excntpairs+PTI->strcntpairs)/rtc.getElapsedSeconds() << std::endl;

	if ( verbose )
		std::cerr << "[V] Checking single fragments...";
	rtc.start();
	::libmaus::bambam::SortedFragDecoder::unique_ptr_type fragDec(fragREC->getDecoder());
	fragREC.reset();
	fragDec->getNext(lfrags);
	while ( fragDec->getNext(nextfrag) )
	{
		if ( !isDupFrag(nextfrag,lfrags.front()) )
		{
			dupcnt += markDuplicateFrags(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	dupcnt += markDuplicateFrags(lfrags,DSCV);
	lfrags.resize(0);
	fragDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << fragcnt/rtc.getElapsedSeconds() << std::endl;		

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
	markDuplicatesInFile(arginfo,verbose,bamheader,maxrank,mod,level,DSCV,tmpfilesnappyreads,rewritebam);
		
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
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "markthreads=<["+::biobambam::Licensing::formatNumber(getDefaultMarkThreads())+"]>", "number of helper threads" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report (default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for each mod'th record/alignment" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebam=<["+::biobambam::Licensing::formatNumber(getDefaultRewriteBam())+"]>", "compression of temporary alignment file when input is via stdin (0=snappy,1=gzip/bam,2=copy)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebamlevel=<["+::biobambam::Licensing::formatNumber(getDefaultRewriteBamLevel())+"]>", "compression settings for temporary alignment file if rewritebam=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rmdup=<["+::biobambam::Licensing::formatNumber(getDefaultRmDup())+"]>", "remove duplicates (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhashbits=<["+::biobambam::Licensing::formatNumber(getDefaultColHashBits())+"]>", "log_2 of size of hash table used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collistsize=<["+::biobambam::Licensing::formatNumber(getDefaultColListSize())+"]>", "output list size for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "fragbufsize=<["+::biobambam::Licensing::formatNumber(getDefaultFragBufSize())+"]>", "size of each fragment/pair file buffer in bytes" ) );

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

