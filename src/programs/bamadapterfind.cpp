/**
    bambam
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013, 2016 Genome Research Limited

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

#include <libmaus2/bambam/AdapterFilter.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/rank/popcnt.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/Histogram.hpp>

#include <libmaus2/parallel/SimpleThreadPool.hpp>
#include <libmaus2/parallel/SimpleThreadPoolWorkPackageFreeList.hpp>

#include <biobambam2/Licensing.hpp>
#include <biobambam2/ClipAdapters.hpp>
#include <biobambam2/KmerPoisson.hpp>
#include <biobambam2/BamBamConfig.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1024*1024; }
static uint64_t getDefaultPCT_MISMATCH() { return 10; }
static unsigned int getDefaultSEED_LENGTH() { return 12; }
static uint64_t getDefaultMIN_OVERLAP() { return 32; }
static uint64_t getDefaultADAPTER_MATCH() { return 12; }
static uint64_t getDefaultMatchMinScore() { return 16; }
static double getDefaultMatchMinFrac() { return 0.75; }
static double getDefaultMatchMinPFrac() { return 0.8; }
static int getDefaultClip() { return 0; }
static uint64_t getDefaultRefLen() { return 3000000000ull; }
static double getDefaultpA() { return 0.25; }
static double getDefaultpC() { return 0.25; }
static double getDefaultpG() { return 0.25; }
static double getDefaultpT() { return 0.25; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

static int getDefaultThreads() { return 8; }
static int getDefaultWriteQueue() { return 100; }
static int getDefaultWorkReads() { return 1000; }


/*
    Next few structures hold data that is passed
    between threads
*/
struct HoldAlign
{
    	libmaus2::bambam::BamAlignment algns[2];
	bool pair;
};


struct WorkBlock
{
	::libmaus2::util::ArgInfo const & arginfo;
	libmaus2::bambam::AdapterFilter::shared_ptr_type AF;
	libmaus2::autoarray::AutoArray<uint8_t> &S;
	libmaus2::autoarray::AutoArray<uint8_t> &R;
	uint64_t const & mmask;
	std::vector<HoldAlign *> entries;
	int wb_no;
	bool ending;
	libmaus2::util::Histogram overlaphist;
	libmaus2::util::Histogram adapterhist;
	uint64_t adptcnt;
	uint64_t paircnt;

	WorkBlock(::libmaus2::util::ArgInfo const & inarg, libmaus2::bambam::AdapterFilter::shared_ptr_type rAF, libmaus2::autoarray::AutoArray<uint8_t> &inS, libmaus2::autoarray::AutoArray<uint8_t> &inR,
		uint64_t const &inmmask) : arginfo(inarg), AF(rAF), S(inS), R(inR), mmask(inmmask), entries(), wb_no(0), ending(false), overlaphist(), adapterhist(), adptcnt(0), paircnt(0) {}

	~WorkBlock() {}
};

libmaus2::bambam::AdapterFilter::unique_ptr_type AF;
struct WorkBlockComparator
{
	bool operator()(WorkBlock const &A, WorkBlock const &B) const
	{
    	    	return A.wb_no > B.wb_no;
	}

	bool operator()(WorkBlock const *A, WorkBlock const *B) const
	{
    	    	return A->wb_no > B->wb_no;
	}
};


struct StatsBlock
{
	uint64_t alcnt;
	uint64_t adptcnt;
	uint64_t lastalcnt;
	uint64_t paircnt;
	uint64_t orphcnt;
	libmaus2::util::Histogram overlaphist;
	libmaus2::util::Histogram adapterhist;

	StatsBlock() : alcnt(0), adptcnt(0), lastalcnt(std::numeric_limits<uint64_t>::max()),
	    paircnt(0), orphcnt(0)
	{
	}
};



/*
    Heap with a passed in mutex to limit size.
    Control of size is in AdapterFindControl below.
*/
struct LimitedLockedHeap
{
	typedef LimitedLockedHeap this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	pthread_mutex_t *mutex;
	pthread_cond_t  *space_available;
	pthread_cond_t  result_available;

	std::priority_queue<WorkBlock*,std::vector<WorkBlock*>, WorkBlockComparator> Q;
	uint64_t max_size;

	LimitedLockedHeap(pthread_mutex_t *rmutex, pthread_cond_t  *rspace_available, const uint64_t rmax_size) :
	    mutex(rmutex), space_available(rspace_available), Q(), max_size(rmax_size)
	{
	    	pthread_cond_init(&result_available, NULL);
	}

	~LimitedLockedHeap()
	{
	    	pthread_cond_destroy(&result_available);
	}

	struct ScopeLimitedMutexLock
	{
		pthread_mutex_t * mutex;

		ScopeLimitedMutexLock(pthread_mutex_t *rmutex) : mutex(rmutex)
		{
			if ( pthread_mutex_lock(mutex) != 0 )
			{
				int const error = errno;
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "ScopeLimitedMutexLock failed pthread_mutex_lock " << strerror(error) << std::endl;
				lme.finish();
				throw lme;
			}
		}

		~ScopeLimitedMutexLock()
		{
			if ( pthread_mutex_unlock(mutex) != 0 )
			{
				int const error = errno;
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "ScopeLimitedMutexLock failed pthread_mutex_unlock " << strerror(error) << std::endl;
				lme.finish();
				throw lme;
			}
		}
	};

	bool empty()
	{
		ScopeLimitedMutexLock llock(mutex);
		return Q.size() == 0;
	}

	bool full()
	{
		ScopeLimitedMutexLock llock(mutex);
		return (Q.size() >= max_size);
	}

	void push(WorkBlock *v)
	{
		ScopeLimitedMutexLock llock(mutex);

		Q.push(v);
    		pthread_cond_signal(&result_available);
	}


	WorkBlock* pop()
	{
		ScopeLimitedMutexLock llock(mutex);
		WorkBlock *p = Q.top();
		Q.pop();

		if (max_size && Q.size() < max_size)
		{
		    	pthread_cond_signal(space_available);
		}

		return p;
	}


	WorkBlock* top()
	{
		ScopeLimitedMutexLock llock(mutex);
		return Q.top();
	}
};


void adapterListMatch(
	libmaus2::autoarray::AutoArray<char> & Aread,
	libmaus2::util::PushBuffer<libmaus2::bambam::AdapterOffsetStrand> & AOSPB,
	libmaus2::fastx::AutoArrayWordPutObject<uint64_t> & AAWPO,
	libmaus2::bambam::BamAlignment & algn,
	libmaus2::bambam::AdapterFilter const & AF,
	int const verbose,
	uint64_t const adpmatchminscore, // = getDefaultMatchMinScore(),
	double const adpmatchminfrac, // = getDefaultMatchMinFrac(),
	double const adpmatchminpfrac, // = getDefaultMatchMinPFrac()
	uint64_t const reflen,
	double const pA,
	double const pC,
	double const pG,
	double const pT
)
{
    	uint64_t len;

    	if ( algn.isReverse() )
	{
	        len = algn.decodeReadRC(Aread);
	}
	else
	{
	    	len = algn.decodeRead(Aread);
	}

	uint8_t const * const ua = reinterpret_cast<uint8_t const *>(Aread.begin());

	bool const matched = AF.searchAdapters(
		ua, len,
		2 /* max mismatches */,
		AOSPB,
		adpmatchminscore /* min score */,
		adpmatchminfrac /* minfrac */,
		adpmatchminpfrac /* minpfrac */,
		AAWPO,
		0 /* verbose */
	);

	if ( matched )
	{
		std::sort(AOSPB.begin(),AOSPB.end(),libmaus2::bambam::AdapterOffsetStrandMatchStartComparator());
		uint64_t const clip = len - AOSPB.begin()[0].getMatchStart();

		uint64_t fA, fC, fG, fT;
		AF.getAdapterMatchFreqs(len,AOSPB.begin()[0],fA,fC,fG,fT);

		// confidence that this is not a random event
		double const randconfid = kmerPoisson(reflen,fA,fC,fG,fT,0/*n*/,pA,pC,pG,pT);

		if ( verbose > 1 )
		{
			std::cerr << "\n" << std::string(80,'-') << "\n\n";

			std::cerr << "read length " << len << " clip " << clip << " randconfig=" << randconfid << " fA=" << fA << " fC=" << fC << " fG=" << fG << " fT=" << fT << std::endl;

			for ( libmaus2::bambam::AdapterOffsetStrand const * it = AOSPB.begin(); it != AOSPB.end(); ++it )
				AF.printAdapterMatch(ua,len,*it);
		}

		algn.putAuxNumber("as",'i',clip);
		algn.putAuxString("aa",AF.getAdapterName(AOSPB.begin()[0]));
		algn.putAuxNumber("af",'f',AOSPB.begin()[0].pfrac);
		algn.putAuxNumber("ar",'f',randconfid);
	}
}


namespace libmaus2
{
    	namespace parallel
	{

    		struct FreeLists;

    		struct AdapterFindWorkPackage : public SimpleThreadWorkPackage
		{
			typedef AdapterFindWorkPackage this_type;
			typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
			typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

			FreeLists *FL;
			WorkBlock *wb;

			AdapterFindWorkPackage() : FL(NULL), wb(NULL) {}

			AdapterFindWorkPackage(
	    		    FreeLists *rFL,
			    WorkBlock *rwb,
	    		    uint64_t const rpriority,
			    uint64_t const rdispatcherid,
			    uint64_t const rpackageid,
			    uint64_t const rsubid) : SimpleThreadWorkPackage(rpriority,rdispatcherid,rpackageid, rsubid), FL(rFL), wb(rwb)
			{
			}

			virtual char const *getPackageName() const
			{
	    		    return "AdapterFindWorkPackage";
			}

		};


    		struct AdapterFindWritePackage : public SimpleThreadWorkPackage
		{
		       typedef AdapterFindWritePackage this_type;
		       typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
		       typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

		       FreeLists *FL;
    		       const ::libmaus2::util::ArgInfo *arginfo;
    		       const std::string *headertext;
		       StatsBlock *sb;

		       AdapterFindWritePackage() : FL(NULL), arginfo(NULL), headertext(NULL) {}

		       AdapterFindWritePackage(
	    		   FreeLists *rFL,
			   const ::libmaus2::util::ArgInfo *rarginfo,
			   const std::string *rheadertext,
			   StatsBlock *rsb,
	       		   uint64_t const rpriority,
			   uint64_t const rdispatcherid,
			   uint64_t const rpackageid) : SimpleThreadWorkPackage(rpriority,rdispatcherid,rpackageid),
		    	    	    			FL(rFL), arginfo(rarginfo), headertext(rheadertext), sb(rsb)
		       {
		       }

		       virtual char const *getPackageName() const
		       {
	    		   return "AdapterFindWritePackage";
		       }
		};


		struct FreeLists
		{
		    SimpleThreadPoolWorkPackageFreeList<AdapterFindWorkPackage> findWorkPackages;
		    SimpleThreadPoolWorkPackageFreeList<AdapterFindWritePackage> findWritePackages;
    	    	    LimitedLockedHeap write_queue;

		    FreeLists(pthread_mutex_t *rmutex, pthread_cond_t  *rspace_available, const uint64_t rmax_size) : write_queue(rmutex, rspace_available, rmax_size)
		    {
		    }
		};


		struct AdapterFindWorkPackageDispatcher : public SimpleThreadWorkPackageDispatcher
		{
			virtual void dispatch(
	    		    SimpleThreadWorkPackage *P,
			    SimpleThreadPoolInterfaceEnqueTermInterface &)
			{
	    			AdapterFindWorkPackage *WP = dynamic_cast<AdapterFindWorkPackage *>(P);
				WorkBlock *wb = WP->wb;
				libmaus2::bambam::AdapterFilter::shared_ptr_type AF = wb->AF;

				uint64_t const adpmatchminscore  = wb->arginfo.getValue<uint64_t>("adpmatchminscore",getDefaultMatchMinScore());
				double   const adpmatchminfrac   = wb->arginfo.getValue<double>("adpmatchminfrac",getDefaultMatchMinFrac());
				double   const adpmatchminpfrac  = wb->arginfo.getValue<double>("adpmatchminpfrac",getDefaultMatchMinPFrac());
				uint64_t const reflen = wb->arginfo.getValue<uint64_t>("reflen",getDefaultRefLen());
				double   const pA = wb->arginfo.getValue<double>("pA",getDefaultpA());
				double   const pC = wb->arginfo.getValue<double>("pC",getDefaultpC());
				double   const pG = wb->arginfo.getValue<double>("pG",getDefaultpG());
				double   const pT = wb->arginfo.getValue<double>("pT",getDefaultpT());

				int const verbose = wb->arginfo.getValue<int>("verbose",getDefaultVerbose());
				int const clip    = wb->arginfo.getValue<int>("clip",getDefaultClip());

				uint64_t const seedlength =
					std::min(
						static_cast<unsigned int>((8*sizeof(uint64_t))/3),
						std::max(wb->arginfo.getValue<unsigned int>("SEED_LENGTH",getDefaultSEED_LENGTH()),1u));
				// maximum mismatch rate in overlap
				double const mismatchrate = std::min(100u,wb->arginfo.getValue<unsigned int>("PCT_MISMATCH",getDefaultPCT_MISMATCH()))/100.0;
				// maximum number of mismatches in seed
				unsigned int const maxseedmismatches = wb->arginfo.getValue<unsigned int>("MAX_SEED_MISMATCHES",static_cast<unsigned int>(std::ceil(seedlength * mismatchrate)));
				// minimum length of overlap between reads
				uint64_t const minoverlap = wb->arginfo.getValue<uint64_t>("MIN_OVERLAP",getDefaultMIN_OVERLAP());
				// maximum number of adapter bases to be compared
				uint64_t const almax = wb->arginfo.getValue<uint64_t>("ADAPTER_MATCH",getDefaultADAPTER_MATCH());

				libmaus2::autoarray::AutoArray<char> Aread;
				libmaus2::util::PushBuffer<libmaus2::bambam::AdapterOffsetStrand> AOSPB;
				libmaus2::autoarray::AutoArray<char> seqs[2];

				libmaus2::bambam::BamAuxFilterVector auxfilter;
				auxfilter.set("a3");
				auxfilter.set("ah");

    				// clip adapter variables
				libmaus2::autoarray::AutoArray<char> CR;
				libmaus2::autoarray::AutoArray<char> CQ;
				libmaus2::bambam::BamSeqEncodeTable const seqenc;
				libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;
				libmaus2::bambam::BamAlignment::D_array_type T;
				libmaus2::fastx::AutoArrayWordPutObject<uint64_t> working_data;

    				for (std::vector<HoldAlign *>::iterator itr = wb->entries.begin(); itr != wb->entries.end(); itr++) {
	    				// find adapters in given list
					adapterListMatch(Aread,AOSPB,working_data,(*itr)->algns[0], *AF, verbose, adpmatchminscore,adpmatchminfrac,adpmatchminpfrac,reflen,pA,pC,pG,pT);

					if (!(*itr)->pair) {
					    if (clip) {
		    				clipAdapters((*itr)->algns[0],CR,CQ,seqenc,cigop,T);
					    }

					    continue;
					}

					adapterListMatch(Aread,AOSPB, working_data, (*itr)->algns[1], *AF, verbose, adpmatchminscore,adpmatchminfrac,adpmatchminpfrac,reflen,pA,pC,pG,pT);

					if (!((*itr)->algns[0].isRead1() && (*itr)->algns[1].isRead2())) {
						std::cerr << "[D] warning: reads are not in the correct order" << std::endl;

						if ( clip )
						{
							clipAdapters((*itr)->algns[0],CR,CQ,seqenc,cigop,T);
							clipAdapters((*itr)->algns[1],CR,CQ,seqenc,cigop,T);
						}

						continue;
					}


					// are the reads both non empty?
					if ( !((*itr)->algns[0].getLseq() && (*itr)->algns[1].getLseq()) )
					{
						std::cerr << "[D] warning: empty read" << std::endl;

						if ( clip )
						{
							clipAdapters((*itr)->algns[0],CR,CQ,seqenc,cigop,T);
							clipAdapters((*itr)->algns[1],CR,CQ,seqenc,cigop,T);
						}

						continue;
					}

					wb->paircnt++;

					unsigned int const rev0 = (*itr)->algns[0].isReverse() ? 1 : 0;
					unsigned int const rev1 = (*itr)->algns[1].isReverse() ? 1 : 0;

					/*
					 * Whether a sequence needs to be reverse complemented for
					 * adaptor detection depends on its strand and in what orientation
					 * it is stored.
					 */

					if ( !rev0 )
					{
    	    					(*itr)->algns[0].decodeRead(seqs[0]);

						if ( !rev1 )
						{
			        			(*itr)->algns[1].decodeReadRC(seqs[1]);
						}
						else
						{
		    	        			(*itr)->algns[1].decodeRead(seqs[1]);
						}
					}
					else
					{
						(*itr)->algns[0].decodeReadRC(seqs[0]);

						if ( !rev1 )
						{
			        			(*itr)->algns[1].decodeReadRC(seqs[1]);
						}
						else
						{
		    	        			(*itr)->algns[1].decodeRead(seqs[1]);
						}
					}

					uint64_t const l0 = (*itr)->algns[0].getLseq();
					uint64_t const l1 = (*itr)->algns[1].getLseq();
					uint64_t const lm = std::min(l0,l1);
					/*
					 * local seed length; this may be smaller then the
					 * global seed length if one of the too reads
					 * is shorter than the seed length
					 */
					uint64_t const lseedlength = std::min(seedlength,lm);
					uint64_t const lseedmask = libmaus2::math::lowbits(3*lseedlength);

					// compute seed from last lseedlength bases of second read
					uint64_t seed = 0;
					uint8_t const * p = reinterpret_cast<uint8_t const *>(seqs[1].begin() + l1);
					for ( uint64_t i = 0; i < lseedlength; ++i )
					{
						seed <<= 3;
						seed  |= wb->S[*(--p)];
					}

					// compute mask of length lseedlength-1 from back of first read
					uint64_t query = 0;
					uint8_t const * const qe = reinterpret_cast<uint8_t const *>(seqs[0].begin());
					uint8_t const * q = qe + l0;
					uint64_t matchpos = l0-lseedlength;

					for ( uint64_t i = 0; i < lseedlength-1; ++i )
					{
						query <<= 3;
						query |= wb->R[*(--q)];
					}


					// try to find seed in first read
					do
					{
						// add next base (step backward from end)
						query <<= 3;
						query |= wb->R[*(--q)];
						query &= lseedmask;

						// compute number of mismatches
						uint64_t dif = (query ^ seed);
						dif = (dif | (dif >> 1) | (dif >> 2)) & wb->mmask;
						unsigned int const difcnt = libmaus2::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(dif);

						#if defined(DIFCNTDEBUG)
						unsigned int debdifcnt = 0;
						for ( unsigned int i = 0; i < lseedlength; ++i )
							if ( wb->S[*(seqs[0].begin()+matchpos+i)] != wb->R[*(seqs[1].begin()+l1-lseedlength + i)] )
								debdifcnt++;
						assert ( debdifcnt == difcnt );
						#endif

						// if the seed matches, then look at the rest
						if ( difcnt <= maxseedmismatches )
						{
							// end of total match on read 1
							uint64_t const end0 = matchpos + lseedlength;
							// end of total match on read 2
							uint64_t const end1 = l1;
							// position of seed on read 1
							uint64_t const seedp0 = end0 - lseedlength;
							// position of seed on read 2
							uint64_t const seedp1 = end1 - lseedlength;
							// length of match between read 1 and read 2 (no indels allowed)
							uint64_t const restoverlap = std::min(seedp0,seedp1);

							// if length of overlap is sufficient
							if ( lseedlength + restoverlap >= minoverlap )
							{
								// iterator range of non seed overlap on read 1
								uint8_t const * check0  = reinterpret_cast<uint8_t const *>(seqs[0].begin()+seedp0-restoverlap);
								uint8_t const * check0e = check0 + restoverlap;
								// start iterator of non seed overlap on read 2
								uint8_t const * check1  = reinterpret_cast<uint8_t const *>(seqs[1].begin()+seedp1-restoverlap);

								// maximum number of mismatches allowed
								uint64_t const maxmis = (restoverlap+lseedlength) * mismatchrate;

								// compute number of mismatches
								uint64_t nummis = difcnt;
								while ( nummis <= maxmis && check0 != check0e  )
									if ( wb->S[*check0++] != wb->R[*check1++] )
										nummis++;

								// if match is within the maximal number of mismatches
								if ( check0 == check0e && nummis <= maxmis )
								{
									// update overlap histogram
									wb->overlaphist(lseedlength + restoverlap);

									// adapter length for read 1
									uint64_t const al0 = (l0-end0);
									// adapter length for read 2
									uint64_t const al1 = end1-(restoverlap+lseedlength);
									// number of bases compared in this case
									uint64_t const alcmp = std::min(almax,std::min(al0,al1));

									// adapter range on read 1
									char const * ap0 = seqs[0].begin()+end0;
									char const * ap0e = ap0 + alcmp;
									// adapter start on read 2 (actually adapter end, we will read it backwards)
									char const * ap1 = seqs[1].begin()+al1;
									// number of mismatches on adapter (up to a length of almax)
									unsigned int aldif = 0;

									// calculate mismatches in adapter
									while ( ap0 != ap0e )
										if ( wb->S[*(ap0++)] != wb->R[libmaus2::fastx::invertUnmapped(*(--ap1))] )
											aldif++;

									// if number of mismatches is within the tolerated range (0 at this time)
									if ( (ap0 == ap0e) && ! aldif )
									{

										if ( verbose > 1 )
										{
											std::cerr
												<< "[V2] overlap: " << (*itr)->algns[0].getName()
												<< " mismatchrate=" <<
													nummis << "/" << (restoverlap+lseedlength) << "=" <<
													static_cast<double>(nummis)/(restoverlap+lseedlength)
												<< " al0=" << al0
												<< " al1=" << al1
												<< " aldif=" << aldif
												<< " alcmp=" << alcmp
												<< "\n" <<
												"[V2] " << std::string(
													seqs[0].begin()+seedp0-restoverlap,
													seqs[0].begin()+end0) << "\n" <<
												"[V2] " << std::string(
													seqs[1].begin()+seedp1-restoverlap,
													seqs[1].begin()+end1) << "\n";

											std::cerr << "[V2] assumed adapter on read 1 ["<<al0<<"]: "
												<< std::string(
													seqs[0].begin()+end0,
													seqs[0].begin()+l0
												) << std::endl;

											std::cerr << "[V2] assumed adapter on read 2 ["<<al1<<"]: "
												<< libmaus2::fastx::reverseComplementUnmapped(std::string(
													seqs[1].begin(),
													seqs[1].begin()+al1)) << std::endl;
										}


										(*itr)->algns[0].filterOutAux(auxfilter);
										(*itr)->algns[1].filterOutAux(auxfilter);

										(*itr)->algns[0].putAuxNumber("ah",'i',1);
										(*itr)->algns[1].putAuxNumber("ah",'i',1);
										(*itr)->algns[0].putAuxNumber("a3",'i',al0);
										(*itr)->algns[1].putAuxNumber("a3",'i',al1);

										wb->adptcnt ++;
										wb->adapterhist(lseedlength + restoverlap);

										break;
									}
								}
							}
						}

						--matchpos;
					} while ( q != qe );

					if ( clip )
					{
						clipAdapters((*itr)->algns[0],CR,CQ,seqenc,cigop,T);
						clipAdapters((*itr)->algns[1],CR,CQ,seqenc,cigop,T);
					}
				}

				WP->FL->write_queue.push(wb);
				WP->FL->findWorkPackages.returnPackage(WP);
			}
		};


		struct AdapterFindWritePackageDispatcher : public SimpleThreadWorkPackageDispatcher
		{
			virtual void dispatch(
	    		    SimpleThreadWorkPackage *P,
			    SimpleThreadPoolInterfaceEnqueTermInterface &tpi)
			{
	    			AdapterFindWritePackage *WP = dynamic_cast<AdapterFindWritePackage *>(P);

				int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(WP->arginfo->getValue<int>("level",getDefaultLevel()));

				// add PG line to header
				std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
					*(WP->headertext),
					"bamadapterfind", // ID
					"bamadapterfind", // PN
					WP->arginfo->commandline, // CL
					::libmaus2::bambam::ProgramHeaderLineSet(*(WP->headertext)).getLastIdInChain(), // PP
					std::string(PACKAGE_VERSION) // VN
				);
				// construct new header
				::libmaus2::bambam::BamHeader uphead(upheadtext);

				/*
				 * start index/md5 callbacks
				 */
				std::string const tmpfilenamebase = WP->arginfo->getValue<std::string>("tmpfile",WP->arginfo->getDefaultTmpFileName());
				std::string const tmpfileindex = tmpfilenamebase + "_index";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

				std::string md5filename;
				std::string indexfilename;

				std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > cbs;
				::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;

				if ( WP->arginfo->getValue<unsigned int>("md5",getDefaultMD5()) )
				{
   	    	    	        	if ( libmaus2::bambam::BamBlockWriterBaseFactory::getMD5FileName(*(WP->arginfo)) != std::string() )
    	    	    		    		md5filename = libmaus2::bambam::BamBlockWriterBaseFactory::getMD5FileName(*(WP->arginfo));
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
				if ( WP->arginfo->getValue<unsigned int>("index",getDefaultIndex()) )
				{
			    		if ( libmaus2::bambam::BamBlockWriterBaseFactory::getIndexFileName(*(WP->arginfo)) != std::string() )
    	    	    		    		md5filename = libmaus2::bambam::BamBlockWriterBaseFactory::getIndexFileName(*(WP->arginfo));

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

				::libmaus2::bambam::BamWriter::unique_ptr_type writer(new ::libmaus2::bambam::BamWriter(std::cout,uphead,level,Pcbs));

 				bool keep_running = true;
				int expected = 0;

    				libmaus2::util::Histogram overlaphist_total;
	    			libmaus2::util::Histogram adapterhist_total;

				while (!WP->FL->write_queue.empty() || keep_running)
				{
					if (!WP->FL->write_queue.empty())
					{
	    					WorkBlock *wb = WP->FL->write_queue.top();

						 if (wb->wb_no == expected)
						 {
							expected++;
							WP->FL->write_queue.pop();

							// add the stats
							WP->sb->adptcnt += wb->adptcnt;
							WP->sb->paircnt += wb->paircnt;
							WP->sb->overlaphist.merge(wb->overlaphist);
							WP->sb->adapterhist.merge(wb->adapterhist);

							for (std::vector<HoldAlign *>::iterator itr = wb->entries.begin(); itr != wb->entries.end(); ++itr)
							{
								(*itr)->algns[0].serialise(writer->getStream());

								if ((*itr)->pair)
								{
		    							(*itr)->algns[1].serialise(writer->getStream());
								}

								delete (*itr);
							}

							if (wb->ending) keep_running = false;

							delete wb;
						}
					}
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

				WP->FL->findWritePackages.returnPackage(WP);
				tpi.terminate();
			}
		};


    		struct AdapterFindControl
		{
			SimpleThreadPool &STP;

			FreeLists free;

			AdapterFindWorkPackageDispatcher  work_dispatcher;
			AdapterFindWritePackageDispatcher write_dispatcher;

    			pthread_mutex_t mutex;
        		pthread_cond_t  space_available;

			uint64_t const max_queued_packages;

			uint64_t work_dispatcher_id;
			uint64_t write_dispatcher_id;

			AdapterFindControl(SimpleThreadPool &rSTP, const uint64_t rmax_size) : STP(rSTP), free(&mutex, &space_available, rmax_size),
	    	    		    max_queued_packages(STP.getNumThreads() * 2)
			{
	    			work_dispatcher_id  = STP.getNextDispatcherId();
				STP.registerDispatcher(work_dispatcher_id, &work_dispatcher);

				write_dispatcher_id = STP.getNextDispatcherId();
				STP.registerDispatcher(write_dispatcher_id, &write_dispatcher);

    	    			pthread_mutexattr_t attr;
    	    			pthread_mutexattr_init(&attr);

    	    			#if defined(HAVE_PTHREAD_MUTEX_RECURSIVE_NP)
    	    			pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE_NP);
    	    			#elif defined(HAVE_PTHREAD_MUTEX_RECURSIVE)
    	    			pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    	    			#else
    	    			#error "Neither PTHREAD_MUTEX_RECURSIVE_NP nor PTHREAD_MUTEX_RECURSIVE are supported for setting up recursive POSIX mutexes"
    	    			#endif

 	    			pthread_mutex_init(&mutex, &attr);
				pthread_cond_init(&space_available, NULL);

    	    			pthread_mutexattr_destroy(&attr);
			}

			~AdapterFindControl()
			{
	    			pthread_mutex_destroy(&mutex);
				pthread_cond_destroy(&space_available);
			}

			AdapterFindWorkPackage *getWorkPackage(void)
			{
	    			return free.findWorkPackages.getPackage();
			}

			AdapterFindWritePackage *getWritePackage(void)
			{
	    		    	return free.findWritePackages.getPackage();
			}

			void enque(SimpleThreadWorkPackage *package)
			{
	    			pthread_mutex_lock(&mutex);

				while ((STP.Q.getFillState() > max_queued_packages) || free.write_queue.full())
				{
				    	pthread_cond_wait(&space_available, &mutex);
				}

	    			STP.enque(package);

				pthread_mutex_unlock(&mutex);
			}

			void join(void) {
	    		    STP.join();
			}
		};
	 }
}


int bamadapterfind(::libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	int const max_threads     = arginfo.getValue<int>("threads", getDefaultThreads());
	int const max_write_queue = arginfo.getValue<int>("max_write_queue", getDefaultWriteQueue());
	int const max_work_reads  = arginfo.getValue<int>("max_work_reads", getDefaultWorkReads());

	libmaus2::bambam::AdapterFilter::shared_ptr_type AF;

	if ( arginfo.hasArg("adaptersbam") )
	{
		libmaus2::aio::InputStreamInstance adapterCIS(arginfo.getUnparsedValue("adaptersbam","adapters.bam"));
		libmaus2::bambam::AdapterFilter::shared_ptr_type tAF(
                                new libmaus2::bambam::AdapterFilter(adapterCIS,12 /* seed length */)
                        );
		AF = tAF;
	}
	else
	{
		std::string const builtinAdapters = libmaus2::bambam::BamDefaultAdapters::getDefaultAdapters();
		std::istringstream builtinAdaptersStr(builtinAdapters);
		libmaus2::bambam::AdapterFilter::shared_ptr_type tAF(
                                new libmaus2::bambam::AdapterFilter(builtinAdaptersStr,12 /* seed length */)
                        );
		AF = tAF;
	}


 	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);

	::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & bamdec = *ppdec;
	::libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
 	std::string const headertext(header.text);
	libmaus2::bambam::BamAlignment & inputalgn = bamdec.getAlignment();

	libmaus2::autoarray::AutoArray<uint8_t> S(256,false);
	std::fill(S.begin(),S.end(),4);
	S['a'] = S['A'] = 0;
	S['c'] = S['C'] = 1;
	S['g'] = S['G'] = 2;
	S['t'] = S['T'] = 3;

	libmaus2::autoarray::AutoArray<uint8_t> R(256,false);
	std::fill(R.begin(),R.end(),5);
	R['a'] = R['A'] = 0;
	R['c'] = R['C'] = 1;
	R['g'] = R['G'] = 2;
	R['t'] = R['T'] = 3;

	static uint64_t const mmask =
		(1ull << 0) |
		(1ull << 3) |
		(1ull << 6) |
		(1ull << 9) |
		(1ull << 12) |
		(1ull << 15) |
		(1ull << 18) |
		(1ull << 21) |
		(1ull << 24) |
		(1ull << 27) |
		(1ull << 30) |
		(1ull << 33) |
		(1ull << 36) |
		(1ull << 39) |
		(1ull << 42) |
		(1ull << 45) |
		(1ull << 48) |
		(1ull << 51) |
		(1ull << 54) |
		(1ull << 57) |
		(1ull << 60) |
		(1ull << 63);

    	libmaus2::parallel::SimpleThreadPool TP(max_threads < 2 ? 2 : max_threads); // 2 thread minimum
	libmaus2::parallel::AdapterFindControl AFC(TP, max_write_queue);

	StatsBlock stats;

	libmaus2::parallel::AdapterFindWritePackage *wp = AFC.getWritePackage();
	*wp = libmaus2::parallel::AdapterFindWritePackage(&(AFC.free), &arginfo, &headertext, &stats, 0, AFC.write_dispatcher_id, 0);
	AFC.enque(wp);

	libmaus2::parallel::AdapterFindWorkPackage *fp;

	WorkBlock *workb      = NULL;
	int wbnum             = 0;
	const int read_max    = max_work_reads;
	int read_count        = 0;
	uint64_t const bshift = libmaus2::math::ilog(libmaus2::math::nextTwoPow(arginfo.getValue<int>("mod",getDefaultMod())));
	int const verbose     = arginfo.getValue<int>("verbose",getDefaultVerbose());

	while ( (bamdec.readAlignment()) )
	{
		if ( verbose && ( (stats.alcnt >> bshift) != (stats.lastalcnt >> bshift) ) )
		{
			std::cerr << "[V]\t" << stats.alcnt << "\t" << stats.paircnt << "\t" << stats.adptcnt << std::endl;
			stats.lastalcnt = stats.alcnt;
		}

	    	if (workb == NULL)
		{
		    	workb = new WorkBlock(arginfo, AF, S, R, mmask);
		    	workb->wb_no = wbnum++;
		}

		HoldAlign *hold = new HoldAlign();

		hold->algns[0].swap(inputalgn);
		stats.alcnt++;
		read_count++;

		if ( !hold->algns[0].isPaired())
		{
		    	hold->pair = false;
		}
		else
		{
		    	if (!bamdec.readAlignment())
			{
			    	stats.orphcnt++;
			}
			else
			{
				// if next is not paired or name does not match
				if ( (! inputalgn.isPaired()) || strcmp(hold->algns[0].getName(), inputalgn.getName()) )
				{
		    			hold->pair = false;
		    			read_count++;
					stats.orphcnt++;
					bamdec.putback();
				}
				else
				{
					// sanity checks
					assert ( hold->algns[0].isPaired() );
					assert ( inputalgn.isPaired() );
					assert ( strcmp(hold->algns[0].getName(),inputalgn.getName()) == 0 );
					hold->algns[1].swap(inputalgn);
					read_count++;
					stats.alcnt++;
					hold->pair = true;
				}
			}
		}

		workb->entries.push_back(hold);

		if (read_count >= read_max)
		{
		        fp = AFC.getWorkPackage();
		        *fp = libmaus2::parallel::AdapterFindWorkPackage(&(AFC.free), workb, 0, AFC.work_dispatcher_id, workb->wb_no, workb->wb_no);
		        AFC.enque(fp);

    	    	    	workb = NULL;
			read_count = 0;
		}
	}

	// mark the last work block
	if (workb)
	{
	    	workb->ending = true;
	}
	else
	{
		// special ending job
		workb = new WorkBlock(arginfo, AF, S, R, mmask);
		workb->wb_no = wbnum++;
		workb->ending = true;
	}

    	fp  = AFC.getWorkPackage();
	*fp = libmaus2::parallel::AdapterFindWorkPackage(&(AFC.free), workb, 0/* workb->wb_no */, AFC.work_dispatcher_id, workb->wb_no, workb->wb_no);
    	AFC.enque(fp);

	AFC.join();

	if ( verbose )
		std::cerr << "[V] processed " << stats.alcnt << std::endl;

	if ( verbose )
		std::cerr << "[V] paircnt=" << stats.paircnt << " adptcnt=" << stats.adptcnt << " alcnt=" << stats.alcnt << " orphcnt=" << stats.orphcnt << std::endl;

	std::map<uint64_t,uint64_t> const overlaphistmap = stats.overlaphist.get();
	std::map<uint64_t,uint64_t> const adapterhistmap = stats.adapterhist.get();

	uint64_t totOverlaps = 0;
	uint64_t totAdapters = 0;

	// print histogram
	for ( std::map<uint64_t,uint64_t>::const_iterator ita = overlaphistmap.begin();
		ita != overlaphistmap.end(); ++ita )
	{
		// number of adapters for overlap
		uint64_t const ladpt = (adapterhistmap.find(ita->first) != adapterhistmap.end()) ?
			adapterhistmap.find(ita->first)->second : 0;

		std::cerr
			<< std::setfill('0') << std::setw(3) << ita->first << std::setw(0)
			<< " "
			<< std::setfill('0') << std::setw(7) << ita->second << std::setw(0)
			<< " "
			<< std::setfill('0') << std::setw(7) << ladpt << std::setw(0)
			<< std::endl;

		totOverlaps += ita->second;
		totAdapters += ladpt;
	}

	double const pctOverlaps = (static_cast<double>(totOverlaps)/stats.paircnt)*100.0;
	double const pctAdapters = (static_cast<double>(totAdapters)/stats.paircnt)*100.0;

	std::cerr << "Processed " << std::setw(7) << std::setfill('0') << stats.paircnt << std::setw(0) << " read pairs." << std::endl;
	std::cerr << "Found     " << std::setw(7) << std::setfill('0') << totOverlaps << std::setw(0) << " overlaps (" << pctOverlaps << ")" << std::endl;
	std::cerr << "Found     " << std::setw(7) << std::setfill('0') << totAdapters << std::setw(0) << " adapters (" << pctAdapters << ")" << std::endl;

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
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam2::Licensing::formatNumber(getDefaultMod())+"]>", "print progress every mod'th line (if verbose>0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adaptersbam=<[]>", "list of adapters/primers stored in a BAM file (use internal list if not given)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "clip=<["+::biobambam2::Licensing::formatNumber(getDefaultClip())+"]>", "clip off adapters (see bamadapterclip program)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "SEED_LENGTH=<["+::biobambam2::Licensing::formatNumber(getDefaultSEED_LENGTH())+"]>", "length of seed for matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "PCT_MISMATCH=<["+::biobambam2::Licensing::formatNumber(getDefaultPCT_MISMATCH())+"]>", "maximum percentage of mismatches in matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "MAX_SEED_MISMATCHES=<[SEED_LENGTH*PCT_MISMATCH]>", "maximum number of mismatches in seed (up to 2)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "MIN_OVERLAP=<["+::biobambam2::Licensing::formatNumber(getDefaultMIN_OVERLAP())+"]>", "minimum overlap between mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "ADAPTER_MATCH=<["+::biobambam2::Licensing::formatNumber(getDefaultADAPTER_MATCH())+"]>", "maximum adapter match check" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adpmatchminscore=<["+::biobambam2::Licensing::formatNumber(getDefaultMatchMinScore())+"]>", "minimum score for adapter list matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adpmatchminfrac=<["+::biobambam2::Licensing::formatFloatingPoint(getDefaultMatchMinFrac())+"]>", "minimum fraction for adapter list matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adpmatchminpfrac=<["+::biobambam2::Licensing::formatFloatingPoint(getDefaultMatchMinPFrac())+"]>", "minimum fraction of overlap for adapter list matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reflen=<["+::biobambam2::Licensing::formatNumber(getDefaultRefLen())+"]>", "length of reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pA=<["+::biobambam2::Licensing::formatFloatingPoint(getDefaultpA())+"]>", "relative frequency of base A in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pC=<["+::biobambam2::Licensing::formatFloatingPoint(getDefaultpC())+"]>", "relative frequency of base C in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pG=<["+::biobambam2::Licensing::formatFloatingPoint(getDefaultpG())+"]>", "relative frequency of base G in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pT=<["+::biobambam2::Licensing::formatFloatingPoint(getDefaultpT())+"]>", "relative frequency of base T in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				V.push_back ( std::pair<std::string,std::string> ( "threads=<["+::biobambam2::Licensing::formatNumber(getDefaultThreads())+"]>", "number of threads to find with" ) );
				V.push_back ( std::pair<std::string,std::string> ( "max_write_queue=<["+::biobambam2::Licensing::formatNumber(getDefaultWriteQueue())+"]>", "max number of waiting write packages" ) );
				V.push_back ( std::pair<std::string,std::string> ( "max_work_reads=<["+::biobambam2::Licensing::formatNumber(getDefaultWorkReads())+"]>", "number of reads per package" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return bamadapterfind(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
