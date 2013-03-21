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

#include <bambam/PutFastQBase.hpp>
#include <bambam/OutputBuffer.hpp>

#include <fstream>
#include <iomanip>
#include <queue>
#include <stdexcept>
#include <libmaus/util/unordered_map.hpp>
#include <csignal>

#include <bambam/BamBamConfig.hpp>


#if defined(BAMBAM_HAVE_SAMTOOLS)
#include <bambam/BamFile.hpp>

#if defined(HAVE_BAM_H)
#include <bam.h>
#include <sam.h>
#elif defined(HAVE_SAMTOOLS_BAM_H)
#include <samtools/bam.h>
#include <samtools/sam.h>
#else
#error "Required bam.h header not available."
#endif
#endif

#define EVAHASH

#include <libmaus/aio/BufferedOutput.hpp>
#include <libmaus/bitio/getBit.hpp>
#include <libmaus/bitio/putBit.hpp>
#include <libmaus/bitbtree/bitbtree.hpp>
#include <libmaus/fastx/FastQReader.hpp>
#if defined(EVAHASH)
#include <libmaus/hashing/hash.hpp>
#endif
#include <libmaus/math/bitsPerNum.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/WriteableString.hpp>
#include <libmaus/fastx/SpaceTable.hpp>
#include <libmaus/util/TempFileNameGenerator.hpp>
#include <libmaus/bambam/BamDecoder.hpp>

#include <bambam/bamToFastQ.hpp>

static std::vector < std::string > removeList;

void removeTempFiles()
{
	for ( uint64_t i = 0; i < removeList.size(); ++i )
		remove ( removeList[i].c_str() );
}

#if defined(__APPLE__) || defined(__FreeBSD__)
typedef sig_t sighandler_t;
#endif

static sighandler_t siginthandler = 0;
static sighandler_t sigtermhandler = 0;
static sighandler_t sighuphandler = 0;
static sighandler_t sigpipehandler = 0;

void sigIntHandler(int arg)
{
	removeTempFiles();
	if ( siginthandler )
	{
		siginthandler(arg);
	}
	else
	{
		signal(SIGINT,SIG_DFL);
		raise(SIGINT);
	}
}

void sigHupHandler(int arg)
{
	removeTempFiles();
	if ( sighuphandler )
	{
		sighuphandler(arg);
	}
	else
	{
		signal(SIGHUP,SIG_DFL);
		raise(SIGHUP);
	}
}

void sigTermHandler(int arg)
{
	removeTempFiles();
	if ( sigtermhandler )
	{
		sigtermhandler(arg);
	}
	else
	{
		signal(SIGTERM,SIG_DFL);
		raise(SIGTERM);
	}
}

void sigPipeHandler(int arg)
{
	removeTempFiles();
	if ( sigpipehandler )
	{
		sigpipehandler(arg);
	}
	else
	{
		signal(SIGPIPE,SIG_DFL);
		raise(SIGPIPE);
	}
}

namespace bambam
{
	std::string getUnmatchedFilename(::libmaus::util::ArgInfo const & arginfo, std::string const pid)
	{
		return arginfo.getValue<std::string>("T",
			std::string("unmatched")
			+
			pid
			+
			std::string(".fq"));
	}

	struct SingleOutputBuffer : public PutFastQBase
	{
		typedef SingleOutputBuffer this_type;
		typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

		::libmaus::autoarray::AutoArray<uint8_t> outbuf;
		uint8_t * opa;
		uint8_t * opc;
		uint8_t * ope;
		uint64_t written;
		::libmaus::fastx::SpaceTable const ST;
			
		SingleOutputBuffer()
		: outbuf(), opa(0), opc(0), ope(0), written(0), ST()
		{}
		virtual ~SingleOutputBuffer() {}
		
		void write(std::ostream & out) const
		{
			out.write ( reinterpret_cast<char const *>(opa), opc-opa );
		}

		uint64_t freeSpace()
		{
			return ope-opc;
		}
		
		inline void prepare(uint64_t)
		{
		}
			
		virtual void flush()
		{
			opc = opa;
		}

		void checkSpace(uint64_t const outlen)
		{
			// buffer overflow?
			if ( freeSpace() < outlen )
			{
				if ( outlen > outbuf.size() )
				{
					::libmaus::autoarray::AutoArray<uint8_t> newbuf(outlen);	
					std::copy( outbuf.begin(), outbuf.end(), newbuf.begin() );
					
					outbuf = newbuf;
					opa = outbuf.begin();
					opc = opa;
					ope = outbuf.end();
				}
			}
			
			assert ( freeSpace() >= outlen );		
		}
		
		inline uint64_t put(CompactFastQEntry const & CFQ)
		{
			flush();
			
			// length of alignment in fastq file		
			uint64_t const outlen = CFQ.entrylen;
			
			checkSpace(outlen);

			uint64_t const putpos = opc-opa;
			
			for ( uint64_t i = 0; i < outlen; ++i )
				*(opc++) = CFQ.T[i];
				
			written += outlen;
				
			return putpos;
		}
		
		#if defined(BAMBAM_HAVE_SAMTOOLS)
		inline uint64_t put(bam1_t const * alignment)
		{
			flush();
		
			bam1_core_t const * alignment_core = &(alignment->core);
			// get query string
			uint8_t const * seq = bam1_seq(alignment);
			// length of query string
			unsigned int const seqlen = alignment_core->l_qseq;
			// flags
			bam_flag_type const flags = alignment_core->flag;
			// name of query
			char const * qname = bam1_qname(alignment);
			unsigned int const qnamelen = strlen(qname);

			// length of alignment in fastq file		
			uint64_t const outlen = getFastqEntryLength(qnamelen,seqlen,flags);
			
			checkSpace(outlen);

			uint64_t const putpos = opc-opa;

			putRead(qname,qnamelen,seq,bam1_qual(alignment),seqlen,flags,opc,ST);
			
			written += outlen;
			
			return putpos;
		}
		#endif

		template<typename ptr_type>
		uint64_t put(ptr_type p, unsigned int const len, bool const = true /* first */)
		{
			flush();
		
			checkSpace(len);
		
			uint64_t const putpos = opc-opa;
		
			ptr_type pe = p + len;
			while ( p != pe )
				*(opc++) = *(p++);
			
			written += len;
				
			return putpos;
		}
	};

	template<typename buffer_pair_type>
	struct SortingOutputBuffer : public OutputBuffer
	{
		typedef OutputBuffer base_type;
		typedef SortingOutputBuffer<buffer_pair_type> this_type;
		typedef typename ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

		::libmaus::autoarray::AutoArray<uint64_t> startbuf;
		uint64_t sortfill;
		std::vector < ::libmaus::fastx::FastInterval > FIV;
		
		buffer_pair_type & matchedoutput;
		uint64_t pairs;
		
		SortingOutputBuffer(
			uint64_t bufsize, 
			std::ostream & rout,
			buffer_pair_type & rmatchedoutput
		)
		: OutputBuffer(bufsize,rout),
		  startbuf( ((bufsize+1+63)/(8*sizeof(uint64_t))), true ),
		  sortfill(0), 
		  matchedoutput(rmatchedoutput), 
		  pairs(0)
		{
		}
		
		struct Comparator
		{
			uint8_t const * const T;
			
			Comparator(uint8_t const * rT) : T(rT) {}
			
			bool operator()(
				std::pair<uint32_t,uint32_t> const & A,
				std::pair<uint32_t,uint32_t> const & B) const
			{
				return strncmp(
					reinterpret_cast<char const *>(T+A.first),
					reinterpret_cast<char const *>(T+B.first),
					std::min(A.second,B.second)) < 0;
			}
		};
		
		std::string getKey(uint32_t p, uint32_t l) const
		{
			uint32_t pp = p;
			while ( l && !isspace(outbuf[pp]) )
				pp++, l--;
			return std::string(outbuf.get()+p,outbuf.get()+pp);
		}
		
		std::string getKeyPrefix(uint32_t p, uint32_t l, std::string & fullkey) const
		{
			fullkey = getKey(p,l);
			uint32_t lastslash = fullkey.size();
			for ( uint32_t i = 0; i < fullkey.size(); ++i )
				if ( fullkey[i] == '/' )
					lastslash = i;
			return fullkey.substr(0,lastslash);
		}
		
		static bool isFirst(std::string const & k)
		{
			return
				k.size() >= 2
				&&
				k[k.size()-2] == '/'
				&&
				k[k.size()-1] == '1';
		}
		static bool isSecond(std::string const & k)
		{
			return
				k.size() >= 2
				&&
				k[k.size()-2] == '/'
				&&
				k[k.size()-1] == '2';
		}
		static bool isFirstAndSecond(std::string const & k0, std::string const & k1)
		{
			return 
				(isFirst(k0)&&isSecond(k1)) ||
				(isFirst(k1)&&isSecond(k0));
		}
		
		virtual void flush()
		{
			if ( opc != opa )
			{
				// set final bit
				::libmaus::bitio::putBit(startbuf.get(),opc-opa,1);
				assert ( ::libmaus::bitio::getBit(startbuf.get(),0) );

				// number of entries
				uint64_t nument = 0;
				// previous start
				uint64_t prev = 0;
				// start and length of entries
				std::vector < std::pair<uint32_t,uint32_t> > P;

				for ( uint64_t i = 1; i < static_cast<uint64_t>((opc-opa) + 1); ++i )
					if ( ::libmaus::bitio::getBit(startbuf.get(),i) )
					{
						P.push_back(std::pair<uint32_t,uint32_t>(prev,i-prev));
						nument++;
						prev = i;
					}
				
				assert ( sortfill == nument );

				// sort entries
				std::sort ( P.begin(), P.end(), Comparator(outbuf.get()) );
				
				::libmaus::aio::BufferedOutput<uint8_t> OB(out,256*1024);

				uint64_t outputbytes = 0, outputal = 0;
				for ( uint64_t i = 0 ; i < P.size(); )
				{
					std::string fullk0, fullk1;
					std::string k0, k1;
					
					if ( 
						i+1 < P.size()
						&&
						((k0=getKeyPrefix(P[i+0].first,P[i+0].second,fullk0))
						==
						 (k1=getKeyPrefix(P[i+1].first,P[i+1].second,fullk1)))
						&&
						isFirstAndSecond(fullk0,fullk1)
					)
					{
						// std::cerr << k0 << "==" << k1 << std::endl;

						if ( isFirst(fullk0) )
						{
							matchedoutput.putPair(
								outbuf.get() + P[i+0].first, P[i+0].second,
								outbuf.get() + P[i+1].first, P[i+1].second
							);
							matchedoutput.finishedPair();
						}
						else
						{
							matchedoutput.putPair(
								outbuf.get() + P[i+1].first, P[i+1].second,
								outbuf.get() + P[i+0].first, P[i+0].second
							);
							matchedoutput.finishedPair();
						}
						
						pairs++;
						
						i += 2;
					}
					else
					{
						OB.put(outbuf.get() + P[i].first, P[i].second );
						outputbytes += P[i].second;
						outputal ++;
						i++;
					}
				}
					
				OB.flush();
				written += OB.totalwrittenbytes;

				uint64_t const readlow = FIV.size() ? FIV.back().high : 0;
				uint64_t const readhigh = readlow + outputal;
				uint64_t const fileoffset = FIV.size() ? FIV.back().fileoffsethigh : 0;
				uint64_t const fileoffsethigh = fileoffset + outputbytes;
				
				if ( readhigh-readlow )
					FIV.push_back(::libmaus::fastx::FastInterval(readlow,readhigh,fileoffset,fileoffsethigh,0/*numsyms*/,0/*minlen*/,0/*maxlen*/));

				// erase start bit vector
				std::fill(startbuf.begin(),startbuf.end(),0ull);	
				sortfill = 0;

				opc = opa;			
			}
		}
		
		uint64_t put(CompactFastQEntry const & CFQ)
		{
			uint64_t const putpos = base_type::put(CFQ);
			// put reads in buffer
			// std::cerr << "Putting len " << outlen << " at " << opc-opa << std::endl;
			::libmaus::bitio::putBit(startbuf.get(),putpos,1);
			sortfill++;
			
			return putpos;
		}
		
		#if defined(BAMBAM_HAVE_SAMTOOLS)
		uint64_t put(bam1_t const * alignment)
		{
			uint64_t const putpos = base_type::put(alignment);

			// std::cerr << "Putting len " << outlen << " at " << opc-opa << std::endl;
			::libmaus::bitio::putBit(startbuf.get(),putpos,1);
			sortfill++;
			
			return putpos;
		}
		#endif

		template<typename ptr_type>
		uint64_t put(ptr_type p, unsigned int const len, bool const first /* first write of entry */)
		{
			uint64_t const putpos = base_type::put(p,len /* ,first */);

			if ( first )
			{
				// std::cerr << "Putting len " << len << " at " << opc-opa << std::endl;
				::libmaus::bitio::putBit(startbuf.get(),putpos,1);
				sortfill++;
			}

			return putpos;
		}
	};

	typedef ::libmaus::fastx::FastQReader orphan_reader_type;
	typedef orphan_reader_type::unique_ptr_type orphan_reader_ptr_type;
	typedef orphan_reader_type::pattern_type orphan_pattern_type;

	struct OrphanQueueElement
	{
		orphan_pattern_type pattern;	
		uint64_t readerid;
		
		OrphanQueueElement() {}
		OrphanQueueElement(orphan_pattern_type const & rpattern, uint64_t const & rreaderid)
		: pattern(rpattern), readerid(rreaderid) {}
		
		bool operator<(OrphanQueueElement const & O) const
		{
			return pattern.sid > O.pattern.sid;
		}
		
		bool isFirst() const
		{
			return 
				pattern.sid.size() >= 2 && 
				pattern.sid[pattern.sid.size()-2] == '/' &&
				pattern.sid[pattern.sid.size()-1] == '1';
		}
		bool isSecond() const
		{
			return 
				pattern.sid.size() >= 2 && 
				pattern.sid[pattern.sid.size()-2] == '/' &&
				pattern.sid[pattern.sid.size()-1] == '2';
		}
		
		static std::string getKeyPrefix(std::string const & key)
		{

			uint32_t keylen = 0;
			while ( keylen < key.size() && (!isspace(key[keylen])) )
				keylen++;
			uint32_t lastslash = keylen;
			for ( uint32_t i = 0; i < keylen; ++i )
				if ( key[i] == '/' )
					lastslash = i;
			return key.substr(0,lastslash);		
		}
		
		std::string getKeyPrefix() const
		{
			return getKeyPrefix(pattern.sid);
		}
	};

	template<typename buffer_type>
	struct OutputBufferPair : protected std::pair < buffer_type * , buffer_type * >
	{
		uint64_t pairs;

		OutputBufferPair(buffer_type * buffera, buffer_type * bufferb)
		: std::pair < buffer_type * , buffer_type * >(buffera,bufferb), pairs(0) {}
		virtual ~OutputBufferPair() {}
		
		buffer_type * getFirst()  { return std::pair < buffer_type * , buffer_type * >::first; }
		buffer_type * getSecond() { return std::pair < buffer_type * , buffer_type * >::second; }
		buffer_type const * getFirst() const  { return std::pair < buffer_type * , buffer_type * >::first; }
		buffer_type const * getSecond() const { return std::pair < buffer_type * , buffer_type * >::second; }
		
		void flush()
		{
			getFirst()->flush();
			getSecond()->flush();
		}

		bool endsOn(std::string const & a, std::string const & b)
		{
			return a.size() >= b.size() && a.substr(a.size()-b.size(),b.size()) == b;
		}

		template<typename ptr_type_a>
		void putFirst(
			ptr_type_a p, unsigned int const len
		)
		{
			getFirst()->prepare(len);
			getFirst()->put(p,len);

		}
		template<typename ptr_type_b>
		void putSecond(
			ptr_type_b p, unsigned int const len
		)
		{
			getSecond()->prepare(len);
			getSecond()->put(p,len);
		}
		
		
		template<typename ptr_type_a, typename ptr_type_b>
		void putPair(
			ptr_type_a pa, unsigned int const lena, 
			ptr_type_b pb, unsigned int const lenb
		)
		{
			putFirst(pa,lena);
			putSecond(pb,lenb);
		}
		uint64_t getWritten() const
		{
			if ( getFirst() != getSecond() )
				return getFirst()->written + getSecond()->written;
			else
				return getFirst()->written;
		}
		virtual void finishedPair()
		{
			pairs++;
		}
	};

	struct FilePairContainer
	{
		typedef OutputBuffer::unique_ptr_type buffer_ptr_type;

		// matched buffer /1
		buffer_ptr_type Pmatchedbuffera;
		OutputBuffer * const matchedbuffera;

		// matched buffer /2
		buffer_ptr_type Pmatchedbufferb;
		OutputBuffer * const matchedbufferb;

		FilePairContainer(std::ostream & rmatchedaostr, std::ostream & rmatchedbostr)
		:
			Pmatchedbuffera(new OutputBuffer(256*1024,rmatchedaostr)),
			matchedbuffera(Pmatchedbuffera.get()),
			Pmatchedbufferb(new OutputBuffer(256*1024,rmatchedbostr)),
			matchedbufferb(Pmatchedbufferb.get())
		{}

		FilePairContainer(std::ostream & rmatchedostr)
		:
			Pmatchedbuffera(new OutputBuffer(256*1024,rmatchedostr)),
			matchedbuffera(Pmatchedbuffera.get()),
			Pmatchedbufferb(),
			matchedbufferb(Pmatchedbuffera.get())
		{}	
	};

	struct FileOutputPair : public FilePairContainer, public OutputBufferPair<OutputBuffer>
	{
		typedef FileOutputPair this_type;
		typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

		FileOutputPair(std::ostream & rmatchedaostr, std::ostream & rmatchedbostr)
		: FilePairContainer(rmatchedaostr,rmatchedbostr), OutputBufferPair<OutputBuffer>(FilePairContainer::matchedbuffera,FilePairContainer::matchedbufferb)
		{}

		FileOutputPair(std::ostream & rmatchedaostr)
		: FilePairContainer(rmatchedaostr), OutputBufferPair<OutputBuffer>(FilePairContainer::matchedbuffera,FilePairContainer::matchedbufferb)
		{}
		
		static bool aligned()
		{
			return true;
		}
	};

	template<typename unmatched_buffer_type>
	struct NamedReadBuffer : public ::libmaus::fastx::SpaceTable
	{
		typedef NamedReadBuffer<unmatched_buffer_type> this_type;
		typedef typename ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

		// character data
		::libmaus::autoarray::AutoArray<uint8_t> B;
		// current offset in character data
		uint64_t p;
		// start of data positions dictionary
		::libmaus::bitbtree::BitBTree<8,8> dict;
			
		// output buffer for unmatched entries
		unmatched_buffer_type & unmatchedbuffer;

		// hash bits (log_2 of table length)
		unsigned int namehashbits;
		// length of hash table (power of two)
		uint64_t namehashlength;
		// length of hash table minus 1
		uint64_t namehashmask;
		// hash table (max value for unset)
		typedef uint32_t hash_pos_type;
		::libmaus::autoarray::AutoArray<hash_pos_type> namehash;
		
		// offset of key in stored entries
		static unsigned int const putkeyoffset = 1;
		// length of a number in bytes
		static unsigned int const numlen = 2;
		
		// portion buffer
		::libmaus::autoarray::AutoArray<uint8_t> portionbuffer;

		#if ! defined(EVAHASH)
		// hash object
		std::hash<std::string> const shash;
		#endif

		NamedReadBuffer(
			uint64_t const rbufsize,
			unmatched_buffer_type & runmatchedbuffer
		)
		:
			B(rbufsize,false), p(0), dict(rbufsize,false), 
			unmatchedbuffer(runmatchedbuffer),
			namehashbits( std::max(static_cast<int>(::libmaus::math::bitsPerNum(rbufsize)) - 2, 0)  ),
			namehashlength(1ul << namehashbits),
			namehashmask ( namehashlength -1 ),
			namehash ( namehashlength, false ),
			portionbuffer( )
			#if ! defined(EVAHASH)
			,shash ()
			#endif
		{
			std::fill ( namehash.begin(), namehash.end(), std::numeric_limits<hash_pos_type>::max() );
		}
		
		void checkPortionBuffer(uint64_t const s)
		{
			if ( s > portionbuffer.size() )
			{
				::libmaus::autoarray::AutoArray<uint8_t> newportionbuffer(s);
				std::copy(portionbuffer.begin(),portionbuffer.end(),newportionbuffer.begin());
				portionbuffer = newportionbuffer;
			}
		}
		
		void rehash()
		{
			namehashbits += 1;
			namehashlength <<= 1;
			namehashmask = namehashlength-1;
			namehash = ::libmaus::autoarray::AutoArray<hash_pos_type>(namehashlength,false);
			std::fill ( namehash.begin(), namehash.end(), std::numeric_limits<hash_pos_type>::max() );

			for ( uint64_t i = 0; i < dict.count1(); ++i )
			{
				uint64_t const q = dict.select1(i);
				namehash [ getHashPos(q) ] = q;
			}
		}
		
		/* free space between current position and next stored entry */
		uint64_t space(uint64_t & no) const
		{
			// do we have an entry?
			if ( dict.hasNext1() )
			{
				// yes, position of next entry
				no = dict.next1(p);
				
				// compute space between current pointer and next used position
				uint64_t r;
				
				if ( no < p )
					r = no+B.size()-p;
				else
					r = no-p;
					
				return r;
			}
			else
			{
				// no entries, buffer is completely unused
				return B.size();
			}
		}
		
		uint32_t getNumber(uint64_t q) const
		{
			// std::cerr << "Getting number from " << q << std::endl;
			uint32_t v = 0;
		
			if ( q + numlen - 1 < B.size() )
				for ( unsigned int i = 0; i < numlen; ++i )
				{
					v <<= 8;
					v |= static_cast<uint32_t>(B[q++]);
				}
			else		
				for ( unsigned int i = 0; i < numlen; ++i )
				{
					v <<= 8;
					v |= static_cast<uint32_t>(B[(q++)%B.size()]);
				}
			
			return v;
		}
		
		void putNumber(uint32_t const q)
		{
			assert ( p < B.size() );
			unsigned int shift = 8*(numlen-1);
			
			if ( p+numlen < B.size() )
				for ( unsigned int i = 0; i < numlen; ++i, shift-=8 )
					B[p++] = (q>>shift)&0xFF;
			else
				for ( unsigned int i = 0; i < numlen; ++i, shift-=8, p %= B.size() )
					B[p++] = (q>>shift)&0xFF;
			
			assert ( p < B.size() );
		}
		

		uint32_t computeHash(uint8_t const * c, uint64_t const keylen) const
		{
			#if defined(EVAHASH)
			return ::libmaus::hashing::EvaHash::hash(c,keylen) & namehashmask;
			#else
			char const * ca = reinterpret_cast<char const *>(c);
			char const * ce = ca + keylen;
			return shash(std::string(ca,ce)) & namehashmask;
			#endif
		}
		
		// write and erase entry starting at position no in buffer
		template<bool use_buffer, typename buffer_type>
		void writePortion(
			uint64_t no, 
			buffer_type & out, 
			uint32_t const h
		)
		{
			dict.set(no,false);
			namehash [ h ] = std::numeric_limits<hash_pos_type>::max();

			uint64_t towrite = getNumber(no);
			// skip meta data
			no += numlen;
			no %= B.size();

			if ( out.freeSpace() < towrite )
				out.flush();
			out.prepare(towrite);

			// std::cerr << "towrite: " << towrite << std::endl;

			if ( ! use_buffer )
			{
				bool first = true;

				while ( towrite )
				{
					assert ( no < B.size() );
					uint32_t const portion = std::min(towrite,static_cast<uint64_t>(B.size()-no));
					out.put ( reinterpret_cast<char const *>(B.begin() + no) , portion, first );
					towrite -= portion;
					no += portion;
					no %= B.size();
					first = false;
				}		
			}
			else
			{
				checkPortionBuffer(towrite);
			
				uint8_t * outp = portionbuffer.get();
				uint64_t const towritec = towrite;
				
				while ( towrite )
				{
					assert ( no < B.size() );
					uint32_t const portion = std::min(towrite,static_cast<uint64_t>(B.size()-no));

					uint8_t const * s = B.begin() + no;
					uint8_t const * se = s + portion;
					
					while ( s != se )
						*(outp++) = *(s++);

					towrite -= portion;
					no += portion;
					no %= B.size();
				}		

				out.put ( reinterpret_cast<char const *>(portionbuffer.get()) , towritec, true /* first */ );
			}
		}
		
		void writePortionUnmatched(uint64_t const no, uint32_t const h)
		{
			writePortion<false>(no,unmatchedbuffer,h);
		}

		template<typename buffer_type>
		void writePortionMatched(
			uint64_t const no, 
			uint32_t const h,
			// bool const first
			buffer_type & out,
			bool const aligned
		)
		{
			if ( aligned )
				writePortion<true>(no,out,h);
			else
				writePortion<false>(no,out,h);
		}
		
		// write all entries present to the unmatched buffer and flush it
		void flush()
		{
			while ( dict.hasNext1() )
			{
				uint64_t const no = dict.next1(p);
				writePortionUnmatched (no, getHashPos(no));		
			}
		}

		// compute length of hash key outside buffer
		// key is up to last slash before first whitespace or up to end of buffer
		unsigned int computeKeyLength(uint8_t const * c, uint64_t n)
		{
			uint8_t const * const ca = c;
			uint8_t const * const ce = ca + n;
			uint8_t const * cp = ce;
			
			// search for space or end
			while ( c != ce && nospacetable[*c] )
			{
				if ( *c == '/' )
					cp = c;
				++c;
			}

			return cp-ca;
		}
		
		// compute length of in buffer string key
		unsigned int getHashLength(uint64_t q) const
		{
			// skip length of entry and offset (@ character)
			q = (q + numlen + putkeyoffset) % B.size();
			
			bool keyset = false;
			unsigned int offset = 0;
			unsigned int keylen = 0;
			while ( q < B.size() && nospacetable[B[q]] )
			{
				if ( B[q] == '/' )
				{
					keyset = true;
					keylen = offset;
				}
				offset++;
				q++;
			}
			if ( q == B.size() )
			{
				q = 0;

				while ( q < B.size() && nospacetable[B[q]] )
				{
					if ( B[q] == '/' )
					{
						keyset = true;
						keylen = offset;
					}
					offset++;
					q++;
				}
			}
			
			if ( ! keyset )
				keylen = offset;
				
			return keylen;
		}
		
		// get hash of entry at position q
		uint32_t getHashPos(uint64_t const q) const
		{
			uint32_t const keylen = getHashLength(q);
			if ( q + numlen + putkeyoffset + keylen <= B.size() )
				return computeHash(B.begin()+q+numlen+putkeyoffset,keylen);
			else
			{
				::libmaus::autoarray::AutoArray<uint8_t> H(keylen,false);
				for ( unsigned int i = 0; i < keylen; ++i )
					H[i] = B[(q+numlen+putkeyoffset+i) % B.size()];
				return computeHash(H.begin(),keylen);
			}
		}
		
		// compute key of given string
		bool haveKey(uint8_t const * c, uint64_t const keylen, uint32_t & kp, uint32_t & h, unsigned int & end) const
		{
			h = computeHash(c,keylen);
			uint32_t keypos = namehash[h];
			if ( keypos == std::numeric_limits<hash_pos_type>::max() )
				return false;

			kp = keypos;
			keypos += numlen /* length of entry */+putkeyoffset;
			
			if ( keypos + keylen + 3 <= B.size() )
			{
				bool eq = true;
				for ( unsigned int i = 0; i < keylen+1; ++i )
					eq = eq && (B[keypos++] == *(c++));
				eq = 
					eq &&
					((B[keypos] == '1' && *c == '2')
					||
					(B[keypos] == '2' && *c == '1'));
				end = B[keypos];
				keypos++;
				c++;
				eq = eq && spacetable[B[keypos]] && spacetable[*c];
				return eq;
			}
			else
			{
				bool eq = true;
				for ( unsigned int i = 0; i < keylen+1; ++i )
					eq = eq && (B[(keypos++)%B.size()] == *(c++));
				eq = 
					eq &&
					((B[keypos % B.size()] == '1' && *c == '2')
					||
					 (B[keypos % B.size()] == '2' && *c == '1'));
				end = B[keypos%B.size()];
				keypos++;
				c++;
				eq = eq && spacetable[B[keypos % B.size()]] && spacetable[*c];
				return eq;		
			}
		}
		
		// put a new entry
		void put(uint8_t const * c, uint64_t n, uint32_t const h)
		{
			// check if entry fits in buffer
			assert ( n+numlen <= B.size() );
			
			/**
			 * if hash slot is in use then free it by writing out the cor. entry
			 **/	
			if ( namehash[h] != std::numeric_limits<hash_pos_type>::max() )
				writePortionUnmatched(namehash[h],h);

			/**
			 * write out entries until we have enough sufficient free space ahead
			 **/
			uint64_t no = 0;
			while ( space(no) < (n+numlen) )
				writePortionUnmatched(no,getHashPos(no));

			/* register new entry */
			namehash[h] = p;
			dict.set(p,true);
			putNumber(n);

			// copy data to buffer
			while ( n )
			{
				uint32_t const portion = std::min(n,static_cast<uint64_t>(B.size()-p));
				std::copy ( c, c+portion, B.begin()+p );
				p = (p + portion)%B.size();
				n -= portion;
				c += portion;
			}
			
			#if defined(BAMTOFASTQDEBUG)
			assert ( h == getHashPos(namehash[h]) );
			#endif
		}
	};

	struct FileSet
	{
		typedef FileSet this_type;
		typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

		typedef std::ofstream out_type;
		typedef ::libmaus::util::unique_ptr<out_type>::type out_ptr_type;

		typedef FileOutputPair pair_buffer_type;
		typedef pair_buffer_type::unique_ptr_type pair_buffer_ptr_type;

		typedef SortingOutputBuffer< pair_buffer_type > unmatched_buffer_type;
		typedef unmatched_buffer_type::unique_ptr_type unmatched_buffer_ptr_type;

		typedef NamedReadBuffer<unmatched_buffer_type> read_buffer_type;
		typedef ::libmaus::util::unique_ptr<read_buffer_type>::type read_buffer_ptr_type;

		::libmaus::util::ArgInfo const & arginfo;
		std::string const unmatchedfilename;
		std::string const orphanafilename, orphanbfilename;

		out_ptr_type unmatchedostr, matchedaostr, matchedbostr, singleostr;
		pair_buffer_ptr_type PB;
		read_buffer_ptr_type NRB;
		out_ptr_type Porphan1, Porphan2;
		out_type * orphan1;
		out_type * orphan2;
		uint64_t orphans;
		
		// single end buffer
		OutputBuffer::unique_ptr_type singlebuffer;

		// output buffer for unmatched entries
		unmatched_buffer_ptr_type unmatchedbuffer;

		FileSet(::libmaus::util::ArgInfo const & rarginfo)
		: arginfo(rarginfo), 
		  unmatchedfilename(getUnmatchedFilename(arginfo,::libmaus::util::NumberSerialisation::formatNumber(getpid(),8))),
		  orphanafilename(arginfo.getValue<std::string>("O",std::string("orphans_1.fq"))),
		  orphanbfilename(arginfo.getValue<std::string>("O2",std::string("orphans_2.fq"))),
		  orphan1(0),
		  orphan2(0),
		  orphans(0)
		{
			std::string const matchedafilename = arginfo.getValue<std::string>("F",std::string("matched_1.fq"));
			std::string const matchedbfilename = arginfo.getValue<std::string>("F2",std::string("matched_2.fq"));
			std::string const singlefilename = arginfo.getValue<std::string>("S",std::string("single.fq"));

			std::map < std::string, ::std::ostream * > filemap;

			if ( filemap.find(matchedafilename) == filemap.end() )
			{
				if ( matchedafilename == "-" )
				{				
					filemap [ matchedafilename ] = &std::cout;
				}
				else
				{
					matchedaostr = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(matchedafilename.c_str(),std::ios::binary)));
					if ( ! matchedaostr->is_open() )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "Failed to open file " << matchedafilename << ": " << strerror(errno) << std::endl;
						se.finish();
						throw se;
					}
					filemap [ matchedafilename ] = matchedaostr.get();
				}
			}
			if ( filemap.find(matchedbfilename) == filemap.end() )
			{
				if ( matchedbfilename == "-" )
				{
					filemap [ matchedbfilename ] = &std::cout;
				}
				else
				{
					matchedbostr = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(matchedbfilename.c_str(),std::ios::binary)));
					if ( ! matchedbostr->is_open() )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "Failed to open file " << matchedbfilename << ": " << strerror(errno) << std::endl;
						se.finish();
						throw se;
					}
					filemap [ matchedbfilename ] = matchedbostr.get();
				}
			}
			if ( filemap.find(singlefilename) == filemap.end() )
			{
				if ( singlefilename == "-" )
				{
					filemap [ singlefilename ] = &std::cout;
				}
				else
				{
					singleostr = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(singlefilename.c_str(),std::ios::binary)));

					if ( ! singleostr->is_open() )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "Failed to open file " << singlefilename << ": " << strerror(errno) << std::endl;
						se.finish();
						throw se;
					}

					filemap [ singlefilename ] = singleostr.get();
				}
			}

			if ( filemap.find(unmatchedfilename) == filemap.end() )
			{
				if ( unmatchedfilename == "-" )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "unmatchedfilename needs to denote a regular file, aborting." << std::endl;
					se.finish();
					throw se;
				}
				else
				{
					unmatchedostr = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(unmatchedfilename.c_str(),std::ios::binary)));

					if ( ! unmatchedostr->is_open() )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "Failed to open file " << unmatchedfilename << ": " << strerror(errno) << std::endl;
						se.finish();
						throw se;
					}

					filemap [ unmatchedfilename ] = unmatchedostr.get();
				}
			}
			else
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "unmatchedfilename is required to be unique, aborting." << std::endl;
				se.finish();
				throw se;
			}
		
			if ( matchedafilename != matchedbfilename )
				PB = UNIQUE_PTR_MOVE(pair_buffer_ptr_type(new pair_buffer_type(*(filemap.find(matchedafilename)->second),*(filemap.find(matchedbfilename)->second))));
			else
				PB = UNIQUE_PTR_MOVE(pair_buffer_ptr_type(new pair_buffer_type(*(filemap.find(matchedafilename)->second))));

			unmatchedbuffer = UNIQUE_PTR_MOVE(unmatched_buffer_ptr_type(new unmatched_buffer_type(128*1024*1024,*(filemap.find(unmatchedfilename)->second),*PB)));

			uint64_t const nrbbufsize = std::max ( static_cast<uint64_t>(128*1024), static_cast<uint64_t>(1ull << (arginfo.getValue<uint64_t>("nrbbufsizelog",20))) );
			NRB = UNIQUE_PTR_MOVE(read_buffer_ptr_type(new read_buffer_type(nrbbufsize,*unmatchedbuffer)));
			
			singlebuffer = UNIQUE_PTR_MOVE(OutputBuffer::unique_ptr_type(new OutputBuffer(256*1024,*(filemap.find(singlefilename)->second))));
		}
		~FileSet()
		{
			remove ( unmatchedfilename.c_str() );	
		}
		
		void flush()
		{
			NRB->flush();
			PB->flush();
			unmatchedbuffer->flush();
			
			if ( unmatchedostr )
				unmatchedostr->flush();
			if ( matchedaostr )
				matchedaostr->flush();
			if ( matchedbostr )
				matchedbostr->flush();
			if ( singleostr )
				singleostr->flush();
			if ( unmatchedostr )
				unmatchedostr->close();
			if ( singlebuffer )
				singlebuffer->flush();
			if ( singleostr )
				singleostr->close();
		}
		
		void checkSingleRemove()
		{
			std::string const singlefilename = arginfo.getValue<std::string>("S",std::string("single.fq"));
			
			if ( 
				(singlefilename != "-") 
				&&
				(::libmaus::util::GetFileSize::fileExists(singlefilename))
				&&
				(::libmaus::util::GetFileSize::getFileSize(singlefilename) == 0) 
			)
				remove (singlefilename.c_str());
		}

		void openOrphanFiles()
		{
			if ( orphanafilename != orphanbfilename )
			{
				Porphan1 = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(orphanafilename.c_str(),std::ios::binary)));

				if ( ! (Porphan1->is_open()) )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "Failed to open file " << orphanafilename << ": " << strerror(errno) << std::endl;
					se.finish();
					throw se;
				}

				Porphan2 = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(orphanbfilename.c_str(),std::ios::binary)));

				if ( ! (Porphan2->is_open()) )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "Failed to open file " << orphanbfilename << ": " << strerror(errno) << std::endl;
					se.finish();
					throw se;
				}

				orphan1 = Porphan1.get();
				orphan2 = Porphan2.get();
			}
			else
			{
				Porphan1 = UNIQUE_PTR_MOVE(out_ptr_type(new out_type(orphanafilename.c_str(),std::ios::binary)));

				if ( ! (Porphan1->is_open()) )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "Failed to open file " << orphanafilename << ": " << strerror(errno) << std::endl;
					se.finish();
					throw se;
				}

				orphan1 = Porphan1.get();
				orphan2 = Porphan1.get();			
			}
		}
		
		void flushOrphan()
		{
			orphan1->flush();
			orphan2->flush();
			orphan1->close();
			orphan2->close();

			if ( ! orphans )
			{
				remove ( orphanafilename.c_str() );
				remove ( orphanbfilename.c_str() );
			}
		}

		void flushMatched()
		{
			if ( PB )
				PB->flush();
				
			if ( matchedaostr )
			{
				matchedaostr->flush();
				matchedaostr->close();
				matchedaostr.reset();
			}
			if ( matchedbostr )
			{
				matchedbostr->flush();
				matchedbostr->close();
				matchedbostr.reset();
			}
		}

		void putOrphan1(orphan_pattern_type const & pattern)
		{
			*orphan1 << pattern;
		}
		void putOrphan2(orphan_pattern_type const & pattern)
		{
			*orphan2 << pattern;
		}
		void putSingle(PutFastQBase::CompactFastQEntry const & CFQ)
		{
			singlebuffer->put(CFQ);
		}
	};

	static uint64_t stringToFlag(std::string const & s)
	{
		if ( s == "PAIRED" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED;
		else if ( s == "PROPER_PAIR" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPROPER_PAIR;
		else if ( s == "UNMAP" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP;
		else if ( s == "MUNMAP" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP;
		else if ( s == "REVERSE" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE;
		else if ( s == "MREVERSE" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMREVERSE;
		else if ( s == "READ1" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1;
		else if ( s == "READ2" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2;
		else if ( s == "SECONDARY" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY;
		else if ( s == "QCFAIL" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL;
		else if ( s == "DUP" )
			return ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP;
		else
		{
			::libmaus::exception::LibMausException se;
			se.getStream() << "Unknown flag " << s << std::endl;
			se.finish();
			throw se;
		}
	}

	static uint64_t stringToFlags(std::string const & s)
	{
		std::deque<std::string> const tokens = ::libmaus::util::stringFunctions::tokenize(s,std::string(","));
		uint64_t flags = 0;
		
		for ( uint64_t i = 0; i < tokens.size(); ++i )
			flags |= stringToFlag(tokens[i]);
			
		return flags;
	}

	template<typename file_set_type>
	void processMainCollate(
		::libmaus::util::ArgInfo const & arginfo, 
		file_set_type * FS,
		bool verbose
	)
	{
		#if defined(BAMBAM_HAVE_SAMTOOLS)
		BamFile BAM(arginfo);
		samfile_t * bamfile = BAM.bamfile;
		bam1_t * alignment = BAM.alignment;
		#else
		::libmaus::bambam::BamDecoder::unique_ptr_type BAM;
		bool const trustinput = arginfo.getValue("trustinput",false);
		
		if ( arginfo.hasArg("filename") && arginfo.getValue<std::string>("filename","-") != "-" )
		{
			BAM = UNIQUE_PTR_MOVE(
				::libmaus::bambam::BamDecoder::unique_ptr_type(
					new ::libmaus::bambam::BamDecoder(arginfo.getValue<std::string>("filename","-"))
				)
			);
		}
		else
		{
			BAM = UNIQUE_PTR_MOVE(
				::libmaus::bambam::BamDecoder::unique_ptr_type(
					new ::libmaus::bambam::BamDecoder(std::cin)
				)
			);
		}
		if ( trustinput )
			BAM->disableValidation();
		
		::libmaus::bambam::BamAlignment const & alignment = BAM->alignment;
		#endif

		uint64_t const exclude = stringToFlags(arginfo.getValue<std::string>("exclude",std::string()));

		typename file_set_type::read_buffer_type * NRB = FS->NRB.get();
		
		uint64_t processed = 0;
		uint64_t pairs = 0;
		uint64_t single = 0;
		::libmaus::timing::RealTimeClock rtc; rtc.start();

		PutFastQBase::CompactFastQEntry CFQ;
		::libmaus::fastx::SpaceTable const ST;

		// read next alignment
		#if defined(BAMBAM_HAVE_SAMTOOLS)
		while ( samread(bamfile,alignment) >= 0 )
		#else
		while ( BAM->readAlignment() )
		#endif
		{
			#if defined(BAMBAM_HAVE_SAMTOOLS)
			if ( (alignment->core.flag & exclude) != 0 )
			#else
			if ( (alignment.getFlags() & exclude) != 0 )			
			#endif
			{
			
			}
			else
			{
				PutFastQBase::putCompact(alignment, CFQ, ST);
				
				if ( ! (CFQ.flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED) )
				{
					FS->putSingle(CFQ);
					single++;
				}
				else
				{
					// pointer to name
					uint8_t      const * name    = CFQ.T.get() + CFQ.namepos;
					// length of name (without /[12])
					unsigned int const   namelen = CFQ.qnamelen;

					// std::cerr << std::string(name,name+namelen+3) << std::endl;
					
					#if 0
					if ( 
						std::string(name,name+namelen) == "HWUSI-EAS610:7:14:1360:1921#0" 
						||
						std::string(name,name+namelen) == "@HWUSI-EAS610:7:14:1360:1921#0" 
					)
					{
						std::cerr << std::string(CFQ.T.get() + CFQ.namepos, CFQ.T.get() + CFQ.namepos + CFQ.qnamelen ) << std::endl;
					}
					#endif

					uint32_t kp = 0;
					// compute length of key outside name hash (up to slash)
					unsigned int hashkeylen = NRB->computeKeyLength(name,namelen);
					uint32_t h;
					// see if we have the key
					unsigned int end = 0;
					bool keyok = NRB->haveKey(name, hashkeylen, kp, h, end);
						
					if ( keyok )
					{
						if ( CFQ.flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 )
						{
							FS->PB->getFirst()->put(CFQ);
							NRB->writePortionMatched(kp,h,*(FS->PB->getSecond()),FS->PB->aligned());
							FS->PB->finishedPair();
						}
						else
						{
							NRB->writePortionMatched(kp,h,*(FS->PB->getFirst()),FS->PB->aligned());					
							FS->PB->getSecond()->put(CFQ);
							FS->PB->finishedPair();
						}
						pairs++;
					}
					// mate not yet found
					else
					{
						NRB->put(CFQ.T.get(),CFQ.entrylen,h);
					}
				}
			}
			
			processed += 1;
			if ( verbose && (processed & (1024*1024-1)) == 0 )
			{
				double const elapsed = rtc.getElapsedSeconds();
				uint64_t const written = 
					(FS->PB->getWritten()+NRB->unmatchedbuffer.written+FS->singlebuffer->written);
				std::cerr << "processed " 
					<< processed/(1024*1024)
					<< "M, "
					<< (processed/elapsed)
					<< " alignments/s"
					<< ", "
					<< written/(elapsed*1024.0*1024.0)
					<< " MB/s"
					<< std::endl;
			}
		}
		
		if ( verbose )
		{
			double const elapsed = rtc.getElapsedSeconds();
			uint64_t const written = 
				(FS->PB->getWritten()+NRB->unmatchedbuffer.written+FS->singlebuffer->written);
			std::cerr << "processed " 
				<< processed/(1024*1024)
				<< "M, "
				<< (processed/elapsed)
				<< " alignments/s"
				<< ", "
				<< written/(elapsed*1024.0*1024.0)
				<< " MB/s"
				<< std::endl;		
		}
		
		FS->flush();
		FS->checkSingleRemove();
		
		std::vector < ::libmaus::fastx::FastInterval > const & FIV = NRB->unmatchedbuffer.FIV;
		if ( FIV.size() )
		{
			if ( verbose )
				std::cerr << "Merging unmatched files...";
			::libmaus::autoarray::AutoArray < orphan_reader_ptr_type > readers(FIV.size());

			std::priority_queue < OrphanQueueElement > queue;

			for ( uint64_t i = 0; i < FIV.size(); ++i )
			{
				readers[i] = UNIQUE_PTR_MOVE(orphan_reader_ptr_type(new orphan_reader_type(FS->unmatchedfilename,FIV[i])));
				
				orphan_pattern_type pattern;
				
				if ( readers[i]->getNextPatternUnlocked(pattern) )
					queue.push(OrphanQueueElement(pattern,i));
			}
			
			FS->openOrphanFiles();

			while ( queue.size() )
			{
				OrphanQueueElement const cur = queue.top();
				queue.pop();
				
				orphan_pattern_type newpat;
				if ( readers[cur.readerid]->getNextPatternUnlocked(newpat) )
					queue.push(OrphanQueueElement(newpat,cur.readerid));
					
				if ( 
					queue.size() 
					&& 
					(queue.top().getKeyPrefix() == cur.getKeyPrefix()) 
					&&
					( (queue.top().isFirst() && cur.isSecond()) ||
					(cur.isFirst() && queue.top().isSecond()) )
				)
				{
					OrphanQueueElement const mate = queue.top();
					queue.pop();
					
					std::ostringstream curostr;
					curostr << cur.pattern;
					std::ostringstream mateostr;
					mateostr << mate.pattern;

					if ( cur.isFirst() )
					{
						FS->PB->putPair(
							curostr.str().c_str(),curostr.str().size(),
							mateostr.str().c_str(),mateostr.str().size()
						);
						FS->PB->finishedPair();
					}
					else
					{
						FS->PB->putPair(
							mateostr.str().c_str(),mateostr.str().size(),
							curostr.str().c_str(),curostr.str().size()
						);
						FS->PB->finishedPair();
					}
					
					// std::cerr << "Found new pair " << cur.pattern.plus << std::endl;

					if ( readers[mate.readerid]->getNextPatternUnlocked(newpat) )
						queue.push(OrphanQueueElement(newpat,mate.readerid));
						
					pairs++;
				}
				else
				{
					// std::cerr << "Unmatched " << cur.pattern.plus << std::endl;
					FS->orphans++;
					
					if ( cur.isSecond() )					
						FS->putOrphan2(cur.pattern);
					else
						FS->putOrphan1(cur.pattern);
				}
			}
			
			FS->flushOrphan();
			
			if ( verbose )
				std::cerr << "done." << std::endl;
		}
		
		NRB->flush();		
		FS->flushMatched();

		pairs += NRB->unmatchedbuffer.pairs;

		if ( verbose )
		{
			std::cerr << "Alignments: " << processed << std::endl;
			std::cerr << "Complete pairs: " << pairs << std::endl;
			std::cerr << "Single: " << single << std::endl;
			std::cerr << "Orphans: " << FS->orphans << std::endl;
			std::cerr << "PB pairs: " << FS->PB->pairs << std::endl;
		}
	}

	void processMainNoCollate(::libmaus::util::ArgInfo const & arginfo, bool const verbose = true)
	{
		#if defined(BAMBAM_HAVE_SAMTOOLS)
		BamFile BAM(arginfo);
		samfile_t * bamfile = BAM.bamfile;
		bam1_t * alignment = BAM.alignment;
		#else
		::libmaus::bambam::BamDecoder::unique_ptr_type BAM;
		bool const trustinput = arginfo.getValue("trustinput",false);
		
		if ( arginfo.hasArg("filename") && arginfo.getValue<std::string>("filename","-") != "-" )
		{
			BAM = UNIQUE_PTR_MOVE(
				::libmaus::bambam::BamDecoder::unique_ptr_type(
					new ::libmaus::bambam::BamDecoder(arginfo.getValue<std::string>("filename","-"))
				)
			);
		}
		else
		{
			BAM = UNIQUE_PTR_MOVE(
				::libmaus::bambam::BamDecoder::unique_ptr_type(
					new ::libmaus::bambam::BamDecoder(std::cin)
				)
			);
		}
		
		if ( trustinput )
			BAM->disableValidation();
		
		::libmaus::bambam::BamAlignment const & alignment = BAM->alignment;
		#endif
		uint64_t const exclude = stringToFlags(arginfo.getValue<std::string>("exclude",std::string()));

		::libmaus::timing::RealTimeClock rtc; rtc.start();
		std::ostream & out = std::cout;
		OutputBuffer OB(1024*1024, out);
		PutFastQBase::CompactFastQEntry CFQ;
		::libmaus::fastx::SpaceTable const ST;
		
		uint64_t processed = 0;
			
		// read next alignment
		#if defined(BAMBAM_HAVE_SAMTOOLS)
		while ( samread(bamfile,alignment) >= 0 )
		#else
		while ( BAM->readAlignment() )
		#endif
		{
			#if defined(BAMBAM_HAVE_SAMTOOLS)
			if ( (alignment->core.flag & exclude) != 0 )
			#else
			if ( (alignment.getFlags() & exclude) != 0 )
			#endif
			{
			
			}
			else
			{
				PutFastQBase::putCompact(alignment, CFQ, ST);
				OB.put(CFQ);
				if ( verbose && ((++processed) & (1024*1024-1)) == 0 )
				{
					double const sec = rtc.getElapsedSeconds();
					std::cerr 
						<< processed/(1024*1024) 
						<< " rate " << (OB.written / sec)/(1024.0*1024.0)
						<< std::endl;
				}
			}
		}

		OB.flush();
	}

	int processMain(::libmaus::util::ArgInfo const & arginfo)
	{
		bool const verbose = arginfo.getValue<uint64_t>("verbose",1);

		if ( arginfo.getValue<uint64_t>("collate",1) )
		{
			FileSet::unique_ptr_type FS(new FileSet(arginfo));
			processMainCollate(arginfo,FS.get(),verbose);
		}
		else
		{
			processMainNoCollate(arginfo,verbose);
		}

		return EXIT_SUCCESS;
	}
}
