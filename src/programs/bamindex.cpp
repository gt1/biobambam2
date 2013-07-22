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

#include <libmaus/aio/SynchronousGenericInput.hpp>
#include <libmaus/bambam/BamHeader.hpp>
#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/lz/BgzfInflate.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/util/ContainerGetObject.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <biobambam/Licensing.hpp>

bool const getDefaultVerbose() { return true; }

struct LinearChunk
{
	uint64_t refid;
	uint64_t pos;
	uint64_t alcmpstart;
	uint64_t alstart;
	
	LinearChunk()
	{
	
	}
	
	LinearChunk(
		uint64_t const rrefid,
		uint64_t const rpos,
		uint64_t const ralcmpstart,
		uint64_t const ralstart
	)
	: refid(rrefid), pos(rpos), alcmpstart(ralcmpstart), alstart(ralstart) {}
	
	bool operator<(LinearChunk const & o) const
	{
		if ( refid != o.refid )
			return refid < o.refid;
		else
			return pos < o.pos;
	}
};

std::ostream & operator<<(std::ostream & out, LinearChunk const & o)
{
	out << "LinearChunk(refid=" << o.refid << ",pos=" << o.pos << ",alcmpstart=" << o.alcmpstart << ",alstart=" << o.alstart << ")";
	return out;
}

struct BinChunk
{
	uint64_t refid;
	uint64_t bin;
	uint64_t alcmpstart;
	uint64_t alstart;
	uint64_t alcmpend;
	uint64_t alend;
	
	BinChunk()
	{
	
	}
	
	BinChunk(
		uint64_t const rrefid,
		uint64_t const rbin,
		uint64_t const ralcmpstart,
		uint64_t const ralstart,
		uint64_t const ralcmpend,
		uint64_t const ralend
	)
	: refid(rrefid), bin(rbin), alcmpstart(ralcmpstart), alstart(ralstart),
	  alcmpend(ralcmpend), alend(ralend)
	{

	}
	
	bool operator<(BinChunk const & o) const
	{
		if ( refid != o.refid )
			return refid < o.refid;
		if ( bin != o.bin )
			return bin < o.bin;
		if ( alcmpstart != o.alcmpstart )
			return alcmpstart < o.alcmpstart;
		if ( alstart != o.alstart )
			return alstart < o.alstart;
		if ( alcmpend != o.alcmpend )
			return alcmpend < o.alcmpend;
		if ( alend != o.alend )
			return alend < o.alend;
			
		return false;
	}
};

std::ostream & operator<<(std::ostream & out, BinChunk const & BC)
{
	out 
		<< "BinChunk("
		<< "refid=" << BC.refid << ","
		<< "bin=" << BC.bin << ","
		<< "alcmpstart=" << BC.alcmpstart << ","
		<< "alstart=" << BC.alstart << ","
		<< "alcmpend=" << BC.alcmpend << ","
		<< "alend=" << BC.alend << ")";
		
	return out;
}

template<typename _type>
struct MergeHeapComparator
{
	typedef _type type;

	bool operator()(
		std::pair<uint64_t,type>  const & A,
		std::pair<uint64_t,type>  const & B)
	{
		bool const A_lt_B = A.second < B.second;
		bool const B_lt_A = B.second < A.second;
		
		// if A.second != B.second
		if ( A_lt_B || B_lt_A )
			return !(A_lt_B);
		else
			return A.first > B.first;
	}
};

template<typename _element_type>
struct Buffer
{
	typedef _element_type element_type;
	libmaus::autoarray::AutoArray<element_type> A;
	element_type * const pa;
	element_type * pc;
	element_type * const pe;
	
	Buffer(uint64_t const bufsize = 8*1024)
	: A(bufsize,false), pa(A.begin()), pc(pa), pe(A.end())
	{

	}
	
	bool put(element_type const & B)
	{
		*(pc++) = B;
		return pc == pe;
	}
	
	template<typename stream_type>
	std::pair<uint64_t,uint64_t> flush(stream_type & out)
	{
		uint64_t const start = out.tellp();
		
		if ( pc != pa )
		{
			std::sort(pa,pc);
			out.write(reinterpret_cast<char const *>(pa),(pc-pa)*sizeof(element_type));
			pc = pa;
		}
		
		uint64_t const end = out.tellp();
				
		return std::pair<uint64_t,uint64_t>(start,end);
	}
	
	bool empty() const
	{
		return pc == pa;
	}
};

template<typename _element_type>
struct Merge
{
	typedef _element_type element_type;
	typedef Merge<element_type> this_type;

	std::string const & fn;
	libmaus::aio::CheckedInputStream CIS;
	std::vector< std::pair<uint64_t,uint64_t> > frags;
	
	std::priority_queue<
		std::pair< uint64_t, element_type >,
		std::vector< std::pair< uint64_t, element_type > >,
		MergeHeapComparator<element_type>
	> Q;
	
	Merge(std::string const & rfn, std::vector< std::pair<uint64_t,uint64_t> > const & rfrags)
	: fn(rfn), CIS(fn), frags(rfrags)
	{
		for ( uint64_t i = 0; i < frags.size(); ++i )
			if ( frags[i].first != frags[i].second )
			{
				CIS.seekg(frags[i].first);
				element_type BC;
				CIS.read(reinterpret_cast<char *>(&BC), sizeof(element_type));
				frags[i].first += sizeof(element_type);
				Q.push(std::pair< uint64_t, element_type >(i,BC));
			}
	}
	
	bool getNext(element_type & BC)
	{
		if ( ! Q.size() )
			return false;
		
		std::pair< uint64_t, element_type > const P = Q.top();
		BC = P.second;
		Q.pop();
		
		uint64_t const fragindex = P.first;
		if ( frags[fragindex].first < frags[fragindex].second )
		{
			CIS.seekg(frags[fragindex].first);
			element_type NBC;

			CIS.read(reinterpret_cast<char *>(&NBC), sizeof(element_type));
			Q.push(std::pair< uint64_t, element_type >(fragindex,NBC));
			
			frags[fragindex].first += sizeof(element_type);
		}
		
		return true;
	}

	static void merge(
		std::string const & infn, std::vector< std::pair<uint64_t,uint64_t> > const & frags
	)
	{
		this_type BCM(infn,frags);
		libmaus::aio::CheckedOutputStream COS(infn + ".tmp");
		element_type BC;
		while ( BCM.getNext(BC) )
			COS.write(reinterpret_cast<char const *>(&BC),sizeof(element_type));
		COS.flush();
		COS.close();
		
		rename((infn + ".tmp").c_str(),infn.c_str());
	}
};

template<typename stream_type>
std::pair<int64_t,uint64_t> countDistinctBins(stream_type & stream)
{
	int64_t prevbin = -1;
	int64_t refid = -1;
	uint64_t bins = 0;
	BinChunk BC;
	
	while ( stream.peek() != stream_type::traits_type::eof() )
	{
		stream.read(reinterpret_cast<char *>(&BC),sizeof(BinChunk));
		
		if ( refid < 0 )
			refid = BC.refid;
		
		if ( refid != static_cast<int64_t>(BC.refid) )
		{
			stream.clear();
			stream.seekg(-static_cast<int64_t>(sizeof(BinChunk)),std::ios::cur);
			break;
		}
		
		if ( static_cast<int64_t>(BC.bin) != prevbin )
			bins++;
			
		prevbin = BC.bin;
	}
		
	return std::pair<uint64_t,uint64_t>(refid,bins);
}

template<typename stream_type>
std::pair<int64_t,uint64_t> countBinChunks(stream_type & stream)
{
	int64_t bin = -1;
	int64_t refid = -1;
	uint64_t chunks = 0;

	BinChunk BC;
	
	while ( stream.peek() != stream_type::traits_type::eof() )
	{
		stream.read(reinterpret_cast<char *>(&BC),sizeof(BinChunk));

		if ( refid < 0 )
			refid = BC.refid;
		if ( bin < 0 )
			bin = BC.bin;
		
		if ( 
			refid != static_cast<int64_t>(BC.refid) 
			||
			bin != static_cast<int64_t>(BC.bin)
		)
		{
			stream.clear();
			stream.seekg(-static_cast<int64_t>(sizeof(BinChunk)),std::ios::cur);
			stream.clear();
			break;
		}
		
		chunks++;
	}
		
	return std::pair<uint64_t,uint64_t>(bin,chunks);
}

template<typename stream_type>
std::pair<int64_t,uint64_t> countLinearChunks(stream_type & stream)
{
	int64_t refid = -1;
	uint64_t cnt = 0;
	
	LinearChunk LC;
	
	while ( stream.peek() != stream_type::traits_type::eof() )
	{
		stream.read(reinterpret_cast<char *>(&LC),sizeof(LinearChunk));
		
		if ( refid == -1 )
			refid = LC.refid;
		
		// put back element, if it has a different refid	
		if ( static_cast<int64_t>(LC.refid) != refid )
		{
			stream.clear();
			stream.seekg(-static_cast<int64_t>(sizeof(LinearChunk)),std::ios::cur);
			break;
		}
		
		cnt++;
	}
	
	return std::pair<int64_t,uint64_t>(refid,cnt);
}

template<typename stream_type>
std::pair<int64_t,uint64_t> getLinearMaxPos(stream_type & stream)
{
	int64_t refid = -1;
	uint64_t maxpos = 0;
	
	LinearChunk LC;
	
	while ( stream.peek() != stream_type::traits_type::eof() )
	{
		stream.read(reinterpret_cast<char *>(&LC),sizeof(LinearChunk));
		
		if ( refid == -1 )
			refid = LC.refid;
		
		// put back element, if it has a different refid	
		if ( static_cast<int64_t>(LC.refid) != refid )
		{
			stream.clear();
			stream.seekg(-static_cast<int64_t>(sizeof(LinearChunk)),std::ios::cur);
			break;
		}
		
		maxpos = std::max(maxpos,LC.pos);
	}
	
	return std::pair<int64_t,uint64_t>(refid,maxpos);
}

template<typename stream_type>
bool peekLinearChunk(stream_type & stream, uint64_t const refid)
{
	LinearChunk LC;

	if ( stream.peek() == stream_type::traits_type::eof() )
		return false;
	
	stream.read(reinterpret_cast<char *>(&LC),sizeof(LinearChunk));
	stream.clear();
	stream.seekg(-static_cast<int64_t>(sizeof(LinearChunk)),std::ios::cur);
	
	return LC.refid == refid;
}

template<typename stream_type>
bool peekLinearChunk(stream_type & stream, uint64_t const refid, uint64_t const pos, unsigned int const posshift)
{
	LinearChunk LC;

	if ( stream.peek() == stream_type::traits_type::eof() )
		return false;
	
	stream.read(reinterpret_cast<char *>(&LC),sizeof(LinearChunk));
	stream.clear();
	stream.seekg(-static_cast<int64_t>(sizeof(LinearChunk)),std::ios::cur);
	
	return (LC.refid == refid) && ((LC.pos >> posshift)==(pos>>posshift));
}

bool checkConsisteny(std::string const & binfn, std::string const & linfn, uint64_t const numrefseq)
{
	libmaus::aio::CheckedInputStream binCIS(binfn);
	libmaus::aio::CheckedInputStream linCIS(linfn);
	
	while ( 
		binCIS.peek() != libmaus::aio::CheckedInputStream::traits_type::eof()
		&&
		linCIS.peek() != libmaus::aio::CheckedInputStream::traits_type::eof()
	)
	{
		std::pair<int64_t,uint64_t> const L = countLinearChunks(linCIS);
		std::pair<int64_t,uint64_t> const B = countDistinctBins(binCIS);
		
		if ( L.first != B.first )
			return false;
			
		if ( L.first >= static_cast<int64_t>(numrefseq) )
			return false;
	}

	return binCIS.peek() == linCIS.peek();
}

template<typename stream_type>
int64_t peekBin(stream_type & stream)
{
	BinChunk BC;
	
	if ( stream.peek() == stream_type::traits_type::eof() )
		return -1;
		
	stream.read(
		reinterpret_cast<char *>(&BC),
		sizeof(BinChunk)
	);
	
	assert ( stream.gcount() == sizeof(BinChunk) );
	
	stream.clear();
	stream.seekg(-static_cast<int64_t>(sizeof(BinChunk)),std::ios::cur);
	
	return BC.refid;
}

int bamindex(libmaus::util::ArgInfo const & arginfo, std::istream & in, std::ostream & out)
{
	bool const debug = arginfo.getValue<unsigned int>("debug",0);
	unsigned int const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());

	libmaus::lz::BgzfInflate<std::istream> rec(in);
	
	std::vector<uint8_t> headerstr;
	bool haveheader = false;
	uint64_t preblocksizes = 0;
	libmaus::autoarray::AutoArray<uint8_t> B(libmaus::lz::BgzfConstants::getBgzfMaxBlockSize());
	
	uint8_t * pa = 0; 
	uint8_t * pc = 0;
	uint64_t gcnt = 0;
	uint64_t cacct = 0; // compressed bytes we have read from file
	std::pair<uint64_t,uint64_t> rinfo;
	
	/* read and copy blocks until we have reached the end of the BAM header */
	while ( 
		(!haveheader) 
		&& 
		(gcnt=(rinfo=rec.readPlusInfo(reinterpret_cast<char *>(B.begin()),B.size())).second) 
	)
	{
		cacct += rinfo.first;
		pa = B.begin();
		pc = B.begin()+gcnt;

		std::copy ( pa, pc, std::back_insert_iterator < std::vector<uint8_t> > (headerstr) );

		try
		{
			libmaus::util::ContainerGetObject< std::vector<uint8_t> > CGO(headerstr);
			libmaus::bambam::BamHeader header;
			header.init(CGO);
			haveheader = true;
			pa += (CGO.i - preblocksizes);
		}
		catch(std::exception const & ex)
		{
			if ( debug )
				std::cerr << "[D] " << ex.what() << std::endl;
		}
	
		// need to read another block to get header, remember size of current block
		if ( ! haveheader )
			preblocksizes += gcnt;
	}
	
	if ( ! haveheader )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "EOF while reading BAM header." << std::endl;
		se.finish();
		throw se;
	}

	libmaus::util::ContainerGetObject< std::vector<uint8_t> > CGO(headerstr);
	libmaus::bambam::BamHeader header;
	header.init(CGO);

	// if header ends on a block boundary, then read the next block
	if ( pa == pc )
	{
		gcnt=(rinfo=rec.readPlusInfo(reinterpret_cast<char *>(B.begin()),B.size())).second;
		cacct += rinfo.first;
		pa = B.begin();
		pc = pa + gcnt;
	}

	/* parser state types and variables */
	enum parsestate { 
		state_reading_blocklen, 
		state_post_skip 
	};
	parsestate state = state_reading_blocklen;

	unsigned int blocklenred = 0;
	uint32_t blocklen = 0;

	uint64_t alcnt = 0;
	uint64_t alcmpstart = cacct - rinfo.first;
	uint64_t alstart = pa - B.begin();
	
	uint64_t binalcmpstart = alcmpstart;
	uint64_t binalstart = alstart;

	libmaus::bambam::BamAlignment algn;
	uint8_t * copyptr = 0;
	
	int64_t prevrefid = -1;
	int64_t prevpos = -1;
	int64_t prevbin = -1;
	
	std::string const tmpfileprefix = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const binchunktmpfilename = tmpfileprefix+".bin";
	std::string const linchunktmpfilename = tmpfileprefix+".lin";
	libmaus::util::TempFileRemovalContainer::addTempFile(binchunktmpfilename);
	libmaus::util::TempFileRemovalContainer::addTempFile(linchunktmpfilename);
	libmaus::aio::CheckedOutputStream::unique_ptr_type chunkCOS(new libmaus::aio::CheckedOutputStream(binchunktmpfilename));
	libmaus::aio::CheckedOutputStream::unique_ptr_type linCOS(new libmaus::aio::CheckedOutputStream(linchunktmpfilename));

	Buffer<BinChunk> BCB;
	std::vector< std::pair<uint64_t,uint64_t> > binchunkfrags;
	Buffer<LinearChunk> BLC;
	std::vector< std::pair<uint64_t,uint64_t> > linearchunkfrags;
	
	/* while we have alignment data blocks */
	while ( gcnt )
	{
		while ( pa != pc )
		{
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
						
						state = state_post_skip;
						
						if ( algn.D.size() < blocklen )
							algn.D = libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
						algn.blocksize = blocklen;
						copyptr = algn.D.begin();
					}
					else
					#endif
					{
						while ( pa != pc && blocklenred < sizeof(uint32_t) )
							blocklen |= static_cast<uint32_t>(*(pa++)) << ((blocklenred++)*8);

						if ( blocklenred == sizeof(uint32_t) )
						{
							state = state_post_skip;

							if ( algn.D.size() < blocklen )
								algn.D = libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
							algn.blocksize = blocklen;
							copyptr = algn.D.begin();
						}
					}
					
					break;
				/* skip data after part we modify */
				case state_post_skip:
				{
					uint32_t const skip = std::min(pc-pa,static_cast<ptrdiff_t>(blocklen));
					std::copy(pa,pa+skip,copyptr);
					copyptr += skip;
					pa += skip;
					blocklen -= skip;
					
					if ( ! blocklen )
					{
						int64_t const thisrefid = algn.getRefID();
						int64_t const thispos = algn.getPos();
						int64_t const thisbin = 
							(thisrefid >= 0 && thispos >= 0) ? algn.computeBin() : -1;
							
						if ( 
							thisrefid != prevrefid ||
							(
								(thisrefid == prevrefid) && (thispos>>14)!=(prevpos>>14)
							)
						)
						{
							if ( thisrefid >= 0 && thispos >= 0 )
							{
								LinearChunk LC(thisrefid,thispos,alcmpstart,alstart);
								if ( BLC.put(LC) )
									linearchunkfrags.push_back(BLC.flush(*linCOS));
							}

							if ( debug )
							{
								std::cerr << "new lin: " << thisrefid << "," << thispos 
									<< " alcmpstart=" << alcmpstart
									<< " alstart=" << alstart
									<< std::endl;
							}
						}	

						if ( 
							thisbin != prevbin 
							||
							thisrefid != prevrefid
						)
						{
							if ( prevbin >= 0 )
							{
								uint64_t binalcmpend = alcmpstart;
								uint64_t binalend = alstart;
							
								if ( debug )
								{
									std::cerr << "end of bin " << prevbin
										<< " on ref id " << prevrefid
										<< " binalcmpstart=" << binalcmpstart
										<< " binalstart=" << binalstart
										<< " binalcmpend=" << binalcmpend
										<< " binalend=" << binalend
										<< std::endl;
								}

								BinChunk BC(prevrefid,prevbin,binalcmpstart,binalstart,binalcmpend,binalend);
								
								if ( debug )
								{
									std::cerr << BC << std::endl;
								}
								
								if ( BCB.put(BC) )
									binchunkfrags.push_back(BCB.flush(*chunkCOS));
							}

							binalcmpstart = alcmpstart;
							binalstart = alstart;
						}
						
						if ( debug )
						{
							std::cerr 
								<< "alcnt=" << alcnt
								<< " alcmpstart=" << alcmpstart
								<< " alstart=" << alstart
								<< " refid=" << thisrefid
								<< " pos=" << thispos
								<< " bin=" << thisbin
								<< std::endl;
						}
					
						// finished an alignment, set up for next one
						state = state_reading_blocklen;
						
						blocklenred = 0;
						blocklen = 0;

						uint64_t const nextalcmpstart = cacct - rinfo.first;
						uint64_t const nextalstart = pa - B.begin();

						alcnt++;
						alcmpstart = nextalcmpstart;
						alstart = nextalstart;

						prevrefid = thisrefid;
						prevpos = thispos;
						prevbin = thisbin;
						
						if ( verbose && (alcnt % (1024*1024) == 0) )
							std::cerr << "[V] " << alcnt/(1024*1024) << std::endl;
					}
					break;
				}
			}
		}

		gcnt=(rinfo=rec.readPlusInfo(reinterpret_cast<char *>(B.begin()),B.size())).second;
		cacct += rinfo.first;
		pa = B.begin();
		pc = pa + gcnt;
	}				

	if ( prevbin >= 0 )
	{
		uint64_t binalcmpend = alcmpstart;
		uint64_t binalend = alstart;
	
		if ( debug )
		{
			std::cerr << "*end of bin " << prevbin
				<< " on ref id " << prevrefid
				<< " binalcmpstart=" << binalcmpstart
				<< " binalstart=" << binalstart
				<< " binalcmpend=" << binalcmpend
				<< " binalend=" << binalend
				<< std::endl;
		}

		BinChunk BC(prevrefid,prevbin,binalcmpstart,binalstart,binalcmpend,binalend);

		if ( BCB.put(BC) )
			binchunkfrags.push_back(BCB.flush(*chunkCOS));
	}

	if ( ! BCB.empty() )
		binchunkfrags.push_back(BCB.flush(*chunkCOS));
	if ( !BLC.empty() )
		linearchunkfrags.push_back(BLC.flush(*linCOS));

	chunkCOS->flush();
	chunkCOS->close();
	chunkCOS.reset();
	
	linCOS->flush();
	linCOS->close();
	linCOS.reset();
	
	if ( verbose )
		std::cerr << "[V] " << alcnt << std::endl;

	Merge<BinChunk>::merge(binchunktmpfilename,binchunkfrags);
	Merge<LinearChunk>::merge(linchunktmpfilename,linearchunkfrags);

	/* check consistency */
	bool const consistent = checkConsisteny(binchunktmpfilename,linchunktmpfilename,header.chromosomes.size());
	
	if ( ! consistent )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "bin index and linear index are inconsistent." << std::endl;
		se.finish();
		throw se;
	}

	/* write index */
	libmaus::aio::CheckedInputStream binCIS(binchunktmpfilename);
	libmaus::aio::CheckedInputStream linCIS(linchunktmpfilename);
		
	out.put('B');
	out.put('A');
	out.put('I');
	out.put('\1');
	libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,header.chromosomes.size());
	
	for ( uint64_t i = 0; i < header.chromosomes.size(); ++i )
	{
		if ( debug )
			std::cerr << "chromosome " << header.chromosomes[i].name << std::endl;
	
		if ( peekBin(binCIS) == static_cast<int64_t>(i) )
		{
			uint64_t const bprepos = binCIS.tellg();
			if ( debug )
				std::cerr << "bprepos=" << bprepos << std::endl;
			std::pair<int64_t,uint64_t> const P=countDistinctBins(binCIS);
			assert ( P.first == static_cast<int64_t>(i) );
			binCIS.clear();
			binCIS.seekg(bprepos,std::ios::beg);
			
			if ( debug )
				std::cerr << "Distinct bins " << P.second << std::endl;
			
			// number of distinct bins
			libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,P.second);
			
			// iterate over bins
			for ( uint64_t b = 0; b < P.second; ++b )
			{
				uint64_t const lprepos = binCIS.tellg();
				std::pair<int64_t,uint64_t> Q = countBinChunks(binCIS);

				binCIS.seekg(lprepos);
				binCIS.clear();

				if ( debug )
					std::cerr << "ref " << P.first << " dist bins " << P.second 
						<< " bin " << Q.first << " chunks " << Q.second
						<< std::endl;

				// bin
				libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,Q.first);	
				// number of chunks
				libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,Q.second);	
	
				// iterate over chunks
				for ( uint64_t c = 0; c < Q.second; ++c )
				{
					BinChunk BC;
					binCIS.read(reinterpret_cast<char *>(&BC),sizeof(BinChunk));
					
					// chunk data
					libmaus::bambam::EncoderBase::putLE<std::ostream,uint64_t>(out,(BC.alcmpstart<<16)|(BC.alstart));	
					libmaus::bambam::EncoderBase::putLE<std::ostream,uint64_t>(out,(BC.alcmpend<<16)|(BC.alend));	
					
					if ( debug )
						std::cerr << BC << std::endl;
				}
			}

			uint64_t const lprepos = linCIS.tellg();
			std::pair<int64_t,uint64_t> const Q = getLinearMaxPos(linCIS);
			assert ( Q.first == static_cast<int64_t>(i) );
			linCIS.seekg(lprepos,std::ios::beg);
			linCIS.clear();
			
			uint64_t const posperchunk = (1ull<<14);
			uint64_t const numchunks = ((Q.second+1) + posperchunk-1)/posperchunk;

			// number of linear chunks bins
			libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,numchunks);	
			
			for ( uint64_t c = 0; c < numchunks; ++c )
			{
				if ( peekLinearChunk(linCIS,P.first,c*posperchunk,14) )
				{
					LinearChunk LC;
					linCIS.read(reinterpret_cast<char *>(&LC),sizeof(LinearChunk));
					libmaus::bambam::EncoderBase::putLE<std::ostream,uint64_t>(out,(LC.alcmpstart<<16)|(LC.alstart));	

					if ( debug )
						std::cerr << "LC[" << c << "]=" << LC << std::endl;
				}
				else
				{
					libmaus::bambam::EncoderBase::putLE<std::ostream,uint64_t>(out,0);					
					if ( debug )
						std::cerr << "LC[" << c << "]=" << "null" << std::endl;
				}
			}
		}
		else
		{
			// number of distinct bins
			libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,0);
			// number of linear index chunks
			libmaus::bambam::EncoderBase::putLE<std::ostream,uint32_t>(out,0);
		}
	}

	out.flush();
	
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

				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report (default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<["+arginfo.getDefaultTmpFileName()+"]>", "temporary file prefix (default: create in current directory)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamindex(arginfo,std::cin,std::cout);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

