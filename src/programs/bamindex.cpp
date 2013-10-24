
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

#include <libmaus/aio/Buffer.hpp>
#include <libmaus/aio/SingleFileFragmentMerge.hpp>
#include <libmaus/aio/SynchronousGenericInput.hpp>
#include <libmaus/bambam/BamHeader.hpp>
#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamIndexLinearChunk.hpp>
#include <libmaus/bambam/BamIndexBinChunk.hpp>
#include <libmaus/lz/BgzfInflate.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/util/ContainerGetObject.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/util/GetObject.hpp>
#include <biobambam/Licensing.hpp>

bool getDefaultVerbose() { return true; }
bool getDefaultDisableValidation() { return false; }

struct BamIndexGenerator
{
	template<typename stream_type>
	static std::pair<int64_t,uint64_t> countDistinctBins(stream_type & stream)
	{
		int64_t prevbin = -1;
		int64_t refid = -1;
		uint64_t bins = 0;
		libmaus::bambam::BamIndexBinChunk BC;
		
		while ( stream.peek() != stream_type::traits_type::eof() )
		{
			stream.read(reinterpret_cast<char *>(&BC),sizeof(libmaus::bambam::BamIndexBinChunk));
			
			if ( refid < 0 )
				refid = BC.refid;
			
			if ( refid != static_cast<int64_t>(BC.refid) )
			{
				stream.clear();
				stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexBinChunk)),std::ios::cur);
				break;
			}
			
			if ( static_cast<int64_t>(BC.bin) != prevbin )
				bins++;
				
			prevbin = BC.bin;
		}
			
		return std::pair<uint64_t,uint64_t>(refid,bins);
	}

	template<typename stream_type>
	static std::pair<int64_t,uint64_t> countBinChunks(stream_type & stream)
	{
		int64_t bin = -1;
		int64_t refid = -1;
		uint64_t chunks = 0;

		libmaus::bambam::BamIndexBinChunk BC;
		
		while ( stream.peek() != stream_type::traits_type::eof() )
		{
			stream.read(reinterpret_cast<char *>(&BC),sizeof(libmaus::bambam::BamIndexBinChunk));

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
				stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexBinChunk)),std::ios::cur);
				stream.clear();
				break;
			}
			
			chunks++;
		}
			
		return std::pair<uint64_t,uint64_t>(bin,chunks);
	}

	template<typename stream_type>
	static std::pair<int64_t,uint64_t> countLinearChunks(stream_type & stream)
	{
		int64_t refid = -1;
		uint64_t cnt = 0;
		
		libmaus::bambam::BamIndexLinearChunk LC;
		
		while ( stream.peek() != stream_type::traits_type::eof() )
		{
			stream.read(reinterpret_cast<char *>(&LC),sizeof(libmaus::bambam::BamIndexLinearChunk));
			
			if ( refid == -1 )
				refid = LC.refid;
			
			// put back element, if it has a different refid	
			if ( static_cast<int64_t>(LC.refid) != refid )
			{
				stream.clear();
				stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexLinearChunk)),std::ios::cur);
				break;
			}
			
			cnt++;
		}
		
		return std::pair<int64_t,uint64_t>(refid,cnt);
	}

	template<typename stream_type>
	static std::pair<int64_t,uint64_t> getLinearMaxPos(stream_type & stream)
	{
		int64_t refid = -1;
		uint64_t maxpos = 0;
		
		libmaus::bambam::BamIndexLinearChunk LC;
		
		while ( stream.peek() != stream_type::traits_type::eof() )
		{
			stream.read(reinterpret_cast<char *>(&LC),sizeof(libmaus::bambam::BamIndexLinearChunk));
			
			if ( refid == -1 )
				refid = LC.refid;
			
			// put back element, if it has a different refid	
			if ( static_cast<int64_t>(LC.refid) != refid )
			{
				stream.clear();
				stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexLinearChunk)),std::ios::cur);
				break;
			}
			
			maxpos = std::max(maxpos,LC.pos);
		}
		
		return std::pair<int64_t,uint64_t>(refid,maxpos);
	}

	template<typename stream_type>
	static bool peekLinearChunk(stream_type & stream, uint64_t const refid)
	{
		libmaus::bambam::BamIndexLinearChunk LC;

		if ( stream.peek() == stream_type::traits_type::eof() )
			return false;
		
		stream.read(reinterpret_cast<char *>(&LC),sizeof(libmaus::bambam::BamIndexLinearChunk));
		stream.clear();
		stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexLinearChunk)),std::ios::cur);
		
		return LC.refid == refid;
	}

	template<typename stream_type>
	static bool peekLinearChunk(stream_type & stream, uint64_t const refid, uint64_t const pos, unsigned int const posshift)
	{
		libmaus::bambam::BamIndexLinearChunk LC;

		if ( stream.peek() == stream_type::traits_type::eof() )
			return false;
		
		stream.read(reinterpret_cast<char *>(&LC),sizeof(libmaus::bambam::BamIndexLinearChunk));
		stream.clear();
		stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexLinearChunk)),std::ios::cur);
		
		return (LC.refid == refid) && ((LC.pos >> posshift)==(pos>>posshift));
	}

	template<typename stream_type>
	static int64_t peekBin(stream_type & stream)
	{
		libmaus::bambam::BamIndexBinChunk BC;
		
		if ( stream.peek() == stream_type::traits_type::eof() )
			return -1;
			
		stream.read(
			reinterpret_cast<char *>(&BC),
			sizeof(libmaus::bambam::BamIndexBinChunk)
		);
		
		assert ( stream.gcount() == sizeof(libmaus::bambam::BamIndexBinChunk) );
		
		stream.clear();
		stream.seekg(-static_cast<int64_t>(sizeof(libmaus::bambam::BamIndexBinChunk)),std::ios::cur);
		
		return BC.refid;
	}

	static bool checkConsisteny(std::string const & binfn, std::string const & linfn, uint64_t const numrefseq)
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


	/* parser state types */
	enum parsestate { state_reading_blocklen,  state_post_skip };
	
	libmaus::bambam::BamHeader header;
	
	bool haveheader;
	
	uint64_t cacct; // accumulated compressed bytes we have read from file
	std::pair<uint64_t,uint64_t> rinfo;

        libmaus::bambam::BamHeader::BamHeaderParserState bamheaderparsestate;

	parsestate state;

	unsigned int blocklenred;
	uint32_t blocklen;

	uint64_t alcnt;
	uint64_t alcmpstart;
	uint64_t alstart;
	
	uint64_t binalcmpstart;
	uint64_t binalstart;

	libmaus::bambam::BamAlignment algn;
	uint8_t * copyptr;
	
	int64_t prevrefid;
	int64_t prevpos;
	int64_t prevbin;
	int64_t prevcheckrefid;
	int64_t prevcheckpos;
	
	std::string const binchunktmpfilename;
	std::string const linchunktmpfilename;
	libmaus::aio::CheckedOutputStream::unique_ptr_type chunkCOS;
	libmaus::aio::CheckedOutputStream::unique_ptr_type linCOS;

	libmaus::aio::Buffer<libmaus::bambam::BamIndexBinChunk> BCB;
	std::vector< std::pair<uint64_t,uint64_t> > binchunkfrags;
	libmaus::aio::Buffer<libmaus::bambam::BamIndexLinearChunk> BLC;
	std::vector< std::pair<uint64_t,uint64_t> > linearchunkfrags;
	
	unsigned int const verbose;
	bool const validate;
	bool const debug;
	
	bool blocksetup;

	BamIndexGenerator(
		std::string const & tmpfileprefix, 
		unsigned int const rverbose = 0,
		bool const rvalidate = true, bool const rdebug = false
	)
	: header(), haveheader(false),
	  cacct(0), rinfo(), bamheaderparsestate(), state(state_reading_blocklen), blocklenred(0),
	  blocklen(0), alcnt(0), alcmpstart(0), alstart(0), binalcmpstart(0), binalstart(0), algn(),
	  copyptr(0), prevrefid(-1), prevpos(-1), prevbin(-1), prevcheckrefid(-1), prevcheckpos(-1),
	  binchunktmpfilename(tmpfileprefix+".bin"), linchunktmpfilename(tmpfileprefix+".lin"),
	  chunkCOS(), linCOS(),
	  BCB(), binchunkfrags(), BLC(), linearchunkfrags(),
	  verbose(rverbose), validate(rvalidate), debug(rdebug),
	  blocksetup(false)
	{
		libmaus::util::TempFileRemovalContainer::addTempFile(binchunktmpfilename);
		libmaus::util::TempFileRemovalContainer::addTempFile(linchunktmpfilename);

		libmaus::aio::CheckedOutputStream::unique_ptr_type tchunkCOS(new libmaus::aio::CheckedOutputStream(binchunktmpfilename));
		libmaus::aio::CheckedOutputStream::unique_ptr_type tlinCOS(new libmaus::aio::CheckedOutputStream(linchunktmpfilename));
		
		chunkCOS = UNIQUE_PTR_MOVE(tchunkCOS);
		linCOS = UNIQUE_PTR_MOVE(tlinCOS);
	}
	
	void addBlock(
		uint8_t const * Bbegin,
		uint64_t const compsize,
		uint64_t const uncompsize
	)
	{
		uint8_t const * pa = Bbegin; // buffer current pointer
		uint8_t const * pc = pa + uncompsize; // buffer end pointer
		
		rinfo.first = compsize;
		rinfo.second = uncompsize;

		cacct += rinfo.first; // accumulator for compressed data bytes read

		if ( (! haveheader) && (pa != pc) )
		{			
			libmaus::util::GetObject<uint8_t const *> G(Bbegin);
			std::pair<bool,uint64_t> const P = libmaus::bambam::BamHeader::parseHeader(G,bamheaderparsestate,uncompsize);

			// header complete?
			if ( P.first )
			{
				header.init(bamheaderparsestate);
				haveheader = true;
				pa = Bbegin + P.second;
				pc = Bbegin + uncompsize;
				
				if ( pa != pc )
				{
					// start of compressed block
					alcmpstart = cacct - rinfo.first;
					// start of alignment in uncompressed block
					alstart = pa - Bbegin;
		
					// start of bin in compressed data
					binalcmpstart = alcmpstart;
					// start of bin in uncompressed block
					binalstart = alstart;
			
					blocksetup = true;	
				}
			}
		}
		
		if ( (haveheader) && (pa != pc) )
		{
			if ( ! blocksetup )
			{
				// start of compressed block
				alcmpstart = cacct - rinfo.first;
				// start of alignment in uncompressed block
				alstart = pa - Bbegin;
		
				// start of bin in compressed data
				binalcmpstart = alcmpstart;
				// start of bin in uncompressed block
				binalstart = alstart;
			
				blocksetup = true;
			}

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
							if ( validate )
								algn.checkAlignment();
						
							int64_t const thisrefid = algn.getRefID();
							int64_t const thispos = algn.getPos();
							int64_t const thisbin = 
								(thisrefid >= 0 && thispos >= 0) ? algn.computeBin() : -1;

							int64_t const thischeckrefid = (thisrefid >= 0) ? thisrefid : std::numeric_limits<int64_t>::max();
							int64_t const thischeckpos   = (thispos >= 0) ? thispos : std::numeric_limits<int64_t>::max();
							
							bool const orderok =
								(thischeckrefid > prevcheckrefid)
								||
								(thischeckrefid == prevcheckrefid && thischeckpos >= prevcheckpos);
								
							if ( ! orderok )
							{
								libmaus::exception::LibMausException se;
								se.getStream() << "File is not sorted by coordinate." << std::endl;
								se.finish();
								throw se;
							}
								
							if ( 
								thisrefid != prevrefid ||
								(
									(thisrefid == prevrefid) && (thispos>>14)!=(prevpos>>14)
								)
							)
							{
								if ( thisrefid >= 0 && thispos >= 0 )
								{
									libmaus::bambam::BamIndexLinearChunk LC(thisrefid,thispos,alcmpstart,alstart);
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

									libmaus::bambam::BamIndexBinChunk BC(prevrefid,prevbin,binalcmpstart,binalstart,binalcmpend,binalend);
									
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
							uint64_t const nextalstart = pa - Bbegin;

							alcnt++;
							alcmpstart = nextalcmpstart;
							alstart = nextalstart;

							prevrefid = thisrefid;
							prevpos = thispos;
							prevbin = thisbin;

							prevcheckrefid = thischeckrefid;
							prevcheckpos = thischeckpos;
		       
							if ( verbose && (alcnt % (1024*1024) == 0) )
								std::cerr << "[V] " << alcnt/(1024*1024) << std::endl;
						}
						break;
					}
				}
			}
		}
	}
	
	template<typename output_stream_type>
	void flush(output_stream_type & out)
	{
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

			libmaus::bambam::BamIndexBinChunk BC(prevrefid,prevbin,binalcmpstart,binalstart,binalcmpend,binalend);

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

		libmaus::aio::SingleFileFragmentMerge<libmaus::bambam::BamIndexBinChunk>::merge(binchunktmpfilename,binchunkfrags);
		libmaus::aio::SingleFileFragmentMerge<libmaus::bambam::BamIndexLinearChunk>::merge(linchunktmpfilename,linearchunkfrags);

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
				std::pair<int64_t,uint64_t> const P = countDistinctBins(binCIS);
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
						libmaus::bambam::BamIndexBinChunk BC;
						binCIS.read(reinterpret_cast<char *>(&BC),sizeof(libmaus::bambam::BamIndexBinChunk));
						
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
						libmaus::bambam::BamIndexLinearChunk LC;
						linCIS.read(reinterpret_cast<char *>(&LC),sizeof(libmaus::bambam::BamIndexLinearChunk));
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
	}
};

int bamindex(libmaus::util::ArgInfo const & arginfo, std::istream & in, std::ostream & out)
{
	bool const debug = arginfo.getValue<unsigned int>("debug",0);
	unsigned int const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	bool const validate = !(arginfo.getValue<unsigned int>("disablevalidation",getDefaultDisableValidation()));
	std::string const tmpfileprefix = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());

	libmaus::lz::BgzfInflate<std::istream> rec(in);
	
	BamIndexGenerator BIG(tmpfileprefix,verbose,validate,debug);

	libmaus::autoarray::AutoArray<uint8_t> B(libmaus::lz::BgzfConstants::getBgzfMaxBlockSize());
	std::pair<uint64_t,uint64_t> rinfo;
	while ((rinfo=rec.readPlusInfo(reinterpret_cast<char *>(B.begin()),B.size())).second)
		BIG.addBlock(B.begin(),rinfo.first,rinfo.second);

	BIG.flush(out);
	
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
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable alignment validation (default: 0)" ) );
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

