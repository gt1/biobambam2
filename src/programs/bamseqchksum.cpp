/**
    bambam
    Copyright (C) 2014-2014 David K. Jackson
    Copyright (C) 2009-2016 German Tischler
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
#include <config.h>

#include <iostream>

#include <libmaus2/bambam/ChecksumsFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/util/ArgInfo.hpp>

#include <biobambam2/Licensing.hpp>
#include <biobambam2/BamBamConfig.hpp>

#if defined(BIOBAMBAM_HAVE_GMP)
#include <libmaus2/math/UnsignedInteger.hpp>
#include <gmp.h>
#endif

static int getDefaultVerbose() { return 0; }
static std::string getDefaultInputFormat() { return "bam"; }
static std::string getDefaultHash() { return "crc32prod"; }
static int getDefaultStats() { return 0; }
static int getDefaultCheckHamming() { return 0; }

struct DComparatorBase
{
	uint64_t const d;
	std::vector<uint8_t> const & VV;

	DComparatorBase(uint64_t const rd, std::vector<uint8_t> const & rVV) : d(rd), VV(rVV) {}

	uint64_t hamming(uint64_t const i, uint64_t const j) const
	{
		uint64_t p0 = i*d;
		uint64_t p1 = j*d;
		uint64_t h = 0;

		for ( uint64_t k = 0; k < d; ++k )
		{
			uint64_t v0 = VV[p0+k];
			uint64_t v1 = VV[p1+k];
			uint64_t v2 = v0^v1;
			h += libmaus2::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(v2);
		}
		return h;
	}

	uint64_t getBits(uint64_t const i, uint64_t from, uint64_t len) const
	{
		assert ( len <= 64 );

		unsigned int const rlen = len;

		uint64_t p0 = i*d;
		uint64_t const byteskip = from / 8;

		if ( byteskip )
		{
			p0 += byteskip;
			from -= byteskip * 8;
		}

		assert ( from < 8 );

		uint64_t w = 0;
		while ( len )
		{
			uint64_t v = VV[p0++];
			// bits available in this word/byte
			unsigned int const wlen = std::min(static_cast<uint64_t>(8-from),len);

			// align on left
			v <<= from;
			v &= 0xFF;
			// align on right
			v >>= (8-wlen);

			w <<= wlen;
			w |= v;

			from = 0;
			len -= wlen;
		}

		assert ( (rlen == 64) || w < (1ull<<rlen) );

		return w;
	}

	void print(std::ostream & out, uint64_t const i) const
	{
		uint64_t p0 = i*d;
		for ( uint64_t k = 0; k < d; ++k )
			out << std::hex << std::setw(2) << std::setfill('0') << static_cast<unsigned int>(VV[p0+k]) << std::setw(0) << std::dec;
	}
};

struct DComparator : public DComparatorBase
{
	DComparator(uint64_t const rd, std::vector<uint8_t> const & rVV) : DComparatorBase(rd,rVV) {}

	bool compare(uint64_t r0, uint64_t r1) const
	{
		uint64_t const p0 = r0*d;
		uint64_t const p1 = r1*d;

		for ( uint64_t k = 0; k < d; ++k )
		{
			uint8_t const v0 = VV[p0+k];
			uint8_t const v1 = VV[p1+k];

			if ( v0 < v1 )
				return true;
			else if ( v1 < v0 )
				return false;
		}

		return false;
	}

	uint64_t find(uint64_t rref, std::vector<uint64_t> const & I) const
	{
		uint64_t low = 0, high = I.size();

		while ( high-low > 1 )
		{
			uint64_t mid = (high+low)/2;

			if ( compare(rref,I[mid]) ) // rref < I[mid]
				high = mid;
			else // rref >= I[mid]
				low = mid;
		}

		return low;
	}

	bool operator()(uint64_t const i, uint64_t const j) const
	{
		return compare(i,j);
	}

	bool eq(uint64_t const i, uint64_t const j) const
	{
		return (!(*this)(i,j)) && (!(*this)(j,i));
	}
};

struct DBitComparator : public DComparatorBase
{
	unsigned int const from;
	unsigned int const len;

	DBitComparator(uint64_t const rd, std::vector<uint8_t> const & rVV, unsigned int const rfrom, unsigned int const rlen) : DComparatorBase(rd,rVV), from(rfrom), len(rlen) {}

	uint64_t get(uint64_t const i) const
	{
		return getBits(i,from,len);
	}

	bool operator()(uint64_t const i, uint64_t const j) const
	{
		uint64_t const v0 = get(i);
		uint64_t const v1 = get(j);
		return v0 < v1;
	}
};

void evaluateStats(std::vector<uint8_t> const & TV, std::vector<uint8_t> const & VV, bool checkHamming)
{
	uint64_t const d = TV.size();
	bool collisionfound = false;

	if ( d )
	{
		assert ( VV.size() % d == 0 );
		uint64_t const n = VV.size() / d;

		// bit counters
		std::vector<uint64_t> CC0(8*d,0);

		for ( uint64_t i = 0; i < n; ++i )
		{
			uint64_t o0 = i*d;
			uint64_t o = 0;

			for ( uint64_t j = 0; j < d; ++j )
			{
				uint8_t const v0 = VV.at(o0++);

				for ( unsigned int k = 0; k < 8; ++k )
				{
					uint64_t const b = (v0 >> (8-1-k)) & 1;
					CC0.at(o++) += b;
				}
			}
		}

		double maxabsbad = 0.0;
		double avgabsbad = 0.0;
		uint64_t absused = 0;
		for ( uint64_t i = 0; i < CC0.size(); ++i )
			if ( CC0[i] )
			{
				uint64_t const v = CC0.at(i);
				double const dv = v;
				double const frac = dv / n;
				double const bad = std::abs(0.5-frac);
				// std::cerr << "[S] abs msb " << i << " bad " << bad << std::endl;
				maxabsbad = std::max(bad,maxabsbad);
				avgabsbad += bad;
				absused += 1;
			}
		avgabsbad = absused ? avgabsbad/absused : avgabsbad;
		std::cerr << "[S] abs max bad " << maxabsbad << " avg bad " << avgabsbad << std::endl;

		// bit difference counters
		std::vector<uint64_t> CC(8*d,0);

		for ( uint64_t i = 1; i < n; ++i )
		{
			uint64_t o0 = (i-1)*d;
			uint64_t o1 = (i*d);
			uint64_t o = 0;

			for ( uint64_t j = 0; j < d; ++j )
			{
				uint8_t const v0 = VV.at(o0++);
				uint8_t const v1 = VV.at(o1++);
				uint8_t const v2 = v0 ^ v1;

				for ( unsigned int k = 0; k < 8; ++k )
				{
					uint64_t const b = (v2 >> (8-1-k)) & 1;
					CC.at(o++) += b;
				}
			}
		}

		if ( n > 1 )
		{
			double maxbad = 0;
			double avgbad = 0;

			uint64_t used = 0;
			for ( uint64_t i = 0; i < CC.size(); ++i )
				if ( CC[i] )
				{
					double const bad = std::abs(0.5 - static_cast<double>(CC[i]) / (n-1));
					maxbad = std::max(maxbad,bad);
					avgbad += bad;
					used += 1;
					//std::cerr << "MSB " << i << " prob " << static_cast<double>(CC[i]) / (n-1) << std::endl;
				}
			std::cerr << "[S] dif max bad " << maxbad << " avg bad " << (used?(avgbad/used) : 0) << std::endl;
		}

		// sort code words
		std::vector<uint64_t> I(n);
		for ( uint64_t i = 0; i < n; ++i )
			I[i] = i;
		DComparator const comp(d,VV);
		std::sort(I.begin(),I.end(),comp);

		// check for collisions
		uint64_t low = 0;
		while ( low < I.size() )
		{
			uint64_t high = low+1;
			while ( high < I.size() && comp.eq(I[low],I[high]) )
				++high;

			if ( high-low > 1 )
			{
				if ( ! collisionfound )
					std::cerr << "[S] collision !" << std::endl;
				collisionfound = true;
			}

			low = high;
		}

		if ( checkHamming )
		{
			// minimum Hamming distance found
			std::vector<uint64_t> M(n,std::numeric_limits<uint64_t>::max());
			// key bit length
			unsigned int const bitlength = d*8;

			// maximum number of errors (Hamming distance) checked for
			unsigned int const maxerr = 7;
			unsigned int const tparts = std::min(maxerr+1,bitlength);
			assert ( tparts <= bitlength );
			unsigned int const tbitsperpart = (bitlength/tparts);
			assert ( tbitsperpart != 0 );
			unsigned int const numparts = (bitlength + tbitsperpart -1)/tbitsperpart;
			unsigned int const bitsperpart = (bitlength + numparts-1)/numparts;

			unsigned int const rest = bitlength % bitsperpart;
			unsigned int const frac = rest ? (bitsperpart-rest) : 0;
			unsigned int const full = numparts - frac;

			for ( unsigned int p = 0; p < numparts; ++p )
			{
				bool const isfull = p < full;
				uint64_t const b_low = isfull ? p*bitsperpart : full*bitsperpart + (p-full)*(bitsperpart-1);
				uint64_t const b_size = isfull ? bitsperpart : (bitsperpart-1);

				DBitComparator const bcomp(d,VV,b_low,b_size);

				std::cerr << "[S] part " << (p+1) << "/" << numparts << " b_low=" << b_low << " b_size=" << b_size << std::endl;

				std::vector<uint64_t> J(n);
				for ( uint64_t i = 0; i < n; ++i )
					J[i] = i;
				std::sort(J.begin(),J.end(),bcomp);

				uint64_t low = 0;
				while ( low < J.size() )
				{
					uint64_t high = low+1;
					while ( high < J.size() && bcomp.get(J[high]) == bcomp.get(J[low]) )
						++high;

					for ( uint64_t i0 = low; i0 < high; ++i0 )
						for ( uint64_t i1 = low; i1 < high; ++i1 )
							if ( i0 != i1 )
							{
								uint64_t const h = bcomp.hamming(J[i0],J[i1]);

								uint64_t const rr = comp.find(J[i0],I); // find rank on I

								// check we found the correct word
								for ( uint64_t zz = 0; zz < comp.d; ++zz )
									assert ( VV [ d * J[i0] + zz ] == VV [ d * I[rr] + zz ]);

								// update minimum distance for code word
								M[rr] = std::min(M[rr],h);
							}

					low = high;
				}
			}

			uint64_t minh = std::numeric_limits<uint64_t>::max();
			for ( uint64_t i = 0; i < M.size(); ++i )
				minh = std::min(minh,M[i]);
			if ( minh <= maxerr )
				std::cerr << "[S] minimum distance " << minh << std::endl;
			else
				std::cerr << "[S] no codeword pair within Hamming distance " << maxerr << std::endl;
		}
	}
}

int bamseqchksumStats(::libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing to read data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	libmaus2::timing::RealTimeClock rtc;
	rtc.start();
	double prevtime = 0;

	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const checkhamming = arginfo.getValue<int>("checkhamming",getDefaultCheckHamming());
	std::string const hash = arginfo.getValue<std::string>("hash",getDefaultHash());

	// input decoder wrapper
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus2::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();

	::libmaus2::bambam::BamAlignment & algn = dec.getAlignment();

	libmaus2::bambam::ChecksumsInterface::unique_ptr_type Pchksums(libmaus2::bambam::ChecksumsFactory::construct(hash,header));
	libmaus2::bambam::ChecksumsInterface & chksums = *Pchksums;

	std::vector<uint8_t> TV;
	std::vector<uint8_t> b_seq_V;

	uint64_t c = 0;
	while ( dec.readAlignment() )
	{
		chksums.update(algn);

		chksums.get_b_seq_all(TV);
		for ( uint64_t i = 0; i < TV.size(); ++i )
			b_seq_V.push_back(TV[i]);

		if ( verbose && (++c & (1024*1024-1)) == 0 )
		{
			double const elapsed = rtc.getElapsedSeconds();
			chksums.printVerbose(std::cerr,c,algn,elapsed-prevtime);
			prevtime = elapsed;

			evaluateStats(TV,b_seq_V,checkhamming);
		}
	}

	chksums.printChecksums(std::cout);

	if ( verbose )
	{
		std::cerr << "[V] run time " << rtc.getElapsedSeconds() << " (" << rtc.formatTime(rtc.getElapsedSeconds()) << ")" << std::endl;
	}

	if ( TV.size() )
	{
		uint64_t const num = b_seq_V.size() / TV.size();
		std::cerr << "[S] num " << num << std::endl;
		evaluateStats(TV,b_seq_V,checkhamming);
	}

	return EXIT_SUCCESS;
}

int bamseqchksum(::libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing to read data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	libmaus2::timing::RealTimeClock rtc;
	rtc.start();
	double prevtime = 0;

	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	std::string const hash = arginfo.getValue<std::string>("hash",getDefaultHash());

	// input decoder wrapper
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus2::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();

	::libmaus2::bambam::BamAlignment & algn = dec.getAlignment();

	libmaus2::bambam::ChecksumsInterface::unique_ptr_type Pchksums(libmaus2::bambam::ChecksumsFactory::construct(hash,header));
	libmaus2::bambam::ChecksumsInterface & chksums = *Pchksums;

	uint64_t c = 0;
	while ( dec.readAlignment() )
	{
		chksums.update(algn);

		if ( verbose && (++c & (1024*1024-1)) == 0 )
		{
			double const elapsed = rtc.getElapsedSeconds();
			chksums.printVerbose(std::cerr,c,algn,elapsed-prevtime);
			prevtime = elapsed;
		}
	}

	chksums.printChecksums(std::cout);

	if ( verbose )
	{
		std::cerr << "[V] run time " << rtc.getElapsedSeconds() << " (" << rtc.formatTime(rtc.getElapsedSeconds()) << ")" << std::endl;
	}

	return EXIT_SUCCESS;
}

#if defined(BIOBAMBAM_HAVE_GMP)
template<size_t k>
static libmaus2::math::UnsignedInteger<k> convertNumber(mpz_t const & gmpnum)
{
	size_t const numbitsperel = 8 * sizeof(uint32_t);
	size_t const numwords = (mpz_sizeinbase(gmpnum,2) + numbitsperel - 1) / numbitsperel;
	libmaus2::autoarray::AutoArray<uint32_t> A(numwords,false);
	size_t countp = numwords;
	mpz_export(A.begin(),&countp,-1,sizeof(uint32_t),0,0,gmpnum);
	libmaus2::math::UnsignedInteger<k> U;
	for ( size_t i = 0; i < std::min(countp,static_cast<size_t>(k)); ++i )
		U[i] = A[i];
	return U;
}
// search for next prime number larger than 2^(32*k) using probabilistic algorithm in gmp library
template<size_t k>
libmaus2::math::UnsignedInteger<k+1> nextPrime()
{
	libmaus2::math::UnsignedInteger<k+1> U(1);
	U <<= (k*32);

	mpz_t gmpU;
	mpz_init(gmpU);
	mpz_import(gmpU,k+1,-1 /* least sign first */,sizeof(uint32_t),0 /* native endianess */,0 /* don't skip bits */, U.getWords());

	while ( mpz_probab_prime_p(gmpU,200) == 0 )
		mpz_add_ui(gmpU,gmpU,1);

	libmaus2::math::UnsignedInteger<k+1> const R = convertNumber<k+1>(gmpU);

	mpz_clear(gmpU);

	return R;
}
#endif

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

				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );
				#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif

				V.push_back ( std::pair<std::string,std::string> ( std::string("hash=<[")+getDefaultHash()+"]>", "hash digest function: " + libmaus2::bambam::ChecksumsFactory::getSupportedHashVariantsList()) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;

				return EXIT_SUCCESS;
			}

		if ( arginfo.getValue<int>("stats",getDefaultStats()) )
			return bamseqchksumStats(arginfo);
		else
			return bamseqchksum(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
