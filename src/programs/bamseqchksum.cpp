/**
    bambam
    Copyright (C) 2014-2014 David K. Jackson
    Copyright (C) 2009-2013 German Tischler
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
#include "config.h"

#include <iostream>

#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/digest/Digests.hpp>
#include <libmaus/util/ArgInfo.hpp>

#include <libmaus/util/I386CacheLineSize.hpp>

#include <biobambam/Licensing.hpp>
#include <biobambam/BamBamConfig.hpp>

#if defined(BIOBAMBAM_HAVE_GMP)
#include <gmp.h>
#endif

static int getDefaultVerbose() { return 0; }
static std::string getDefaultInputFormat()
{
	return "bam";
}
static std::vector<std::string> getDefaultAuxTags() {
	static const char * defaults[] = {"BC","FI","QT","RT","TC"}; //lexical order
	std::vector<std::string> v (defaults, defaults + sizeof(defaults)/sizeof(*defaults));
	return v;
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getMersenne31()
{
	return libmaus::math::UnsignedInteger<k>(0x7FFFFFFFull);
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getMersenne521()
{
	libmaus::math::UnsignedInteger< (521 + 31)/32 > M;

	M[0] = 1;
	M <<= 521;
	M -= 1;
	
	return libmaus::math::UnsignedInteger<k>(M);
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime64()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 59UL;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	return M; // 2^64-59 prime according to https://primes.utm.edu/lists/2small/0bit.html
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime96()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 17;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	if ( 2 < k )
		M[2] = 0xFFFFFFFFUL;
	return M; // 2^96-17 prime according to https://primes.utm.edu/lists/2small/0bit.html
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime128()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 159;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	if ( 2 < k )
		M[2] = 0xFFFFFFFFUL;
	if ( 3 < k )
		M[3] = 0xFFFFFFFFUL;
	return M; // 2^128-159 prime according to https://primes.utm.edu/lists/2small/100bit.html
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime160()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 47;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	if ( 2 < k )
		M[2] = 0xFFFFFFFFUL;
	if ( 3 < k )
		M[3] = 0xFFFFFFFFUL;
	if ( 4 < k )
		M[4] = 0xFFFFFFFFUL;
	return M; // 2^160-47 prime according to https://primes.utm.edu/lists/2small/100bit.html
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime192()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 237;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	if ( 2 < k )
		M[2] = 0xFFFFFFFFUL;
	if ( 3 < k )
		M[3] = 0xFFFFFFFFUL;
	if ( 4 < k )
		M[4] = 0xFFFFFFFFUL;
	if ( 5 < k )
		M[5] = 0xFFFFFFFFUL;
	return M; // 2^192-237 prime according to https://primes.utm.edu/lists/2small/100bit.html
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime224()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 63;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	if ( 2 < k )
		M[2] = 0xFFFFFFFFUL;
	if ( 3 < k )
		M[3] = 0xFFFFFFFFUL;
	if ( 4 < k )
		M[4] = 0xFFFFFFFFUL;
	if ( 5 < k )
		M[5] = 0xFFFFFFFFUL;
	if ( 6 < k )
		M[6] = 0xFFFFFFFFUL;
	return M; // 2^224-63 prime according to https://primes.utm.edu/lists/2small/200bit.html
}

template<size_t k>
static libmaus::math::UnsignedInteger<k> getPrime256()
{
	libmaus::math::UnsignedInteger<k> M;
	if ( 0 < k )
		M[0] = 0xFFFFFFFFUL - 189;
	if ( 1 < k )
		M[1] = 0xFFFFFFFFUL;
	if ( 2 < k )
		M[2] = 0xFFFFFFFFUL;
	if ( 3 < k )
		M[3] = 0xFFFFFFFFUL;
	if ( 4 < k )
		M[4] = 0xFFFFFFFFUL;
	if ( 5 < k )
		M[5] = 0xFFFFFFFFUL;
	if ( 6 < k )
		M[6] = 0xFFFFFFFFUL;
	if ( 7 < k )
		M[7] = 0xFFFFFFFFUL;
	return M; // 2^256-189 prime according to https://primes.utm.edu/lists/2small/200bit.html
}

template<typename _digest_type>
struct UpdateContext
{
	typedef _digest_type digest_type;

	digest_type ctx_name_flags_seq;
	uint8_t ctx_name_flags_seq_digest[digest_type::digestlength];
	libmaus::math::UnsignedInteger<digest_type::digestlength/4> name_flags_seq_digest;
	
	digest_type ctx_flags_seq;
	uint8_t ctx_flags_seq_digest[digest_type::digestlength];
	libmaus::math::UnsignedInteger<digest_type::digestlength/4> flags_seq_digest;

	digest_type ctx_flags_seq_qual;
	uint8_t ctx_flags_seq_qual_digest[digest_type::digestlength];
	libmaus::math::UnsignedInteger<digest_type::digestlength/4> flags_seq_qual_digest;
	
	digest_type ctx_flags_seq_tags;
	uint8_t ctx_flags_seq_tags_digest[digest_type::digestlength];
	libmaus::math::UnsignedInteger<digest_type::digestlength/4> flags_seq_tags_digest;
	
	bool pass;
	bool valid;
};

typedef UpdateContext<libmaus::digest::Null> NullUpdateContext;
typedef UpdateContext<libmaus::digest::CRC32> CRC32UpdateContext;
typedef UpdateContext<libmaus::util::MD5> MD5UpdateContext;
typedef UpdateContext<libmaus::digest::SHA1> SHA1UpdateContext;
typedef UpdateContext<libmaus::digest::SHA2_224> SHA2_224_UpdateContext;
typedef UpdateContext<libmaus::digest::SHA2_256> SHA2_256_UpdateContext;
typedef UpdateContext<libmaus::digest::SHA2_256_sse4> SHA2_256_sse4_UpdateContext;
typedef UpdateContext<libmaus::digest::SHA2_384> SHA2_384_UpdateContext;
typedef UpdateContext<libmaus::digest::SHA2_512> SHA2_512_UpdateContext;
typedef UpdateContext<libmaus::digest::SHA2_512_sse4> SHA2_512_sse4_UpdateContext;

template<typename _context_type>
struct SimpleSums
{
	typedef _context_type context_type;
	typedef SimpleSums<context_type> this_type;
	
	private:
	uint64_t count;
	libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> b_seq;
	libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> name_b_seq;
	libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> b_seq_qual;
	libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> b_seq_tags;
	
	public:
	void push(context_type const & context)
	{
		count += 1;
		b_seq += context.flags_seq_digest;
		name_b_seq += context.name_flags_seq_digest;
		b_seq_qual += context.flags_seq_qual_digest;
		b_seq_tags += context.flags_seq_tags_digest;
	}

	void push (this_type const & subsetproducts)
	{
		count += subsetproducts.count;
		b_seq += subsetproducts.b_seq;
		name_b_seq += subsetproducts.name_b_seq;
		b_seq_qual += subsetproducts.b_seq_qual;
		b_seq_tags += subsetproducts.b_seq_tags;
	};
	
	bool operator== (this_type const & other) const
	{
		return count==other.count &&
			b_seq==other.b_seq && name_b_seq==other.name_b_seq &&
			b_seq_qual==other.b_seq_qual && b_seq_tags==other.b_seq_tags;
	};

	std::string get_b_seq() const
	{
		std::ostringstream ostr;
		ostr << std::hex << b_seq;
		return ostr.str();
	}

	std::string get_name_b_seq() const
	{
		std::ostringstream ostr;
		ostr << std::hex << name_b_seq;
		return ostr.str();
	}
	
	std::string get_b_seq_qual() const
	{
		std::ostringstream ostr;
		ostr << std::hex << b_seq_qual;
		return ostr.str();
	}
	
	std::string get_b_seq_tags() const
	{
		std::ostringstream ostr;
		ostr << std::hex << b_seq_tags;
		return ostr.str();
	}
	
	uint64_t get_count() const
	{
		return count;
	}
};

typedef SimpleSums<CRC32UpdateContext> CRC32SimpleSums;
typedef SimpleSums<MD5UpdateContext> MD5SimpleSums;
typedef SimpleSums<SHA1UpdateContext> SHA1SimpleSums;
typedef SimpleSums<SHA2_224_UpdateContext> SHA2_224_SimpleSums;
typedef SimpleSums<SHA2_256_UpdateContext> SHA2_256_SimpleSums;
typedef SimpleSums<SHA2_256_sse4_UpdateContext> SHA2_256_sse4_SimpleSums;
typedef SimpleSums<SHA2_384_UpdateContext> SHA2_384_SimpleSums;
typedef SimpleSums<SHA2_512_UpdateContext> SHA2_512_SimpleSums;
typedef SimpleSums<SHA2_512_sse4_UpdateContext> SHA2_512_sse4_SimpleSums;

template<size_t a, size_t b>
struct TemplateMax
{
	static size_t const value = a > b ? a : b;
};

template<typename _context_type, size_t _primeWidth>
struct PrimeProduct
{
	typedef _context_type context_type;
	typedef PrimeProduct<context_type,_primeWidth> this_type;
	
	// width of the prime in bits
	static size_t const primeWidth = _primeWidth;
	// width of the products in bits
	static size_t const productWidth = 2 * TemplateMax<primeWidth/32, context_type::digest_type::digestlength/4>::value;
	// folding prime
	static libmaus::math::UnsignedInteger<productWidth> const foldprime;
	// prime
	static libmaus::math::UnsignedInteger<productWidth> const prime;
		
	private:
	uint64_t count;

	#if defined(BIOBAMBAM_HAVE_GMP)
	mpz_t gmpprime;
	mpz_t gmpfoldprime;	
	mpz_t gmpb_seq;
	mpz_t gmpname_b_seq;
	mpz_t gmpb_seq_qual;
	mpz_t gmpb_seq_tags;
	mpz_t gmpcrc;
	mpz_t gmpfoldcrc;
	mpz_t gmpprod;
	#else
	libmaus::math::UnsignedInteger<productWidth> b_seq;
	libmaus::math::UnsignedInteger<productWidth> name_b_seq;
	libmaus::math::UnsignedInteger<productWidth> b_seq_qual;
	libmaus::math::UnsignedInteger<productWidth> b_seq_tags;	
	#endif
	
	template<size_t k>
	void updateProduct(libmaus::math::UnsignedInteger<productWidth> & product, libmaus::math::UnsignedInteger<k> const & rcrc)
	{
		// convert crc to different size
		libmaus::math::UnsignedInteger<productWidth> crc(rcrc);
		// break down crc to range < prime
		crc = crc % foldprime;
		// make sure it is not null
		if ( crc.isNull() )
			crc = libmaus::math::UnsignedInteger<productWidth>(1);
		product = product * crc;
		product = product % prime;
	}

	#if defined(BIOBAMBAM_HAVE_GMP)
	template<size_t k>
	void updateProduct(mpz_t & gmpproduct, libmaus::math::UnsignedInteger<k> const & rcrc)
	{
		mpz_import(gmpcrc,k,-1 /* least sign first */,sizeof(uint32_t),0 /* native endianess */,0 /* don't skip bits */, rcrc.getWords());
		mpz_mod(gmpfoldcrc,gmpcrc,gmpfoldprime); // foldcrc = crc % foldprime
		if ( mpz_cmp_ui(gmpfoldcrc,0) == 0 )
			mpz_set_ui(gmpfoldcrc,1);
		mpz_mul(gmpprod,gmpproduct,gmpfoldcrc); // prod = product * foldcrc
		mpz_mod(gmpproduct,gmpprod,gmpprime); // product = prod % prime
	}
	#endif
	
	#if defined(BIOBAMBAM_HAVE_GMP)
	static libmaus::math::UnsignedInteger<productWidth> convertNumber(mpz_t const & gmpnum)
	{
		size_t const numbitsperel = 8 * sizeof(uint32_t);
		size_t const numwords = (mpz_sizeinbase(gmpnum,2) + numbitsperel - 1) / numbitsperel;
		libmaus::autoarray::AutoArray<uint32_t> A(numwords,false);
		size_t countp = numwords;
		mpz_export(A.begin(),&countp,-1,sizeof(uint32_t),0,0,gmpnum);
		libmaus::math::UnsignedInteger<productWidth> U;
		for ( size_t i = 0; i < std::min(countp,static_cast<size_t>(productWidth)); ++i )
			U[i] = A[i];
		return U;
	}
	#endif
	
	public:
	PrimeProduct()
	: count(0)
		#if !defined(BIOBAMBAM_HAVE_GMP)
		, b_seq(1), name_b_seq(1), b_seq_qual(1), b_seq_tags(1) 
		#endif
	{
		#if defined(BIOBAMBAM_HAVE_GMP)
		mpz_init(gmpprime);
		mpz_import(gmpprime,productWidth,-1 /* least sign first */,sizeof(uint32_t),0 /* native endianess */,0 /* don't skip bits */, prime.getWords());		
		mpz_init(gmpfoldprime);
		mpz_import(gmpfoldprime,productWidth,-1 /* least sign first */,sizeof(uint32_t),0 /* native endianess */,0 /* don't skip bits */, foldprime.getWords());
		mpz_init_set_ui(gmpb_seq,1);
		mpz_init_set_ui(gmpname_b_seq,1);
		mpz_init_set_ui(gmpb_seq_qual,1);
		mpz_init_set_ui(gmpb_seq_tags,1);
		mpz_init(gmpcrc);
		mpz_init(gmpfoldcrc);
		mpz_init(gmpprod);
		#endif
	}
	~PrimeProduct()
	{
		#if defined(BIOBAMBAM_HAVE_GMP)
		mpz_clear(gmpprime);
		mpz_clear(gmpfoldprime);
		mpz_clear(gmpb_seq);
		mpz_clear(gmpname_b_seq);
		mpz_clear(gmpb_seq_qual);
		mpz_clear(gmpb_seq_tags);
		mpz_clear(gmpfoldcrc);
		mpz_clear(gmpprod);
		#endif		
	}
	
	void push(context_type const & context)
	{
		count += 1;		
		#if defined(BIOBAMBAM_HAVE_GMP)
		updateProduct(gmpb_seq,context.flags_seq_digest);
		updateProduct(gmpname_b_seq,context.name_flags_seq_digest);
		updateProduct(gmpb_seq_qual,context.flags_seq_qual_digest);
		updateProduct(gmpb_seq_tags,context.flags_seq_tags_digest);
		#else
		updateProduct(b_seq,context.flags_seq_digest);
		updateProduct(name_b_seq,context.name_flags_seq_digest);
		updateProduct(b_seq_qual,context.flags_seq_qual_digest);
		updateProduct(b_seq_tags,context.flags_seq_tags_digest);
		#endif
	}

	void push (this_type const & subsetproducts)
	{
		count += subsetproducts.count;
		#if defined(BIOBAMBAM_HAVE_GMP)
		updateProduct(gmpb_seq,convertNumber(subsetproducts.gmpb_seq));
		updateProduct(gmpname_b_seq,convertNumber(subsetproducts.gmpname_b_seq));
		updateProduct(gmpb_seq_qual,convertNumber(subsetproducts.gmpb_seq_qual));
		updateProduct(gmpb_seq_tags,convertNumber(subsetproducts.gmpb_seq_tags));
		#else
		updateProduct(b_seq,subsetproducts.b_seq);
		updateProduct(name_b_seq,subsetproducts.name_b_seq);
		updateProduct(b_seq_qual,subsetproducts.b_seq_qual);
		updateProduct(b_seq_tags,subsetproducts.b_seq_tags);
		#endif
	};
	
	bool operator== (this_type const & other) const
	{
		return count==other.count 
			#if defined(BIOBAMBAM_HAVE_GMP)
			&& (mpz_cmp(gmpb_seq,other.gmpb_seq) == 0)
			&& (mpz_cmp(gmpname_b_seq,other.gmpname_b_seq) == 0)
			&& (mpz_cmp(gmpb_seq_qual,other.gmpb_seq_qual) == 0)
			&& (mpz_cmp(gmpb_seq_tags,other.gmpb_seq_tags) == 0)
			#else
			&& b_seq==other.b_seq 
			&& name_b_seq==other.name_b_seq 
			&& b_seq_qual==other.b_seq_qual 
			&& b_seq_tags==other.b_seq_tags
			#endif
			;
	};

	std::string get_b_seq() const
	{
		std::ostringstream ostr;
		ostr << std::hex 
			#if defined(BIOBAMBAM_HAVE_GMP)
			<< libmaus::math::UnsignedInteger<primeWidth/32>(convertNumber(gmpb_seq))
			#else
			<< libmaus::math::UnsignedInteger<primeWidth/32>(b_seq)
			#endif
			;
		return ostr.str();
	}

	std::string get_name_b_seq() const
	{
		std::ostringstream ostr;
		ostr << std::hex 
			#if defined(BIOBAMBAM_HAVE_GMP)
			<< libmaus::math::UnsignedInteger<primeWidth/32>(convertNumber(gmpname_b_seq))
			#else
			<< libmaus::math::UnsignedInteger<primeWidth/32>(name_b_seq)
			#endif
			;
		return ostr.str();
	}
	
	std::string get_b_seq_qual() const
	{
		std::ostringstream ostr;
		ostr << std::hex 
			#if defined(BIOBAMBAM_HAVE_GMP)
			<< libmaus::math::UnsignedInteger<primeWidth/32>(convertNumber(gmpb_seq_qual))
			#else
			<< libmaus::math::UnsignedInteger<primeWidth/32>(b_seq_qual)
			#endif
			;
		return ostr.str();
	}
	
	std::string get_b_seq_tags() const
	{
		std::ostringstream ostr;
		ostr << std::hex 
			#if defined(BIOBAMBAM_HAVE_GMP)
			<< libmaus::math::UnsignedInteger<primeWidth/32>(convertNumber(gmpb_seq_tags))
			#else
			<< libmaus::math::UnsignedInteger<primeWidth/32>(b_seq_tags)
			#endif
			;
		return ostr.str();
	}
	
	uint64_t get_count() const
	{
		return count;
	}
};

typedef PrimeProduct<CRC32UpdateContext,32> CRC32PrimeProduct32;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct32::productWidth> const CRC32PrimeProduct32::prime = getMersenne31<CRC32PrimeProduct32::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct32::productWidth> const CRC32PrimeProduct32::foldprime = 
	getMersenne31<CRC32PrimeProduct32::productWidth>() + libmaus::math::UnsignedInteger<CRC32PrimeProduct32::productWidth>(1);

typedef PrimeProduct<CRC32UpdateContext,64> CRC32PrimeProduct64;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct64::productWidth> const CRC32PrimeProduct64::prime = getPrime64<CRC32PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct64::productWidth> const CRC32PrimeProduct64::foldprime = getPrime64<CRC32PrimeProduct64::productWidth>();

typedef PrimeProduct<MD5UpdateContext,64> MD5PrimeProduct64;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct64::productWidth> const MD5PrimeProduct64::prime = getPrime64<MD5PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct64::productWidth> const MD5PrimeProduct64::foldprime = getPrime64<MD5PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,64> SHA1PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct64::productWidth> const SHA1PrimeProduct64::prime = getPrime64<SHA1PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct64::productWidth> const SHA1PrimeProduct64::foldprime = getPrime64<SHA1PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,64> SHA2_224_PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct64::productWidth> const SHA2_224_PrimeProduct64::prime = getPrime64<SHA2_224_PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct64::productWidth> const SHA2_224_PrimeProduct64::foldprime = getPrime64<SHA2_224_PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,64> SHA2_256_PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct64::productWidth> const SHA2_256_PrimeProduct64::prime = getPrime64<SHA2_256_PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct64::productWidth> const SHA2_256_PrimeProduct64::foldprime = getPrime64<SHA2_256_PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,64> SHA2_256_sse4_PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct64::productWidth> const SHA2_256_sse4_PrimeProduct64::prime = getPrime64<SHA2_256_sse4_PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct64::productWidth> const SHA2_256_sse4_PrimeProduct64::foldprime = getPrime64<SHA2_256_sse4_PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,64> SHA2_384_PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct64::productWidth> const SHA2_384_PrimeProduct64::prime = getPrime64<SHA2_384_PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct64::productWidth> const SHA2_384_PrimeProduct64::foldprime = getPrime64<SHA2_384_PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,64> SHA2_512_PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct64::productWidth> const SHA2_512_PrimeProduct64::prime = getPrime64<SHA2_512_PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct64::productWidth> const SHA2_512_PrimeProduct64::foldprime = getPrime64<SHA2_512_PrimeProduct64::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,64> SHA2_512_sse4_PrimeProduct64;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct64::productWidth> const SHA2_512_sse4_PrimeProduct64::prime = getPrime64<SHA2_512_sse4_PrimeProduct64::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct64::productWidth> const SHA2_512_sse4_PrimeProduct64::foldprime = getPrime64<SHA2_512_sse4_PrimeProduct64::productWidth>();

typedef PrimeProduct<CRC32UpdateContext,96> CRC32PrimeProduct96;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct96::productWidth> const CRC32PrimeProduct96::prime = getPrime96<CRC32PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct96::productWidth> const CRC32PrimeProduct96::foldprime = getPrime96<CRC32PrimeProduct96::productWidth>();

typedef PrimeProduct<MD5UpdateContext,96> MD5PrimeProduct96;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct96::productWidth> const MD5PrimeProduct96::prime = getPrime96<MD5PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct96::productWidth> const MD5PrimeProduct96::foldprime = getPrime96<MD5PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,96> SHA1PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct96::productWidth> const SHA1PrimeProduct96::prime = getPrime96<SHA1PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct96::productWidth> const SHA1PrimeProduct96::foldprime = getPrime96<SHA1PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,96> SHA2_224_PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct96::productWidth> const SHA2_224_PrimeProduct96::prime = getPrime96<SHA2_224_PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct96::productWidth> const SHA2_224_PrimeProduct96::foldprime = getPrime96<SHA2_224_PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,96> SHA2_256_PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct96::productWidth> const SHA2_256_PrimeProduct96::prime = getPrime96<SHA2_256_PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct96::productWidth> const SHA2_256_PrimeProduct96::foldprime = getPrime96<SHA2_256_PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,96> SHA2_256_sse4_PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct96::productWidth> const SHA2_256_sse4_PrimeProduct96::prime = getPrime96<SHA2_256_sse4_PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct96::productWidth> const SHA2_256_sse4_PrimeProduct96::foldprime = getPrime96<SHA2_256_sse4_PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,96> SHA2_384_PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct96::productWidth> const SHA2_384_PrimeProduct96::prime = getPrime96<SHA2_384_PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct96::productWidth> const SHA2_384_PrimeProduct96::foldprime = getPrime96<SHA2_384_PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,96> SHA2_512_PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct96::productWidth> const SHA2_512_PrimeProduct96::prime = getPrime96<SHA2_512_PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct96::productWidth> const SHA2_512_PrimeProduct96::foldprime = getPrime96<SHA2_512_PrimeProduct96::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,96> SHA2_512_sse4_PrimeProduct96;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct96::productWidth> const SHA2_512_sse4_PrimeProduct96::prime = getPrime96<SHA2_512_sse4_PrimeProduct96::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct96::productWidth> const SHA2_512_sse4_PrimeProduct96::foldprime = getPrime96<SHA2_512_sse4_PrimeProduct96::productWidth>();

typedef PrimeProduct<CRC32UpdateContext,128> CRC32PrimeProduct128;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct128::productWidth> const CRC32PrimeProduct128::prime = getPrime128<CRC32PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct128::productWidth> const CRC32PrimeProduct128::foldprime = getPrime128<CRC32PrimeProduct128::productWidth>();

typedef PrimeProduct<MD5UpdateContext,128> MD5PrimeProduct128;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct128::productWidth> const MD5PrimeProduct128::prime = getPrime128<MD5PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct128::productWidth> const MD5PrimeProduct128::foldprime = getPrime128<MD5PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,128> SHA1PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct128::productWidth> const SHA1PrimeProduct128::prime = getPrime128<SHA1PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct128::productWidth> const SHA1PrimeProduct128::foldprime = getPrime128<SHA1PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,128> SHA2_224_PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct128::productWidth> const SHA2_224_PrimeProduct128::prime = getPrime128<SHA2_224_PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct128::productWidth> const SHA2_224_PrimeProduct128::foldprime = getPrime128<SHA2_224_PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,128> SHA2_256_PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct128::productWidth> const SHA2_256_PrimeProduct128::prime = getPrime128<SHA2_256_PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct128::productWidth> const SHA2_256_PrimeProduct128::foldprime = getPrime128<SHA2_256_PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,128> SHA2_256_sse4_PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct128::productWidth> const SHA2_256_sse4_PrimeProduct128::prime = getPrime128<SHA2_256_sse4_PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct128::productWidth> const SHA2_256_sse4_PrimeProduct128::foldprime = getPrime128<SHA2_256_sse4_PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,128> SHA2_384_PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct128::productWidth> const SHA2_384_PrimeProduct128::prime = getPrime128<SHA2_384_PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct128::productWidth> const SHA2_384_PrimeProduct128::foldprime = getPrime128<SHA2_384_PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,128> SHA2_512_PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct128::productWidth> const SHA2_512_PrimeProduct128::prime = getPrime128<SHA2_512_PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct128::productWidth> const SHA2_512_PrimeProduct128::foldprime = getPrime128<SHA2_512_PrimeProduct128::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,128> SHA2_512_sse4_PrimeProduct128;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct128::productWidth> const SHA2_512_sse4_PrimeProduct128::prime = getPrime128<SHA2_512_sse4_PrimeProduct128::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct128::productWidth> const SHA2_512_sse4_PrimeProduct128::foldprime = getPrime128<SHA2_512_sse4_PrimeProduct128::productWidth>();

typedef PrimeProduct<CRC32UpdateContext,160> CRC32PrimeProduct160;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct160::productWidth> const CRC32PrimeProduct160::prime = getPrime160<CRC32PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct160::productWidth> const CRC32PrimeProduct160::foldprime = getPrime160<CRC32PrimeProduct160::productWidth>();

typedef PrimeProduct<MD5UpdateContext,160> MD5PrimeProduct160;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct160::productWidth> const MD5PrimeProduct160::prime = getPrime160<MD5PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct160::productWidth> const MD5PrimeProduct160::foldprime = getPrime160<MD5PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,160> SHA1PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct160::productWidth> const SHA1PrimeProduct160::prime = getPrime160<SHA1PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct160::productWidth> const SHA1PrimeProduct160::foldprime = getPrime160<SHA1PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,160> SHA2_224_PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct160::productWidth> const SHA2_224_PrimeProduct160::prime = getPrime160<SHA2_224_PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct160::productWidth> const SHA2_224_PrimeProduct160::foldprime = getPrime160<SHA2_224_PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,160> SHA2_256_PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct160::productWidth> const SHA2_256_PrimeProduct160::prime = getPrime160<SHA2_256_PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct160::productWidth> const SHA2_256_PrimeProduct160::foldprime = getPrime160<SHA2_256_PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,160> SHA2_256_sse4_PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct160::productWidth> const SHA2_256_sse4_PrimeProduct160::prime = getPrime160<SHA2_256_sse4_PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct160::productWidth> const SHA2_256_sse4_PrimeProduct160::foldprime = getPrime160<SHA2_256_sse4_PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,160> SHA2_384_PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct160::productWidth> const SHA2_384_PrimeProduct160::prime = getPrime160<SHA2_384_PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct160::productWidth> const SHA2_384_PrimeProduct160::foldprime = getPrime160<SHA2_384_PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,160> SHA2_512_PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct160::productWidth> const SHA2_512_PrimeProduct160::prime = getPrime160<SHA2_512_PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct160::productWidth> const SHA2_512_PrimeProduct160::foldprime = getPrime160<SHA2_512_PrimeProduct160::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,160> SHA2_512_sse4_PrimeProduct160;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct160::productWidth> const SHA2_512_sse4_PrimeProduct160::prime = getPrime160<SHA2_512_sse4_PrimeProduct160::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct160::productWidth> const SHA2_512_sse4_PrimeProduct160::foldprime = getPrime160<SHA2_512_sse4_PrimeProduct160::productWidth>();

typedef PrimeProduct<CRC32UpdateContext,192> CRC32PrimeProduct192;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct192::productWidth> const CRC32PrimeProduct192::prime = getPrime192<CRC32PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct192::productWidth> const CRC32PrimeProduct192::foldprime = getPrime192<CRC32PrimeProduct192::productWidth>();

typedef PrimeProduct<MD5UpdateContext,192> MD5PrimeProduct192;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct192::productWidth> const MD5PrimeProduct192::prime = getPrime192<MD5PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct192::productWidth> const MD5PrimeProduct192::foldprime = getPrime192<MD5PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,192> SHA1PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct192::productWidth> const SHA1PrimeProduct192::prime = getPrime192<SHA1PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct192::productWidth> const SHA1PrimeProduct192::foldprime = getPrime192<SHA1PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,192> SHA2_224_PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct192::productWidth> const SHA2_224_PrimeProduct192::prime = getPrime192<SHA2_224_PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct192::productWidth> const SHA2_224_PrimeProduct192::foldprime = getPrime192<SHA2_224_PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,192> SHA2_256_PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct192::productWidth> const SHA2_256_PrimeProduct192::prime = getPrime192<SHA2_256_PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct192::productWidth> const SHA2_256_PrimeProduct192::foldprime = getPrime192<SHA2_256_PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,192> SHA2_256_sse4_PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct192::productWidth> const SHA2_256_sse4_PrimeProduct192::prime = getPrime192<SHA2_256_sse4_PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct192::productWidth> const SHA2_256_sse4_PrimeProduct192::foldprime = getPrime192<SHA2_256_sse4_PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,192> SHA2_384_PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct192::productWidth> const SHA2_384_PrimeProduct192::prime = getPrime192<SHA2_384_PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct192::productWidth> const SHA2_384_PrimeProduct192::foldprime = getPrime192<SHA2_384_PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,192> SHA2_512_PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct192::productWidth> const SHA2_512_PrimeProduct192::prime = getPrime192<SHA2_512_PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct192::productWidth> const SHA2_512_PrimeProduct192::foldprime = getPrime192<SHA2_512_PrimeProduct192::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,192> SHA2_512_sse4_PrimeProduct192;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct192::productWidth> const SHA2_512_sse4_PrimeProduct192::prime = getPrime192<SHA2_512_sse4_PrimeProduct192::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct192::productWidth> const SHA2_512_sse4_PrimeProduct192::foldprime = getPrime192<SHA2_512_sse4_PrimeProduct192::productWidth>();

typedef PrimeProduct<CRC32UpdateContext,224> CRC32PrimeProduct224;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct224::productWidth> const CRC32PrimeProduct224::prime = getPrime224<CRC32PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct224::productWidth> const CRC32PrimeProduct224::foldprime = getPrime224<CRC32PrimeProduct224::productWidth>();

typedef PrimeProduct<MD5UpdateContext,224> MD5PrimeProduct224;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct224::productWidth> const MD5PrimeProduct224::prime = getPrime224<MD5PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct224::productWidth> const MD5PrimeProduct224::foldprime = getPrime224<MD5PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,224> SHA1PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct224::productWidth> const SHA1PrimeProduct224::prime = getPrime224<SHA1PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct224::productWidth> const SHA1PrimeProduct224::foldprime = getPrime224<SHA1PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,224> SHA2_224_PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct224::productWidth> const SHA2_224_PrimeProduct224::prime = getPrime224<SHA2_224_PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct224::productWidth> const SHA2_224_PrimeProduct224::foldprime = getPrime224<SHA2_224_PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,224> SHA2_256_PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct224::productWidth> const SHA2_256_PrimeProduct224::prime = getPrime224<SHA2_256_PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct224::productWidth> const SHA2_256_PrimeProduct224::foldprime = getPrime224<SHA2_256_PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,224> SHA2_256_sse4_PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct224::productWidth> const SHA2_256_sse4_PrimeProduct224::prime = getPrime224<SHA2_256_sse4_PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct224::productWidth> const SHA2_256_sse4_PrimeProduct224::foldprime = getPrime224<SHA2_256_sse4_PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,224> SHA2_384_PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct224::productWidth> const SHA2_384_PrimeProduct224::prime = getPrime224<SHA2_384_PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct224::productWidth> const SHA2_384_PrimeProduct224::foldprime = getPrime224<SHA2_384_PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,224> SHA2_512_PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct224::productWidth> const SHA2_512_PrimeProduct224::prime = getPrime224<SHA2_512_PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct224::productWidth> const SHA2_512_PrimeProduct224::foldprime = getPrime224<SHA2_512_PrimeProduct224::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,224> SHA2_512_sse4_PrimeProduct224;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct224::productWidth> const SHA2_512_sse4_PrimeProduct224::prime = getPrime224<SHA2_512_sse4_PrimeProduct224::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct224::productWidth> const SHA2_512_sse4_PrimeProduct224::foldprime = getPrime224<SHA2_512_sse4_PrimeProduct224::productWidth>();

typedef PrimeProduct<CRC32UpdateContext,256> CRC32PrimeProduct256;
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct256::productWidth> const CRC32PrimeProduct256::prime = getPrime256<CRC32PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<CRC32PrimeProduct256::productWidth> const CRC32PrimeProduct256::foldprime = getPrime256<CRC32PrimeProduct256::productWidth>();

typedef PrimeProduct<MD5UpdateContext,256> MD5PrimeProduct256;
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct256::productWidth> const MD5PrimeProduct256::prime = getPrime256<MD5PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<MD5PrimeProduct256::productWidth> const MD5PrimeProduct256::foldprime = getPrime256<MD5PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA1UpdateContext,256> SHA1PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct256::productWidth> const SHA1PrimeProduct256::prime = getPrime256<SHA1PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA1PrimeProduct256::productWidth> const SHA1PrimeProduct256::foldprime = getPrime256<SHA1PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA2_224_UpdateContext,256> SHA2_224_PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct256::productWidth> const SHA2_224_PrimeProduct256::prime = getPrime256<SHA2_224_PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_224_PrimeProduct256::productWidth> const SHA2_224_PrimeProduct256::foldprime = getPrime256<SHA2_224_PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA2_256_UpdateContext,256> SHA2_256_PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct256::productWidth> const SHA2_256_PrimeProduct256::prime = getPrime256<SHA2_256_PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_PrimeProduct256::productWidth> const SHA2_256_PrimeProduct256::foldprime = getPrime256<SHA2_256_PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA2_256_sse4_UpdateContext,256> SHA2_256_sse4_PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct256::productWidth> const SHA2_256_sse4_PrimeProduct256::prime = getPrime256<SHA2_256_sse4_PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_256_sse4_PrimeProduct256::productWidth> const SHA2_256_sse4_PrimeProduct256::foldprime = getPrime256<SHA2_256_sse4_PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA2_384_UpdateContext,256> SHA2_384_PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct256::productWidth> const SHA2_384_PrimeProduct256::prime = getPrime256<SHA2_384_PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_384_PrimeProduct256::productWidth> const SHA2_384_PrimeProduct256::foldprime = getPrime256<SHA2_384_PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA2_512_UpdateContext,256> SHA2_512_PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct256::productWidth> const SHA2_512_PrimeProduct256::prime = getPrime256<SHA2_512_PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeProduct256::productWidth> const SHA2_512_PrimeProduct256::foldprime = getPrime256<SHA2_512_PrimeProduct256::productWidth>();

typedef PrimeProduct<SHA2_512_sse4_UpdateContext,256> SHA2_512_sse4_PrimeProduct256;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct256::productWidth> const SHA2_512_sse4_PrimeProduct256::prime = getPrime256<SHA2_512_sse4_PrimeProduct256::productWidth>();
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeProduct256::productWidth> const SHA2_512_sse4_PrimeProduct256::foldprime = getPrime256<SHA2_512_sse4_PrimeProduct256::productWidth>();

template<typename _context_type, size_t _primeWidth, size_t truncate>
struct PrimeSums
{
	typedef _context_type context_type;
	typedef PrimeSums<context_type,_primeWidth,truncate> this_type;
	
	// width of the prime in bits
	static size_t const primeWidth = _primeWidth;
	// width of the sums in bits
	static size_t const sumWidth = 1 + TemplateMax<primeWidth/32, context_type::digest_type::digestlength/4>::value;
	// prime
	static libmaus::math::UnsignedInteger<sumWidth> const prime;

	private:
	uint64_t count;
	uint64_t modcount;

	libmaus::math::UnsignedInteger<sumWidth> b_seq;
	libmaus::math::UnsignedInteger<sumWidth> name_b_seq;
	libmaus::math::UnsignedInteger<sumWidth> b_seq_qual;
	libmaus::math::UnsignedInteger<sumWidth> b_seq_tags;
	
	// get prime width in bits
	unsigned int getPrimeWidth() const
	{
		for ( unsigned int i = 0; i < sumWidth; ++i )
		{
			// most significant to least significant word
			uint32_t v = prime[sumWidth-i-1];
			
			// if word is not null
			if ( v )
			{
				// get word shift
				unsigned int s = (sumWidth-i-1)*32;
				
				if ( (v & 0xFFFF0000UL) )
					v >>= 16, s += 16;
				if ( (v & 0xFF00UL) )
					v >>= 8, s += 8;
				if ( (v & 0xF0UL) )
					v >>= 4, s += 4;
				if ( (v & 0xCUL) )
					v >>= 2, s += 2;
								
				// add bit shift
				while ( v )
				{
					s++;
					v >>= 1;
				}
				
				return s;
			}
		}
		
		return 0;
	}
	
	// get prime width in hex digits
	unsigned int getPrimeHexWidth() const
	{
		unsigned int const primewidth = getPrimeWidth();
		return (primewidth + 3)/4;
	}
	
	std::string shortenHexDigest(std::string const & s) const
	{
		if ( truncate )
		{
			return s.substr(s.size() - (truncate/32)*8);
		}
		else
		{
			assert ( s.size() >= getPrimeHexWidth() );
			return s.substr(s.size()-getPrimeHexWidth());
		}
	}

	template<size_t k>
	void updateSum(libmaus::math::UnsignedInteger<sumWidth> & sum, libmaus::math::UnsignedInteger<k> const & crc)
	{
		// prime must be larger than crc to guarantee that (crc%prime) != 0 assuming crc!=0
		if ( expect_false(crc.isNull()) )
			sum += libmaus::math::UnsignedInteger<sumWidth>(1);
		else
			sum += crc;

		if ( expect_false ( (++modcount & ((1ul<<30)-1)) == 0 ) )
		{
			sum %= prime;
			modcount = 0;
		}
	}

	public:
	PrimeSums() : count(0), modcount(0), b_seq(0), name_b_seq(0), b_seq_qual(0), b_seq_tags(0) {}
	
	void push(context_type const & context)
	{
		count += 1;		
		updateSum(b_seq     ,context.flags_seq_digest);
		updateSum(name_b_seq,context.name_flags_seq_digest);
		updateSum(b_seq_qual,context.flags_seq_qual_digest);
		updateSum(b_seq_tags,context.flags_seq_tags_digest);
	}

	void push (this_type const & subsetproducts)
	{
		count += subsetproducts.count;
		b_seq      = ((b_seq      % prime) + (subsetproducts.b_seq      % prime)) % prime;
		name_b_seq = ((name_b_seq % prime) + (subsetproducts.name_b_seq % prime)) % prime;
		b_seq_qual = ((b_seq_qual % prime) + (subsetproducts.b_seq_qual % prime)) % prime;
		b_seq_tags = ((b_seq_tags % prime) + (subsetproducts.b_seq_tags % prime)) % prime;		
	};
	
	bool operator== (this_type const & other) const
	{
		return count==other.count 
			&& (b_seq%prime)==(other.b_seq%prime)
			&& (name_b_seq%prime)==(other.name_b_seq%prime)
			&& (b_seq_qual%prime)==(other.b_seq_qual%prime)
			&& (b_seq_tags%prime)==(other.b_seq_tags%prime)
			;
	};

	std::string get_b_seq() const
	{
		std::ostringstream ostr;
		ostr << std::hex << libmaus::math::UnsignedInteger<primeWidth/32>(b_seq % prime);
		return shortenHexDigest(ostr.str());
	}

	std::string get_name_b_seq() const
	{
		std::ostringstream ostr;
		ostr << std::hex << libmaus::math::UnsignedInteger<primeWidth/32>(name_b_seq % prime);
		return shortenHexDigest(ostr.str());
	}
	
	std::string get_b_seq_qual() const
	{
		std::ostringstream ostr;
		ostr << std::hex << libmaus::math::UnsignedInteger<primeWidth/32>(b_seq_qual % prime);
		return shortenHexDigest(ostr.str());
	}
	
	std::string get_b_seq_tags() const
	{
		std::ostringstream ostr;
		ostr << std::hex << libmaus::math::UnsignedInteger<primeWidth/32>(b_seq_tags % prime);
		return shortenHexDigest(ostr.str());
	}
	
	uint64_t get_count() const
	{
		return count;
	}
};

template<size_t k>
static libmaus::math::UnsignedInteger<k> getNextPrime512()
{
	libmaus::math::UnsignedInteger<k> U(1);
	U <<= 512;
	U += libmaus::math::UnsignedInteger<k>(75);
	return U;
}

typedef PrimeSums<SHA2_512_UpdateContext,544,512> SHA2_512_PrimeSums512;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeSums512::sumWidth> const SHA2_512_PrimeSums512::prime = getNextPrime512<SHA2_512_PrimeSums512::sumWidth>();

typedef PrimeSums<SHA2_512_UpdateContext,544,0> SHA2_512_PrimeSums;
template<> libmaus::math::UnsignedInteger<SHA2_512_PrimeSums::sumWidth> const SHA2_512_PrimeSums::prime = getMersenne521<SHA2_512_PrimeSums::sumWidth>();

typedef PrimeSums<SHA2_512_sse4_UpdateContext,544,512> SHA2_512_sse4_PrimeSums512;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeSums512::sumWidth> const SHA2_512_sse4_PrimeSums512::prime = getNextPrime512<SHA2_512_sse4_PrimeSums512::sumWidth>();

typedef PrimeSums<SHA2_512_sse4_UpdateContext,544,0> SHA2_512_sse4_PrimeSums;
template<> libmaus::math::UnsignedInteger<SHA2_512_sse4_PrimeSums::sumWidth> const SHA2_512_sse4_PrimeSums::prime = getMersenne521<SHA2_512_sse4_PrimeSums::sumWidth>();

/**
* Product checksums calculated based on basecalls and (multi segment, first
* and last) bit info, and this combined with the query name, or the basecall
* qualities, or certain BAM auxilary fields.
**/
struct CRC32Products 
{
	typedef CRC32UpdateContext context_type;

	private:
	uint64_t count;
	uint64_t b_seq;
	uint64_t name_b_seq;
	uint64_t b_seq_qual;
	uint64_t b_seq_tags;

	/**
	* Multiply existing product by new checksum, having altered checksum ready for
	* multiplication in a finite field i.e. 0 < chksum < (2^31 -1)
	* @param product to be updated
	* @param chksum to update product with
	**/
	static void product_munged_chksum_multiply (uint64_t & product, uint32_t chksum) {
		static uint64_t const MERSENNE31 = 0x7FFFFFFFull; // Mersenne Prime 2^31 - 1
		chksum &= MERSENNE31;
		if (!chksum || chksum == MERSENNE31 ) chksum = 1;
		product = ( product * chksum ) % MERSENNE31;
	}

	void push_count(uint64_t const add)
	{
		count += add;
	}
	
	void push_b_seq(uint32_t const chksum)
	{
		product_munged_chksum_multiply(b_seq,chksum);
	}
	
	void push_name_b_seq(uint32_t const chksum)
	{
		product_munged_chksum_multiply(name_b_seq,chksum);
	}
	
	void push_b_seq_qual(uint32_t const chksum)
	{
		product_munged_chksum_multiply(b_seq_qual,chksum);
	}
	
	void push_b_seq_tags(uint32_t const chksum)
	{
		product_munged_chksum_multiply(b_seq_tags,chksum);
	}

	void push(
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const & name_b_seq,
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const & b_seq,
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const & b_seq_qual,
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const & b_seq_tags
	)
	{
		push_count(1);
		push_name_b_seq(name_b_seq[0]);
		push_b_seq(b_seq[0]);
		push_b_seq_qual(b_seq_qual[0]);
		push_b_seq_tags(b_seq_tags[0]);
	}
	
	public:
	CRC32Products() : count(0), b_seq(1), name_b_seq(1), b_seq_qual(1), b_seq_tags(1)
	{
	}

	uint64_t get_b_seq() const
	{
		return b_seq;
	}

	uint64_t get_name_b_seq() const
	{
		return name_b_seq;
	}
	
	uint64_t get_b_seq_qual() const
	{
		return b_seq_qual;
	}
	
	uint64_t get_b_seq_tags() const
	{
		return b_seq_tags;
	}
	
	uint64_t get_count() const
	{
		return count;
	}
	
	void push(context_type const & context)
	{
		push(
			context.name_flags_seq_digest,
			context.flags_seq_digest,
			context.flags_seq_qual_digest,
			context.flags_seq_tags_digest
		);
	}

	void push (CRC32Products const & subsetproducts)
	{
		count += subsetproducts.count;
		product_munged_chksum_multiply(b_seq, subsetproducts.b_seq);
		product_munged_chksum_multiply(name_b_seq, subsetproducts.name_b_seq);
		product_munged_chksum_multiply(b_seq_qual, subsetproducts.b_seq_qual);
		product_munged_chksum_multiply(b_seq_tags, subsetproducts.b_seq_tags);
	};
	bool operator== (CRC32Products const & other) const
	{
		return count==other.count &&
			b_seq==other.b_seq && name_b_seq==other.name_b_seq &&
			b_seq_qual==other.b_seq_qual && b_seq_tags==other.b_seq_tags;
	};
};

/**
 * null checksums for performance testing
**/
struct NullChecksums
{
	typedef NullUpdateContext context_type;

	private:
	uint64_t count;

	void push_count(uint64_t const add)
	{
		count += add;
	}
	
	void push(
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const &,
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const &,
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const &,
		libmaus::math::UnsignedInteger<context_type::digest_type::digestlength/4> const &
	)
	{
		push_count(1);
	}
	
	public:
	NullChecksums() : count(0)
	{
	}

	uint64_t get_b_seq() const
	{
		return 0;
	}

	uint64_t get_name_b_seq() const
	{
		return 0;
	}
	
	uint64_t get_b_seq_qual() const
	{
		return 0;
	}
	
	uint64_t get_b_seq_tags() const
	{
		return 0;
	}
	
	uint64_t get_count() const
	{
		return count;
	}
	
	void push(context_type const & context)
	{
		push(
			context.name_flags_seq_digest,
			context.flags_seq_digest,
			context.flags_seq_qual_digest,
			context.flags_seq_tags_digest
		);
	}

	void push (NullChecksums const & subsetproducts)
	{
		count += subsetproducts.count;
	};
	bool operator== (NullChecksums const & other) const
	{
		return count==other.count;
	};
};

/**
* Finite field products of CRC32 checksums of primary/source sequence data
*
* Checksum products should remain unchanged if primary data is retained no matter what
* alignment information is added or altered, or the ordering of the records.
*
* Products calulated for all records and for those with pass QC bit.
**/
template<typename _crc_container_type>
struct OrderIndependentSeqDataChecksums {
	typedef _crc_container_type crc_container_type;
	typedef typename crc_container_type::context_type context_type;

	private:
	::libmaus::autoarray::AutoArray<char> A; // check with German: can we change/treat underlying data block to unsigned char / uint8_t?
	::libmaus::autoarray::AutoArray<char> B; //separate A & B can allow reordering speed up?
	public:
	std::vector<std::string> const auxtags;
	// to provide fast checking of aux fields when processing and validation at startup
	::libmaus::bambam::BamAuxFilterVector const auxtagsfilter;
	crc_container_type all;
	crc_container_type pass;
	OrderIndependentSeqDataChecksums() : A(), B(), auxtags(getDefaultAuxTags()), auxtagsfilter(auxtags), all(), pass() { };
	/**
	* Combine primary sequence data from alignment record into checksum products
	*
	* Ignores Supplementary and auxilary alignment records so primary data in not 
	* included twice.
	*
	* @param algn BAM alignment record
	**/
	void push(libmaus::bambam::BamAlignment const & algn, typename crc_container_type::context_type & context)
	{
		if 
		( ! (
			algn.getFlags() &
			(
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY
			)
		) )
		{
			static uint16_t const maskflags =
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 | 
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2;
				
			uint8_t const flags = (algn.getFlags() & maskflags) & 0xFF;
				
			context.pass = ! algn.isQCFail();
			context.valid = true;
			uint64_t const len = algn.isReverse() ? algn.decodeReadRC(A) : algn.decodeRead(A);
			
			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t const CRC32_INITIAL = crc32(0L, Z_NULL, 0);
			#endif

			context.ctx_name_flags_seq.init();
			// read name (algn.getLReadName() bytes)
			context.ctx_name_flags_seq.update(reinterpret_cast<uint8_t const *>(algn.getName()),algn.getLReadName());
			// flags (1 byte)
			context.ctx_name_flags_seq.update(reinterpret_cast<uint8_t const *>(&flags),1);
			// query sequence
			context.ctx_name_flags_seq.update(reinterpret_cast<uint8_t const *>(A.begin()),len);
						
			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t chksum_name_flags_seq = crc32(CRC32_INITIAL,reinterpret_cast<const unsigned char *>(algn.getName()),algn.getLReadName());
			chksum_name_flags_seq = crc32(chksum_name_flags_seq,reinterpret_cast<const unsigned char*>( &flags), 1);
			chksum_name_flags_seq = crc32(chksum_name_flags_seq,reinterpret_cast<const unsigned char *>( A.begin()), len);
			#endif
						
			context.ctx_flags_seq.init();
			// flags
			context.ctx_flags_seq.update(reinterpret_cast<uint8_t const *>(&flags),1);
			// query sequence
			context.ctx_flags_seq.update(reinterpret_cast<uint8_t const *>(A.begin()),len);
			
			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t chksum_flags_seq = crc32(CRC32_INITIAL,reinterpret_cast<const unsigned char*>(&flags), 1);
			chksum_flags_seq = crc32(chksum_flags_seq,reinterpret_cast<const unsigned char *>( A.begin()), len);
			#endif

			// add quality
			uint64_t const lenq = algn.isReverse() ? algn.decodeQualRC(B) : algn.decodeQual(B);
			
			context.ctx_flags_seq_qual.copyFrom(context.ctx_flags_seq);
			context.ctx_flags_seq_qual.update(reinterpret_cast<uint8_t *>(B.begin()), lenq);

			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t chksum_flags_seq_qual = chksum_flags_seq;
			chksum_flags_seq_qual = crc32(chksum_flags_seq_qual,reinterpret_cast<const unsigned char *>( B.begin()), lenq);
			#endif

			// set aux to start pointer of aux area
			uint8_t const * aux = ::libmaus::bambam::BamAlignmentDecoderBase::getAux(algn.D.begin());
			// end of algn data block (and so of aux area)
			uint8_t const * const auxend = algn.D.begin() + algn.blocksize;
			
			// start of each aux entry for auxtags
			std::vector<uint8_t const *> paux(auxtags.size(),NULL);
			// size of each aux entry for auxtags
			std::vector<uint64_t> saux(auxtags.size(),0);
			
			while(aux < auxend)
			{
				// get length of aux field in bytes
				uint64_t const auxlen = ::libmaus::bambam::BamAlignmentDecoderBase::getAuxLength(aux);
				
				// skip over aux tag data we're not interested in
				if( auxtagsfilter(aux[0],aux[1]) )
				{
					std::vector<std::string>::const_iterator it_auxtags = auxtags.begin();
					std::vector<uint8_t const *>::iterator it_paux = paux.begin();
					std::vector<uint64_t>::iterator it_saux = saux.begin();
					// consider each tag in auxtag
					while ( it_auxtags != auxtags.end())
					{
						// until we match the tag for the current aux data
						if(! it_auxtags->compare(0,2,reinterpret_cast<const char *>(aux),2) )
						{
							*it_paux = aux;
							*it_saux = auxlen;
							// stop this loop
							it_auxtags = auxtags.end();
						}
						else
						{
							++it_auxtags;
							++it_paux;
							++it_saux;
						}
					}
				}
				aux += auxlen;
			}
			std::vector<uint8_t const *>::const_iterator it_paux = paux.begin();
			std::vector<uint64_t>::const_iterator it_saux = saux.begin();
			
			// copy context
			context.ctx_flags_seq_tags.copyFrom(context.ctx_flags_seq);

			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t chksum_flags_seq_tags = chksum_flags_seq;
			#endif
			
			//loop over the chunks of data corresponding to the auxtags
			while (it_paux != paux.end())
			{
				//if data exists push into running checksum
				if(*it_paux)
				{
					context.ctx_flags_seq_tags.update(reinterpret_cast<uint8_t const *>(*it_paux), *it_saux);
					#if defined(BAM_SEQ_CHKSUM_DEBUG)
					chksum_flags_seq_tags = crc32(chksum_flags_seq_tags,reinterpret_cast<const unsigned char *>(*it_paux), *it_saux);
					#endif
				}
				++it_paux;
				++it_saux;
			}
	

			context.flags_seq_qual_digest = context.ctx_flags_seq_qual.digestui();
			context.name_flags_seq_digest = context.ctx_name_flags_seq.digestui();
			context.flags_seq_digest = context.ctx_flags_seq.digestui();
			context.flags_seq_tags_digest = context.ctx_flags_seq_tags.digestui();

			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			assert ( context.flags_seq_qual_digest[0] == chksum_flags_seq_qual );			
			assert ( context.name_flags_seq_digest[0] == chksum_name_flags_seq );
			assert ( context.flags_seq_digest[0] == chksum_flags_seq );
			assert ( context.flags_seq_tags_digest[0] == chksum_flags_seq_tags );
			#endif

			push(context);
		}
		else
		{
			context.valid = false;
		}
	};
	
	void push(typename crc_container_type::context_type const & context)
	{
		if ( context.valid )
		{
			all.push(context);
			
			if ( context.pass )
				pass.push(context);	
		}
	}
	
	void push(OrderIndependentSeqDataChecksums const & subsetchksum)
	{
		pass.push(subsetchksum.pass);
		all.push(subsetchksum.all);
	};
	bool operator== (OrderIndependentSeqDataChecksums const & other) const
	{
		return pass==other.pass && all==other.all && auxtags==other.auxtags;
	}
};

namespace libmaus
{
	namespace autoarray
	{
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<NullChecksums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<NullChecksums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32Products> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32Products> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_SimpleSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_SimpleSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct32> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct32> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct64> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct64> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct96> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct96> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct128> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct128> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct160> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct160> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct192> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct192> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct224> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct224> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<CRC32PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<CRC32PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<MD5PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<MD5PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA1PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA1PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_224_PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_256_sse4_PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_384_PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct256> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeProduct256> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeSums> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeSums> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_PrimeSums512> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_PrimeSums512> *, uint64_t const)
			{
			
			}
		};
		template<>
		struct ArrayErase<OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeSums512> >
		{
			static void erase(OrderIndependentSeqDataChecksums<SHA2_512_sse4_PrimeSums512> *, uint64_t const)
			{
			
			}
		};
	}
}

template<typename container_type>
int bamseqchksumTemplate(::libmaus::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	libmaus::timing::RealTimeClock rtc;
	rtc.start();
	double prevtime = 0;
	
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();

	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	::libmaus::bambam::BamAlignment & algn = dec.getAlignment();
	OrderIndependentSeqDataChecksums<container_type> chksums;
	libmaus::autoarray::AutoArray< OrderIndependentSeqDataChecksums<container_type> > readgroup_chksums(1 + header.getNumReadGroups(),false);
	typename OrderIndependentSeqDataChecksums<container_type>::context_type updatecontext;

	uint64_t c = 0;
	while ( dec.readAlignment() )
	{
		chksums.push(algn,updatecontext);
		readgroup_chksums[algn.getReadGroupId(header)+1].push(updatecontext);
		
		if ( verbose && (++c & (1024*1024-1)) == 0 )
		{
			double const elapsed = rtc.getElapsedSeconds();
			std::cerr << "[V] " << c/(1024*1024) << " " << chksums.all.get_count() << " " << algn.getName() << " " << algn.isRead1()
			<< " " << algn.isRead2() << " " << ( algn.isReverse() ? algn.getReadRC() : algn.getRead() ) << " "
			<< ( algn.isReverse() ? algn.getQualRC() : algn.getQual() ) << " " << std::hex << (0x0 + algn.getFlags())
			<< std::dec << " " << chksums.all.get_b_seq() << " " << chksums.all.get_b_seq_tags() << " " 
			<< " " << (elapsed-prevtime)
			<< std::endl;
			prevtime = elapsed;
		}
	}

	std::cout << "###\tset\t" << "count" << "\t" << "\t" << "b_seq" << "\t" << "name_b_seq" << "\t" << "b_seq_qual" << "\t"
		<< "b_seq_tags(";
	for (std::vector<std::string>::const_iterator it_aux = chksums.auxtags.begin(); it_aux != chksums.auxtags.end(); ++it_aux)
	{
		if (chksums.auxtags.begin() != it_aux)
			std::cout << ",";
		std::cout << *it_aux;
	}
	std::cout << ")" << std::endl;
	std::cout << "all\tall\t" << chksums.all.get_count() << "\t" << std::hex << "\t" << chksums.all.get_b_seq() << "\t"
		<< chksums.all.get_name_b_seq() << "\t" << chksums.all.get_b_seq_qual() << "\t" << chksums.all.get_b_seq_tags() << std::dec << std::endl;
	std::cout << "all\tpass\t" << chksums.pass.get_count() << "\t" << std::hex << "\t" << chksums.pass.get_b_seq() << "\t"
		<< chksums.pass.get_name_b_seq() << "\t" << chksums.pass.get_b_seq_qual() << "\t" << chksums.pass.get_b_seq_tags() << std::dec << std::endl;

	if(header.getNumReadGroups())
	{
		OrderIndependentSeqDataChecksums<container_type> chksumschk;
		for(unsigned int i=0; i<=header.getNumReadGroups(); i++)
		{
			chksumschk.push(readgroup_chksums[i]);
			std::cout << (i>0 ? header.getReadGroups().at(i-1).ID : "") << "\tall\t" << readgroup_chksums[i].all.get_count() << "\t"
				<< std::hex << "\t" << readgroup_chksums[i].all.get_b_seq() << "\t" << readgroup_chksums[i].all.get_name_b_seq() << "\t"
				<< readgroup_chksums[i].all.get_b_seq_qual() << "\t" << readgroup_chksums[i].all.get_b_seq_tags() << std::dec << std::endl;
			std::cout << (i>0 ? header.getReadGroups().at(i-1).ID : "") << "\tpass\t" << readgroup_chksums[i].pass.get_count() << "\t"
				<< std::hex << "\t" << readgroup_chksums[i].pass.get_b_seq() << "\t" << readgroup_chksums[i].pass.get_name_b_seq() << "\t"
				<< readgroup_chksums[i].pass.get_b_seq_qual() << "\t" << readgroup_chksums[i].pass.get_b_seq_tags() << std::dec << std::endl;
		}
		assert(chksumschk == chksums);
	}
	
	if ( verbose )
	{
		std::cerr << "[V] run time " << rtc.getElapsedSeconds() << " (" << rtc.formatTime(rtc.getElapsedSeconds()) << ")" << std::endl;
	}

	return EXIT_SUCCESS;
}

static std::string getDefaultHash()
{
	return "crc32prod";
}

static std::vector<std::string> getSupportedHashVariants()
{
	std::vector<std::string> V;
	V.push_back("crc32prod");
	V.push_back("crc32");
	V.push_back("md5");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1");
	V.push_back("sha224");
	V.push_back("sha256");
	V.push_back("sha384");
	V.push_back("sha512");
	#endif
	V.push_back("crc32prime32");
	V.push_back("crc32prime64");
	V.push_back("md5prime64");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime64");
	V.push_back("sha224prime64");
	V.push_back("sha256prime64");
	V.push_back("sha384prime64");
	V.push_back("sha512prime64");
	#endif
	V.push_back("crc32prime96");
	V.push_back("md5prime96");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime96");
	V.push_back("sha224prime96");
	V.push_back("sha256prime96");
	V.push_back("sha384prime96");
	V.push_back("sha512prime96");
	#endif
	V.push_back("crc32prime128");
	V.push_back("md5prime128");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime128");
	V.push_back("sha224prime128");
	V.push_back("sha256prime128");
	V.push_back("sha384prime128");
	V.push_back("sha512prime128");
	#endif
	V.push_back("crc32prime160");
	V.push_back("md5prime160");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime160");
	V.push_back("sha224prime160");
	V.push_back("sha256prime160");
	V.push_back("sha384prime160");
	V.push_back("sha512prime160");
	#endif
	V.push_back("crc32prime192");
	V.push_back("md5prime192");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime192");
	V.push_back("sha224prime192");
	V.push_back("sha256prime192");
	V.push_back("sha384prime192");
	V.push_back("sha512prime192");
	#endif
	V.push_back("crc32prime224");
	V.push_back("md5prime224");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime224");
	V.push_back("sha224prime224");
	V.push_back("sha256prime224");
	V.push_back("sha384prime224");
	V.push_back("sha512prime224");
	#endif
	V.push_back("crc32prime256");
	V.push_back("md5prime256");
	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha1prime256");
	V.push_back("sha224prime256");
	V.push_back("sha256prime256");
	V.push_back("sha384prime256");
	V.push_back("sha512prime256");
	#endif
	V.push_back("null");

	#if defined(LIBMAUS_HAVE_NETTLE)
	V.push_back("sha512primesums");
	V.push_back("sha512primesums512");
	#endif

	#if (!defined(LIBMAUS_HAVE_NETTLE)) && (defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386) && defined(LIBMAUS_HAVE_SHA2_ASSEMBLY))
	V.push_back("sha256");
	V.push_back("sha512");
	V.push_back("sha256prime64");
	V.push_back("sha512prime64");
	V.push_back("sha256prime96");
	V.push_back("sha512prime96");
	V.push_back("sha256prime128");
	V.push_back("sha512prime128");
	V.push_back("sha256prime160");
	V.push_back("sha512prime160");
	V.push_back("sha256prime192");
	V.push_back("sha512prime192");
	V.push_back("sha256prime224");
	V.push_back("sha512prime224");
	V.push_back("sha256prime256");
	V.push_back("sha512prime256");
	V.push_back("sha512primesums");
	V.push_back("sha512primesums512");
	#endif
	
	return V;
}

static std::string getSupportedHashVariantsList()
{
	std::ostringstream ostr;
	std::vector<std::string> const V = getSupportedHashVariants();
	
	if ( V.size() )
	{
		ostr << V[0];
		for ( uint64_t i = 1; i < V.size(); ++i )
			ostr << "," << V[i];
	}
	
	return ostr.str();
}

int bamseqchksum(::libmaus::util::ArgInfo const & arginfo)
{
	std::string const hash = arginfo.getValue<std::string>("hash",getDefaultHash());
	
	if ( hash == "crc32prod" )
	{
		return bamseqchksumTemplate<CRC32Products>(arginfo);
	}
	else if ( hash == "crc32" )
	{
		return bamseqchksumTemplate<CRC32SimpleSums>(arginfo);
	}
	else if ( hash == "md5" )
	{
		return bamseqchksumTemplate<MD5SimpleSums>(arginfo);
	}
	#if defined(LIBMAUS_HAVE_NETTLE)
	else if ( hash == "sha1" )
	{
		return bamseqchksumTemplate<SHA1SimpleSums>(arginfo);
	}
	else if ( hash == "sha224" )
	{
		return bamseqchksumTemplate<SHA2_224_SimpleSums>(arginfo);
	}
	else if ( hash == "sha256" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_SimpleSums>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_SimpleSums>(arginfo);
	}
	else if ( hash == "sha384" )
	{
		return bamseqchksumTemplate<SHA2_384_SimpleSums>(arginfo);
	}
	else if ( hash == "sha512" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_SimpleSums>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_SimpleSums>(arginfo);
	}
	#endif
	else if ( hash == "crc32prime32" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct32>(arginfo);
	}
	else if ( hash == "crc32prime64" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct64>(arginfo);
	}
	else if ( hash == "md5prime64" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct64>(arginfo);
	}
	else if ( hash == "crc32prime96" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct96>(arginfo);
	}
	else if ( hash == "md5prime96" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct96>(arginfo);
	}
	else if ( hash == "crc32prime128" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct128>(arginfo);
	}
	else if ( hash == "md5prime128" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct128>(arginfo);
	}
	else if ( hash == "crc32prime160" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct160>(arginfo);
	}
	else if ( hash == "md5prime160" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct160>(arginfo);
	}
	else if ( hash == "crc32prime192" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct192>(arginfo);
	}
	else if ( hash == "md5prime192" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct192>(arginfo);
	}
	else if ( hash == "crc32prime224" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct224>(arginfo);
	}
	else if ( hash == "md5prime224" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct224>(arginfo);
	}
	else if ( hash == "crc32prime256" )
	{
		return bamseqchksumTemplate<CRC32PrimeProduct256>(arginfo);
	}
	else if ( hash == "md5prime256" )
	{
		return bamseqchksumTemplate<MD5PrimeProduct256>(arginfo);
	}
	#if defined(LIBMAUS_HAVE_NETTLE)
	else if ( hash == "sha1prime64" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha224prime64" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha256prime64" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			//std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct64>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha384prime64" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha512prime64" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			//std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct64>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha1prime96" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha224prime96" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha256prime96" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			//std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct96>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha384prime96" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha512prime96" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct96>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha1prime128" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha224prime128" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha256prime128" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct128>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha384prime128" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha512prime128" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct128>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha1prime160" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha224prime160" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha256prime160" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct160>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha384prime160" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha512prime160" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct160>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha1prime192" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha224prime192" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha256prime192" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct192>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha384prime192" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha512prime192" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct192>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha1prime224" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha224prime224" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha256prime224" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct224>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha384prime224" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha512prime224" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct224>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha1prime256" )
	{
		return bamseqchksumTemplate<SHA1PrimeProduct256>(arginfo);
	}
	else if ( hash == "sha224prime256" )
	{
		return bamseqchksumTemplate<SHA2_224_PrimeProduct256>(arginfo);
	}
	else if ( hash == "sha256prime256" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_256" << std::endl;
			return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct256>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_256_PrimeProduct256>(arginfo);
	}
	else if ( hash == "sha384prime256" )
	{
		return bamseqchksumTemplate<SHA2_384_PrimeProduct256>(arginfo);
	}
	else if ( hash == "sha512prime256" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct256>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeProduct256>(arginfo);
	}
	#endif
	else if ( hash == "null" )
	{
		return bamseqchksumTemplate<NullChecksums>(arginfo);
	}
	#if defined(LIBMAUS_HAVE_NETTLE)
	else if ( hash == "sha512primesums" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeSums>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeSums>(arginfo);
	}
	else if ( hash == "sha512primesums512" )
	{
		#if defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386)	&& defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
		if ( libmaus::util::I386CacheLineSize::hasSSE41() )
		{
			// std::cerr << "[V] running sse4 SHA2_512" << std::endl;
			return bamseqchksumTemplate<SHA2_512_sse4_PrimeSums512>(arginfo);
		}
		else
		#endif
			return bamseqchksumTemplate<SHA2_512_PrimeSums512>(arginfo);
	}
	#endif
	#if (! defined(NETTLE)) && defined(LIBMAUS_USE_ASSEMBLY) && defined(LIBMAUS_HAVE_i386) && defined(LIBMAUS_HAVE_SHA2_ASSEMBLY)
	else if ( hash == "sha256" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_SimpleSums>(arginfo);	
	}
	else if ( hash == "sha512" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_SimpleSums>(arginfo);	
	}
	else if ( hash == "sha256prime64" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha256prime96" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha256prime128" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha256prime160" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha256prime192" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha256prime224" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha256prime256" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_256_sse4_PrimeProduct256>(arginfo);
	}
	else if ( hash == "sha512prime64" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct64>(arginfo);
	}
	else if ( hash == "sha512prime96" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct96>(arginfo);
	}
	else if ( hash == "sha512prime128" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct128>(arginfo);
	}
	else if ( hash == "sha512prime160" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct160>(arginfo);
	}
	else if ( hash == "sha512prime192" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct192>(arginfo);
	}
	else if ( hash == "sha512prime224" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct224>(arginfo);
	}
	else if ( hash == "sha512prime256" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeProduct256>(arginfo);
	}
	else if ( hash == "sha512primesums" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeSums>(arginfo);	
	}
	else if ( hash == "sha512primesums512" && libmaus::util::I386CacheLineSize::hasSSE41() )
	{
		return bamseqchksumTemplate<SHA2_512_sse4_PrimeSums512>(arginfo);
	}
	#endif
	else
	{
		libmaus::exception::LibMausException lme;
		lme.getStream() << "unsupported hash type " << hash << std::endl;
		lme.finish();
		throw lme;
	}
}

#if defined(BIOBAMBAM_HAVE_GMP)
template<size_t k>
static libmaus::math::UnsignedInteger<k> convertNumber(mpz_t const & gmpnum)
{
	size_t const numbitsperel = 8 * sizeof(uint32_t);
	size_t const numwords = (mpz_sizeinbase(gmpnum,2) + numbitsperel - 1) / numbitsperel;
	libmaus::autoarray::AutoArray<uint32_t> A(numwords,false);
	size_t countp = numwords;
	mpz_export(A.begin(),&countp,-1,sizeof(uint32_t),0,0,gmpnum);
	libmaus::math::UnsignedInteger<k> U;
	for ( size_t i = 0; i < std::min(countp,static_cast<size_t>(k)); ++i )
		U[i] = A[i];
	return U;
}
// search for next prime number larger than 2^(32*k) using probabilistic algorithm in gmp library
template<size_t k>
libmaus::math::UnsignedInteger<k+1> nextPrime()
{
	libmaus::math::UnsignedInteger<k+1> U(1);
	U <<= (k*32);
	
	mpz_t gmpU;
	mpz_init(gmpU);
	mpz_import(gmpU,k+1,-1 /* least sign first */,sizeof(uint32_t),0 /* native endianess */,0 /* don't skip bits */, U.getWords());
	
	while ( mpz_probab_prime_p(gmpU,200) == 0 )
		mpz_add_ui(gmpU,gmpU,1);
	
	libmaus::math::UnsignedInteger<k+1> const R = convertNumber<k+1>(gmpU);
	
	mpz_clear(gmpU);
	
	return R;
}
#endif


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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif

				V.push_back ( std::pair<std::string,std::string> ( std::string("hash=<[")+getDefaultHash()+"]>", "hash digest function: " + getSupportedHashVariantsList()) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
								
				return EXIT_SUCCESS;
			}
			
		return bamseqchksum(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

