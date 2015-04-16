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
			
		return bamseqchksum(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

