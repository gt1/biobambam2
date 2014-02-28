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
//#include <queue>
#include <zlib.h>



#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamHeader.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>


#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GetObject.hpp>
#include <libmaus/util/PutObject.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 1; }


int bamseqchksum(::libmaus::util::ArgInfo const & arginfo)
{

	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	
	::libmaus::bambam::BamDecoder dec(std::cin,false);
//	::libmaus::bambam::BamHeader const & header = dec.getHeader();

//	std::string const headertext(header.text);



//	::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type LHTsnofailure;
	
//	std::deque<std::string> qreadgroups = ::libmaus::util::stringFunctions::tokenize(readgroups,std::string(","));
//	std::vector<std::string> vreadgroups = std::vector<std::string>(qreadgroups.begin(),qreadgroups.end());
//	::libmaus::trie::Trie<char> trienofailure;
//	trienofailure.insertContainer(vreadgroups);
//	::libmaus::trie::LinearHashTrie<char,uint32_t>::unique_ptr_type LHTnofailure(trienofailure.toLinearHashTrie<uint32_t>());
//	LHTsnofailure = ::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type(LHTnofailure.release());

	libmaus::bambam::BamAlignment & algn = dec.getAlignment();

	const static uint64_t MERSENNE31 = 0x7FFFFFFFull;
	struct OrderIndependentSeqDataChecksums {
		private:
		::libmaus::autoarray::AutoArray<char> A;
		::libmaus::autoarray::AutoArray<char> B; //separate A & B can allow reordering speed up?
                static void munge_chksum (uint32_t & chksum) {
			//alter chksum value used for multiplication in a finite field
			// i.e. 0 < chksum < (2^31 -1)
			chksum &= MERSENNE31;
			if (!chksum || chksum == MERSENNE31 ) chksum = 1;
		}
		public:
		uint64_t count;
        	uint64_t name_b_seq_qual;
        	uint64_t b_seq_qual;
        	uint64_t name_b_seq;
        	uint64_t b_seq;
		OrderIndependentSeqDataChecksums() : A(), B() {
			count = 0;
        		name_b_seq_qual = 1;
        		b_seq_qual = 1;
        		name_b_seq = 1;
        		b_seq = 1;
		};
		void push(libmaus::bambam::BamAlignment const & algn) {
			if ( ! (algn.getFlags() & (::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY | ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY)) ) {
				++count;
				uint8_t flags = ( ( algn.getFlags() & ( //following flags are in the least significant byte!
					::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
					::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 | 
					::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2 ) ) >> 0 )  & 0xFF;
				uint32_t chksum = crc32(0,0,0);
				chksum = crc32(chksum,(const unsigned char *)algn.getName(),algn.getLReadName());
				chksum = crc32(chksum,(const unsigned char*) &flags, 1);
				uint64_t const len = algn.isReverse() ? algn.decodeReadRC(A) : algn.decodeRead(A);
				chksum = crc32(chksum,(const unsigned char *) A.begin(), len);
				uint32_t chksum_copy = chksum;
				munge_chksum(chksum_copy);
				name_b_seq = ( name_b_seq * chksum_copy ) % MERSENNE31;
				uint64_t const len2 = algn.isReverse() ? algn.decodeQualRC(B) : algn.decodeQual(B);
				chksum = crc32(chksum,(const unsigned char *) B.begin(), len2);
				munge_chksum(chksum);
				name_b_seq_qual = ( name_b_seq_qual * chksum ) % MERSENNE31;
			}
		};
	};

	OrderIndependentSeqDataChecksums chksums;
	
		uint64_t c = 0;
		while ( dec.readAlignment() )
		{
			chksums.push(algn);
			if ( verbose && (++c & (1024*1024-1)) == 0 ) {
				std::cerr << "[V] " << c/(1024*1024) << " " << chksums.count << " " << algn.getName() << " " << algn.isRead1() << " " << algn.isRead2() << " " << ( algn.isReverse() ? algn.getReadRC() : algn.getRead() ) << " " << ( algn.isReverse() ? algn.getQualRC() : algn.getQual() ) << " " << std::hex << (0x0 + algn.getFlags()) << std::dec << " " << chksums.name_b_seq << " " << chksums.name_b_seq_qual << " " << std::endl;
			}
		}

	std::cout << std::dec << chksums.count << " " << std::hex << " " << chksums.name_b_seq << " " << chksums.name_b_seq_qual << std::endl;
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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );

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

