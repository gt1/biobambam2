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

static int getDefaultVerbose() { return 0; }


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
	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	libmaus::bambam::BamAlignment & algn = dec.getAlignment();

	const static uint64_t MERSENNE31 = 0x7FFFFFFFull;
	struct OrderIndependentSeqDataChecksums {
		private:
		::libmaus::autoarray::AutoArray<char> A;
		::libmaus::autoarray::AutoArray<char> B; //separate A & B can allow reordering speed up?
                static void product_munged_chksum_multiply (uint64_t & product, uint32_t chksum) {
			// alter chksum passed value ready for multiplication in a finite field
			// i.e. 0 < chksum < (2^31 -1)
			chksum &= MERSENNE31;
			if (!chksum || chksum == MERSENNE31 ) chksum = 1;
			// and multiply into existing product
			product = ( product * chksum ) % MERSENNE31;
		}
		public:
		struct Products {
			uint64_t count;
        		uint64_t b_seq;
        		uint64_t name_b_seq;
        		uint64_t b_seq_qual;
        		uint64_t name_b_seq_qual;
			Products(){
				count = 0;
        			name_b_seq_qual = 1;
        			b_seq_qual = 1;
        			name_b_seq = 1;
        			b_seq = 1;
			}
		};
		Products all;
		Products pass;
		OrderIndependentSeqDataChecksums() : A(), B(), all(), pass() { };
		void push(libmaus::bambam::BamAlignment const & algn) {
			if ( ! (algn.getFlags() & (::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY | ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY)) ) {
				++all.count;
				const uint8_t flags = ( ( algn.getFlags() & ( //following flags are in the least significant byte!
					::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
					::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 | 
					::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2 ) ) >> 0 )  & 0xFF;
				const bool is_pass = ! algn.isQCFail();
				const uint32_t CRC32_INITIAL = crc32(0L, Z_NULL, 0);
				uint32_t chksum = crc32(CRC32_INITIAL,(const unsigned char *)algn.getName(),algn.getLReadName());
				chksum = crc32(chksum,(const unsigned char*) &flags, 1);
				const uint64_t len = algn.isReverse() ? algn.decodeReadRC(A) : algn.decodeRead(A);
				chksum = crc32(chksum,(const unsigned char *) A.begin(), len);
				product_munged_chksum_multiply(all.name_b_seq, chksum);
				uint32_t chksumnn = crc32(CRC32_INITIAL,(const unsigned char*) &flags, 1);
				chksumnn = crc32(chksumnn,(const unsigned char *) A.begin(), len);
				product_munged_chksum_multiply(all.b_seq, chksumnn);
				if(is_pass) {
					++pass.count;
					product_munged_chksum_multiply(pass.name_b_seq, chksum);
					product_munged_chksum_multiply(pass.b_seq, chksumnn);
				}
				const uint64_t len2 = algn.isReverse() ? algn.decodeQualRC(B) : algn.decodeQual(B);
				chksum = crc32(chksum,(const unsigned char *) B.begin(), len2);
				product_munged_chksum_multiply(all.name_b_seq_qual, chksum);
				chksumnn = crc32(chksumnn,(const unsigned char *) B.begin(), len2);
				product_munged_chksum_multiply(all.b_seq_qual, chksumnn);
				if(is_pass) {
					product_munged_chksum_multiply(pass.name_b_seq_qual, chksum);
					product_munged_chksum_multiply(pass.b_seq_qual, chksumnn);
				}
			}
		};
	};

	OrderIndependentSeqDataChecksums chksums;
	OrderIndependentSeqDataChecksums readgroup_chksums[1 + header.getNumReadGroups()];
	
		uint64_t c = 0;
		while ( dec.readAlignment() )
		{
			chksums.push(algn);
			readgroup_chksums[algn.getReadGroupId(header)+1].push(algn);
			if ( verbose && (++c & (1024*1024-1)) == 0 ) {
				std::cerr << "[V] " << c/(1024*1024) << " " << chksums.all.count << " " << algn.getName() << " " << algn.isRead1() << " " << algn.isRead2() << " " << ( algn.isReverse() ? algn.getReadRC() : algn.getRead() ) << " " << ( algn.isReverse() ? algn.getQualRC() : algn.getQual() ) << " " << std::hex << (0x0 + algn.getFlags()) << std::dec << " " << chksums.all.b_seq << " " << chksums.all.name_b_seq_qual << " " << std::endl;
			}
		}

	std::cout << "###\tset\t" << "count" << "\t" << std::hex << "\t" << "b_seq" << "\t" << "name_b_seq" << "\t" << "b_seq_qual" << "\t" << "name_b_seq_qual" << std::dec << std::endl;
	std::cout << "all\tall\t" << chksums.all.count << "\t" << std::hex << "\t" << chksums.all.b_seq << "\t" << chksums.all.name_b_seq << "\t" << chksums.all.b_seq_qual << "\t" << chksums.all.name_b_seq_qual << std::dec << std::endl;
	std::cout << "all\tpass\t" << chksums.pass.count << "\t" << std::hex << "\t" << chksums.pass.b_seq << "\t" << chksums.pass.name_b_seq << "\t" << chksums.pass.b_seq_qual << "\t" << chksums.pass.name_b_seq_qual << std::dec << std::endl;
	if(header.getNumReadGroups()){
		for(unsigned int i=0; i<=header.getNumReadGroups(); i++){
			std::cout << (i>0 ? header.getReadGroups().at(i-1).ID : "") << "\tall\t" << readgroup_chksums[i].all.count << "\t" << std::hex << "\t" << readgroup_chksums[i].all.b_seq << "\t" << readgroup_chksums[i].all.name_b_seq << "\t" << readgroup_chksums[i].all.b_seq_qual << "\t" << readgroup_chksums[i].all.name_b_seq_qual << std::dec << std::endl;
			std::cout << (i>0 ? header.getReadGroups().at(i-1).ID : "") << "\tpass\t" << readgroup_chksums[i].pass.count << "\t" << std::hex << "\t" << readgroup_chksums[i].pass.b_seq << "\t" << readgroup_chksums[i].pass.name_b_seq << "\t" << readgroup_chksums[i].pass.b_seq_qual << "\t" << readgroup_chksums[i].pass.name_b_seq_qual << std::dec << std::endl;
		}
	}
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

