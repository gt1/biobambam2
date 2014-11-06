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
#include <libmaus/digest/CRC32.hpp>
#include <libmaus/util/ArgInfo.hpp>

#include <biobambam/Licensing.hpp>
#include <biobambam/BamBamConfig.hpp>

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

struct CRC32UpdateContext
{
	typedef libmaus::digest::CRC32 digest_type;

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
};

/**
* Product checksums calculated based on basecalls and (multi segment, first
* and last) bit info, and this combined with the query name, or the basecall
* qualities, or certain BAM auxilary fields.
**/
struct Products 
{
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
		libmaus::math::UnsignedInteger<CRC32UpdateContext::digest_type::digestlength/4> const & name_b_seq,
		libmaus::math::UnsignedInteger<CRC32UpdateContext::digest_type::digestlength/4> const & b_seq,
		libmaus::math::UnsignedInteger<CRC32UpdateContext::digest_type::digestlength/4> const & b_seq_qual,
		libmaus::math::UnsignedInteger<CRC32UpdateContext::digest_type::digestlength/4> const & b_seq_tags
	)
	{
		push_count(1);
		push_name_b_seq(name_b_seq[0]);
		push_b_seq(b_seq[0]);
		push_b_seq_qual(b_seq_qual[0]);
		push_b_seq_tags(b_seq_tags[0]);
	}
	
	public:
	Products() : count(0), b_seq(1), name_b_seq(1), b_seq_qual(1), b_seq_tags(1)
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
	
	void push(CRC32UpdateContext const & context)
	{
		push(
			context.name_flags_seq_digest,
			context.flags_seq_digest,
			context.flags_seq_qual_digest,
			context.flags_seq_tags_digest
		);
	}

	void push (Products const & subsetproducts)
	{
		count += subsetproducts.count;
		product_munged_chksum_multiply(b_seq, subsetproducts.b_seq);
		product_munged_chksum_multiply(name_b_seq, subsetproducts.name_b_seq);
		product_munged_chksum_multiply(b_seq_qual, subsetproducts.b_seq_qual);
		product_munged_chksum_multiply(b_seq_tags, subsetproducts.b_seq_tags);
	};
	bool operator== (Products const & other) const
	{
		return count==other.count &&
			b_seq==other.b_seq && name_b_seq==other.name_b_seq &&
			b_seq_qual==other.b_seq_qual && b_seq_tags==other.b_seq_tags;
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
struct OrderIndependentSeqDataChecksums {
	private:
	::libmaus::autoarray::AutoArray<char> A; // check with German: can we change/treat underlying data block to unsigned char / uint8_t?
	::libmaus::autoarray::AutoArray<char> B; //separate A & B can allow reordering speed up?
	public:
	std::vector<std::string> const auxtags;
	// to provide fast checking of aux fields when processing and validation at startup
	::libmaus::bambam::BamAuxFilterVector const auxtagsfilter;
	Products all;
	Products pass;
	OrderIndependentSeqDataChecksums() : A(), B(), auxtags(getDefaultAuxTags()), auxtagsfilter(auxtags), all(), pass() { };
	/**
	* Combine primary sequence data from alignment record into checksum products
	*
	* Ignores Supplementary and auxilary alignment records so primary data in not 
	* included twice.
	*
	* @param algn BAM alignment record
	**/
	void push(libmaus::bambam::BamAlignment const & algn, CRC32UpdateContext & context)
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
			uint8_t const flags = ( ( algn.getFlags() & ( //following flags are in the least significant byte!
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 | 
				::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2 ) ) >> 0 )  & 0xFF;
				
			context.pass = ! algn.isQCFail();
			uint64_t const len = algn.isReverse() ? algn.decodeReadRC(A) : algn.decodeRead(A);
			
			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t const CRC32_INITIAL = crc32(0L, Z_NULL, 0);
			#endif

			context.ctx_name_flags_seq.init();
			context.ctx_name_flags_seq.update(reinterpret_cast<uint8_t const *>(algn.getName()),algn.getLReadName());
			context.ctx_name_flags_seq.update(reinterpret_cast<uint8_t const *>(&flags),1);
			context.ctx_name_flags_seq.update(reinterpret_cast<uint8_t const *>(A.begin()),len);
						
			#if defined(BAM_SEQ_CHKSUM_DEBUG)
			uint32_t chksum_name_flags_seq = crc32(CRC32_INITIAL,reinterpret_cast<const unsigned char *>(algn.getName()),algn.getLReadName());
			chksum_name_flags_seq = crc32(chksum_name_flags_seq,reinterpret_cast<const unsigned char*>( &flags), 1);
			chksum_name_flags_seq = crc32(chksum_name_flags_seq,reinterpret_cast<const unsigned char *>( A.begin()), len);
			#endif
						
			// flags + sequence
			context.ctx_flags_seq.init();
			context.ctx_flags_seq.update(reinterpret_cast<uint8_t const *>(&flags),1);
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
	};
	
	void push(CRC32UpdateContext const & context)
	{
		all.push(context);
			
		if ( context.pass )
			pass.push(context);	
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
		struct ArrayErase<OrderIndependentSeqDataChecksums>
		{
			static void erase(OrderIndependentSeqDataChecksums *, uint64_t const)
			{
			
			}
		};
	}
}

int bamseqchksum(::libmaus::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();

	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	::libmaus::bambam::BamAlignment & algn = dec.getAlignment();
	OrderIndependentSeqDataChecksums chksums;
	libmaus::autoarray::AutoArray<OrderIndependentSeqDataChecksums> readgroup_chksums(1 + header.getNumReadGroups(),false);
	CRC32UpdateContext updatecontext;
	
	uint64_t c = 0;
	while ( dec.readAlignment() )
	{
		chksums.push(algn,updatecontext);
		readgroup_chksums[algn.getReadGroupId(header)+1].push(updatecontext);
		
		if ( verbose && (++c & (1024*1024-1)) == 0 )
		{
			std::cerr << "[V] " << c/(1024*1024) << " " << chksums.all.get_count() << " " << algn.getName() << " " << algn.isRead1()
			<< " " << algn.isRead2() << " " << ( algn.isReverse() ? algn.getReadRC() : algn.getRead() ) << " "
			<< ( algn.isReverse() ? algn.getQualRC() : algn.getQual() ) << " " << std::hex << (0x0 + algn.getFlags())
			<< std::dec << " " << chksums.all.get_b_seq() << " " << chksums.all.get_b_seq_tags() << " " << std::endl;
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
		OrderIndependentSeqDataChecksums chksumschk;
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
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif

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

