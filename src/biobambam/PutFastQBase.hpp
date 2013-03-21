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


#if ! defined(PUTFASTQBASE_HPP)
#define PUTFASTQBASE_HPP

#include <biobambam/BamBamConfig.hpp>

#if defined(BAMBAM_HAVE_SAMTOOLS)

#if defined(HAVE_BAM_H)
#include <bam.h>
#elif defined(HAVE_SAMTOOLS_BAM_H)
#include <samtools/bam.h>
#else
#error "Required bam.h header not available."
#endif

#endif

#include <libmaus/util/shared_ptr.hpp>
#include <libmaus/util/unique_ptr.hpp>
#include <libmaus/autoarray/AutoArray.hpp>
#include <libmaus/fastx/SpaceTable.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>
#include <libmaus/bambam/BamAlignment.hpp>

namespace biobambam
{
	struct PutFastQBase
	{
		typedef uint32_t bam_flag_type;

		struct CompactFastQEntry
		{
			typedef CompactFastQEntry this_type;
			typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
			typedef ::libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

			::libmaus::autoarray::AutoArray<uint8_t> T;
			unsigned int namepos;
			unsigned int seqpos;
			unsigned int qualpos;
			unsigned int qnamelen;
			unsigned int seqlen;
			unsigned int entrylen;
			bam_flag_type flags;
			
			CompactFastQEntry() : T(), namepos(1), seqpos(0), qualpos(0), qnamelen(0), seqlen(0), entrylen(0), flags(0)
			{
			
			}
			
			void checkSize(uint64_t const n)
			{
				while ( n > T.size() )
					T = ::libmaus::autoarray::AutoArray<uint8_t>( std::max(static_cast<uint64_t>(1ull),static_cast<uint64_t>(2*T.size())) );
			}
		};

		/**
		 * put fastq @ line
		 **/
		template<typename opc_type>	
		static inline void putAtLine(
			char const * qname, unsigned int const qnamelen, bam_flag_type const & flags, opc_type & opc,
			::libmaus::fastx::SpaceTable const & ST
		)
		{
			*(opc++) = '@';
			char const * qnamee = qname+qnamelen;
			
			// paired? add /1 or /2 before first space or at end of line
			if ( flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED )
			{
				while ( qname != qnamee && !ST.spacetable[*qname] )
					*(opc++) = *(qname++);

				*(opc++) = '/';
				if ( flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 )
					*(opc++) = '1';
				else
					*(opc++) = '2';

				while ( qname != qnamee )
					*(opc++) = *(qname++);
			}
			else
			{
				while ( qname != qnamee )
					*(opc++) = *(qname++);
			}

			
			*(opc++) = '\n';
		}
		
		/**
		 * put fastq + line
		 **/
		template<typename opc_type>	
		static inline void putPlusLine(
			#if defined(FULLPLUS)
			char const * qname, 
			unsigned int const qnamelen, 
			#else
			char const *,
			unsigned int const,
			#endif
			opc_type & opc)
		{
			// plus line
			*(opc++) = '+';
			#if defined(FULLPLUS)
			char const * qnamee = qname+qnamelen;
			while ( qname != qnamee )
				*(opc++) = *(qname++);
			#endif
			*(opc++) = '\n';		
		}

		#if defined(BAMBAM_HAVE_SAMTOOLS)
		/**
		 * pust fastq sequence
		 **/
		template<typename opc_type>	
		static inline void putSequence(uint8_t const * seq, unsigned int const seqlen, bam_flag_type const & flags, opc_type & opc)
		{
			// character mapping table
			static uint8_t const T[16] = { 
				4 /* 0 */, 0 /* 1 */, 1 /* 2 */, 4 /* 3 */,
				2 /* 4 */, 4 /* 5 */, 4 /* 6 */, 4 /* 7 */,
				3 /* 8 */, 4 /* 9 */, 4 /* 10 */, 4 /* 11 */,
				4 /* 12 */, 4 /* 13 */, 4 /* 14 */, 4 /* 15 */
			};
			// reverse complement
			static uint8_t const I[5] = { 3, 2, 1, 0, 4 };
			// remap
			static uint8_t const R[5] = { 'A', 'C', 'G', 'T', 'N' };
			// character mapping table
			static uint8_t const F[16] = { 
				R[T[0]] /* 0 */, R[T[1]] /* 1 */, R[T[2]] /* 2 */, R[T[3]] /* 3 */,
				R[T[4]] /* 4 */, R[T[5]] /* 5 */, R[T[6]] /* 6 */, R[T[7]] /* 7 */,
				R[T[8]] /* 8 */, R[T[9]] /* 9 */, R[T[10]] /* 10 */, R[T[11]] /* 11 */,
				R[T[12]] /* 12 */, R[T[13]] /* 13 */, R[T[14]] /* 14 */, R[T[15]] /* 15 */
			};
			// character mapping table
			static uint8_t const C[16] = { 
				R[I[T[0]]] /* 0 */, R[I[T[1]]] /* 1 */, R[I[T[2]]] /* 2 */, R[I[T[3]]] /* 3 */,
				R[I[T[4]]] /* 4 */, R[I[T[5]]] /* 5 */, R[I[T[6]]] /* 6 */, R[I[T[7]]] /* 7 */,
				R[I[T[8]]] /* 8 */, R[I[T[9]]] /* 9 */, R[I[T[10]]] /* 10 */, R[I[T[11]]] /* 11 */,
				R[I[T[12]]] /* 12 */, R[I[T[13]]] /* 13 */, R[I[T[14]]] /* 14 */, R[I[T[15]]] /* 15 */
			};

			// reverse complement?
			if ( flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE )
			{
				unsigned int i = seqlen;
				while ( i-- )
					*(opc++) = C[bam1_seqi(seq,i)];
			}
			else
			{
				for ( unsigned int i = 0; i < seqlen; ++i )
					*(opc++) = F[bam1_seqi(seq,i)];
			}

			// newline
			*(opc++) = '\n';
		
		}

		template<typename opc_type>
		static inline void putQuality(uint8_t const * qual, unsigned int const seqlen, bam_flag_type const & flags, opc_type & opc)
		{
			if ( seqlen )
			{
				if ( qual[0] == 0xFF )
				{
					for ( uint64_t i = 0; i < seqlen; ++i )				
						*(opc++) = '*';
				}
				else
				{
					if ( flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE )
					{
						uint8_t const * qualp = qual + seqlen;
						
						while ( qualp != qual )
							*(opc++) = (*(--qualp)) + 33;
					}
					else
					{
						uint8_t const * const quale = qual+seqlen;

						while ( qual != quale )
							*(opc++) = *(qual++) + 33;
					}
				}
			}
			*(opc++) = '\n';
		}

		template<typename opc_type>
		static inline void putRead(
			char const * qname,
			unsigned int const qnamelen,
			uint8_t const * seq,
			uint8_t const * qual, 
			unsigned int const seqlen, 
			bam_flag_type const & flags, 
			opc_type & opc,
			::libmaus::fastx::SpaceTable const & ST
			)
		{
			// at line
			putAtLine(qname,qnamelen,flags,opc,ST);
			// sequence
			putSequence(seq,seqlen,flags,opc);
			// plus line
			putPlusLine(qname,qnamelen,opc);
			// quality
			putQuality(qual,seqlen,flags,opc);
		}

		template<typename opc_type>
		static inline void put(bam1_t const * alignment, opc_type & opc, ::libmaus::fastx::SpaceTable const & ST)
		{
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
			
			// put reads in buffer
			putRead(qname,qnamelen,seq,bam1_qual(alignment),seqlen,flags,opc,ST);
		}

		static void putCompact(bam1_t const * alignment, CompactFastQEntry & CFQ, ::libmaus::fastx::SpaceTable const & ST)
		{
			bam1_core_t const * alignment_core = &(alignment->core);
			// length of query string
			unsigned int const seqlen = alignment_core->l_qseq;
			// flags
			bam_flag_type const flags = alignment_core->flag;
			// name of query
			char const * qname = bam1_qname(alignment);
			unsigned int const qnamelen = strlen(qname);

			unsigned int const entrylen = getFastqEntryLength(qnamelen,seqlen,flags);
			CFQ.checkSize(entrylen);
			uint8_t * T = CFQ.T.get();
			put(alignment,T,ST);

			CFQ.namepos = 1;
			CFQ.seqpos = getFastqNameLineLength(qnamelen,flags);
			CFQ.qualpos = getFastqNameLineLength(qnamelen,flags) + getFastqSeqLineLength(seqlen) + getFastqPlusLineLength(qnamelen);
			CFQ.qnamelen = qnamelen;
			CFQ.seqlen = seqlen;
			CFQ.flags = flags;
			CFQ.entrylen = entrylen;
		}
		#endif

		static uint64_t getFastqNameLineLength(
			unsigned int const qnamelen,
			bam_flag_type const flags
			)
		{
			return 1 /* @ */ + qnamelen + (( flags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED ) ? 2 : 0) /* /[12] */ + 1 /* \n */;
		}
		static uint64_t getFastqSeqLineLength(
			unsigned int const seqlen
			)
		{
			return seqlen + 1;
		}
		static uint64_t getFastqPlusLineLength(
			#if defined(FULLPLUS)
			unsigned int const qnamelen
			#else
			unsigned int const
			#endif
			)
		{
			#if defined(FULLPLUS)
			return 1 /* + */ + qnamelen + 1 /* \n */;
			#else
			return 1 /* + */ + 1 /* \n */;
			#endif
		}
		static uint64_t getFastqQualLineLength(
			unsigned int const seqlen
			)
		{
			return seqlen + 1;
		}

		static uint64_t getFastqEntryLength(
			unsigned int const qnamelen,
			unsigned int const seqlen,
			bam_flag_type const flags
		)
		{
			return
				getFastqNameLineLength(qnamelen,flags) +
				getFastqSeqLineLength(seqlen) +
				getFastqPlusLineLength(qnamelen) +
				getFastqQualLineLength(seqlen);
		}

		static void putCompact(::libmaus::bambam::BamAlignment const & alignment, CompactFastQEntry & CFQ, ::libmaus::fastx::SpaceTable const & ST)
		{
			// length of query sequence
			unsigned int const seqlen = alignment.getLseq();
			// flags
			bam_flag_type const flags = alignment.getFlags();			
			// name of query
			char const * qname = alignment.getName();
			// length of query name (-1 to remove terminating null byte)
			unsigned int const qnamelen = alignment.getLReadName()-1;

			unsigned int const entrylen = getFastqEntryLength(qnamelen,seqlen,flags);
			CFQ.checkSize(entrylen);
			uint8_t * T = CFQ.T.get();

			// at line
                        putAtLine(qname,qnamelen,flags,T,ST);
                        // put read
                        if ( alignment.isReverse() )
                        	T = alignment.decodeReadRCIt(T,seqlen);
			else
				T = alignment.decodeRead(T,seqlen);
			*(T++) = '\n';
                        // put plus line
			putPlusLine(qname,qnamelen,T);
			// quality
			if ( alignment.isReverse() )
				T = alignment.decodeQualRcIt(T,seqlen);
			else
				T = alignment.decodeQual(T,seqlen);
			*(T++) = '\n';

			CFQ.namepos = 1;
			CFQ.seqpos = getFastqNameLineLength(qnamelen,flags);
			CFQ.qualpos = getFastqNameLineLength(qnamelen,flags) + getFastqSeqLineLength(seqlen) + getFastqPlusLineLength(qnamelen);
			CFQ.qnamelen = qnamelen;
			CFQ.seqlen = seqlen;
			CFQ.flags = flags;
			CFQ.entrylen = entrylen;
		}

	};
}
#endif
