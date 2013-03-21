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

#if ! defined(BAMBAMOUTPUTBUFFER_HPP)
#define BAMBAMOUTPUTBUFFER_HPP

#include <bambam/PutFastQBase.hpp>

namespace bambam
{
	struct OutputBuffer : public PutFastQBase
	{
		typedef OutputBuffer this_type;
		typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

		::libmaus::autoarray::AutoArray<uint8_t> outbuf;
		uint8_t * opa;
		uint8_t * opc;
		uint8_t * ope;
		std::ostream & out;
		uint64_t written;
		::libmaus::fastx::SpaceTable const ST;
			
		OutputBuffer(uint64_t bufsize, std::ostream & rout)
		: outbuf(bufsize,false), opa(outbuf.begin()), opc(opa), ope(outbuf.end()), out(rout), written(0), ST()
		{}
		virtual ~OutputBuffer() {}

		uint64_t freeSpace()
		{
			return ope-opc;
		}
		
		inline void prepare(uint64_t)
		{
		
		}
			
		virtual void flush()
		{
			if ( opc != opa )
			{
				out.write ( reinterpret_cast<char const *>(opa),opc-opa );
					
				if ( ! out )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "Failed to write " << (opc-opa) << " bytes." << std::endl;
					se.finish();
					throw se;
				}

				written += opc-opa;
			}

			opc = opa;			
		}

		void checkSpace(uint64_t const outlen)
		{
			// buffer overflow?
			if ( freeSpace() < outlen )
			{
				flush();
				assert ( opc == opa );
			
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
			// length of alignment in fastq file		
			uint64_t const outlen = CFQ.entrylen;
			
			checkSpace(outlen);

			uint64_t const putpos = opc-opa;
			
			for ( uint64_t i = 0; i < outlen; ++i )
				*(opc++) = CFQ.T[i];
				
			return putpos;
		}
		
		template<typename ptr_type>
		uint64_t put(ptr_type p, unsigned int const len, bool const = true /* first */)
		{
			checkSpace(len);
		
			uint64_t const putpos = opc-opa;
		
			ptr_type pe = p + len;
			while ( p != pe )
				*(opc++) = *(p++);
				
			return putpos;
		}
	};
}
#endif
