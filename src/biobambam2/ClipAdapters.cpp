/**
    bambam
    Copyright (C) 2013 German Tischler
    Copyright (C) 2013, 2015 Genome Research Limited

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

    ----------------------------------------------------------------------

    Cigar type code taken from htslib/sam.h under the MIT/Expat License
    as follows:

     Copyright (C) 2008, 2009, 2013-2014 Genome Research Ltd.
     Copyright (C) 2010, 2012, 2013 Broad Institute.

     Author: Heng Li <lh3@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE
**/
#include <biobambam2/ClipAdapters.hpp>

/* Start of htslib code */

#define BAM_CIGAR_TYPE  0x3C1A7
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3)

/*
    The above is a useful bit of code taken from htslib.  Their
    explanation is below.

    bam_cigar_type returns a bit flag with:
    bit 1 set if the cigar operation consumes the query
    bit 2 set if the cigar operation consumes the reference

    For reference, the unobfuscated truth table for this function is:
    BAM_CIGAR_TYPE  QUERY  REFERENCE
    --------------------------------
    BAM_CMATCH      1      1
    BAM_CINS        1      0
    BAM_CDEL        0      1
    BAM_CREF_SKIP   0      1
    BAM_CSOFT_CLIP  1      0
    BAM_CHARD_CLIP  0      0
    BAM_CPAD        0      0
    BAM_CEQUAL      1      1
    BAM_CDIFF       1      1
    BAM_CBACK       0      0
    --------------------------------

    End of htslib code.
*/


bool clipAdapters(
	libmaus2::bambam::BamAlignment & algn,
	libmaus2::autoarray::AutoArray<char> & R,
	libmaus2::autoarray::AutoArray<char> & Q,
	libmaus2::bambam::BamSeqEncodeTable const & seqenc,
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> & cigop,
	libmaus2::bambam::BamAlignment::D_array_type & T
)
{
	// a3,as
	uint64_t const asclip = algn.hasAux("as") ? algn.getAuxAsNumber<int>("as") : 0;
	uint64_t const a3clip = algn.hasAux("a3") ? algn.getAuxAsNumber<int>("a3") : 0;
	uint64_t const aclip = std::max(asclip,a3clip);
	bool     const reverse = algn.isReverse();

	if ( aclip )
	{
		uint64_t const len = algn.decodeRead(R);
		algn.decodeQual(Q);

		if ( (len - aclip) > 1 )
		{
			if ( algn.isMapped() )
			{
				uint32_t const numcigop = algn.getCigarOperations(cigop);

				if ( numcigop == cigop.size() )
					cigop.resize(numcigop+1);

				if ( reverse )
				{
				    std::reverse(cigop.begin(),cigop.begin()+numcigop);
				}


				// can't just add a HC to the cigar
				uint32_t index;
				uint32_t hardclip = 0;
				uint32_t cig_type;
				int32_t  left     = aclip;
				int32_t  repos    = 0;

				for ( index = numcigop - 1; index > 0; index-- )
				{
				    	cig_type = bam_cigar_type(cigop[index].first);

					if ( cig_type == 0 )
					{
				    	    	hardclip += cigop[index].second;
					}
					else
					{
					    	if ( cig_type & 1 )
						{
						    	if ( cigop[index].second < left )
							{
					    	    	    	left -= cigop[index].second;
						    	}
							else
							{
							    	break;
						    	}
					    	}

				    	    	if ( cig_type & 2 )
						{
				    		    	// move pos if reversed
						    	repos += cigop[index].second;
					    	}
					}
				}

				cig_type = bam_cigar_type(cigop[index].first);

				if ( cigop[index].second != left )
				{
				    	cigop[index++].second -= left;
				}

 				cigop[index] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP, aclip + hardclip);

			    	if ( numcigop > index + 1 )
				    	cigop.resize(index + 1);

				if ( reverse )
				{
				    std::reverse(cigop.begin(),cigop.begin() + index + 1);

				    // account for the last possible pos move
    	    	    	    	    if ( cig_type & 2 )
				    	    repos += left;

				    if ( repos )
				    {
				    	// clipping has moved the pos point
					algn.putPos(algn.getPos() + repos);
				    }

				}

				algn.replaceCigarString(cigop.begin(),index + 1,T);
			}


		    	if ( !reverse )
			{
				algn.replaceSequence(seqenc,R.begin(),Q.begin(),len-aclip,T);
				algn.putAuxString("qs",std::string(R.begin()+(len-aclip),R.begin()+len));
				algn.putAuxString("qq",std::string(Q.begin()+(len-aclip),Q.begin()+len));
			}
			else
			{
				algn.replaceSequence(seqenc, (R.begin() + aclip), (Q.begin() + aclip), len - aclip, T);
				algn.putAuxString("qs", std::string(R.begin(), R.begin() + aclip));
				algn.putAuxString("qq", std::string(Q.begin(), Q.begin() + aclip));
			}
		}
	}

	return true;
}
