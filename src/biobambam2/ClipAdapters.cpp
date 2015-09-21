#include <biobambam2/ClipAdapters.hpp>

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
		
		if ( len - aclip )
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
				
 				cigop[index] = libmaus2::bambam::cigar_operation(5,aclip + hardclip);
				
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
