#include <biobambam/ClipAdapters.hpp>

bool clipAdapters(
	libmaus::bambam::BamAlignment & algn,
	libmaus::autoarray::AutoArray<char> & R,
	libmaus::autoarray::AutoArray<char> & Q,
	libmaus::bambam::BamSeqEncodeTable const & seqenc,
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> & cigop,
	libmaus::bambam::BamAlignment::D_array_type & T
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
		    	if ( !reverse )
			{
				if ( algn.isMapped() )
				{
					uint32_t const numcigop = algn.getCigarOperations(cigop);

					if ( numcigop == cigop.size() )
						cigop.resize(numcigop+1);

					cigop[numcigop] = libmaus::bambam::cigar_operation(5,aclip);
					algn.replaceCigarString(cigop.begin(),numcigop+1,T);
				}

				algn.replaceSequence(seqenc,R.begin(),Q.begin(),len-aclip,T);
				algn.putAuxString("qs",std::string(R.begin()+(len-aclip),R.begin()+len));
				algn.putAuxString("qq",std::string(Q.begin()+(len-aclip),Q.begin()+len));
			}
			else
			{
				if ( algn.isMapped() )
				{
			    		uint32_t const numcigop = algn.getCigarOperations(cigop);

					if ( numcigop == cigop.size() )
					    cigop.resize(numcigop+1);


					for (uint32_t i = numcigop; i > 0; i--) {
					    cigop[i] = cigop[i - 1];
					}

					cigop[0] = libmaus::bambam::cigar_operation(5,aclip);
					algn.replaceCigarString(cigop.begin(),numcigop+1,T);
				}

				algn.replaceSequence(seqenc, (R.begin() + aclip), (Q.begin() + aclip), len - aclip, T);
				algn.putAuxString("qs", std::string(R.begin(), R.begin() + aclip));
				algn.putAuxString("qq", std::string(Q.begin(), Q.begin() + aclip));
			}
		}
	}

	return true;
}
