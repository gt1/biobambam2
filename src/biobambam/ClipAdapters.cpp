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
	
	if ( aclip )
	{
		uint64_t const len = algn.decodeRead(R);
		/* uint64_t const qlen = */ algn.decodeQual(Q);
		
		if ( len - aclip )
		{
			// FIXME: if read is on the reverse strand, clip front instead of back
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
	}

	return true;
}
