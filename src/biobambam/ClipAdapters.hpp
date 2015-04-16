#if !defined(BIOBAMBAM_CLIPADAPTERS_HPP)
#define BIOBAMBAM_CLIPADAPTERS_HPP

#include <libmaus/bambam/BamAlignment.hpp>

bool clipAdapters(
	libmaus::bambam::BamAlignment & algn,
	libmaus::autoarray::AutoArray<char> & R,
	libmaus::autoarray::AutoArray<char> & Q,
	libmaus::bambam::BamSeqEncodeTable const & seqenc,
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> & cigop,
	libmaus::bambam::BamAlignment::D_array_type & T
);
#endif
