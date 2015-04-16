#if !defined(BIOBAMBAM_CLIPADAPTERS_HPP)
#define BIOBAMBAM_CLIPADAPTERS_HPP

#include <libmaus2/bambam/BamAlignment.hpp>

bool clipAdapters(
	libmaus2::bambam::BamAlignment & algn,
	libmaus2::autoarray::AutoArray<char> & R,
	libmaus2::autoarray::AutoArray<char> & Q,
	libmaus2::bambam::BamSeqEncodeTable const & seqenc,
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> & cigop,
	libmaus2::bambam::BamAlignment::D_array_type & T
);
#endif
