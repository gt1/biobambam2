/**
    bambam
    Copyright (C) 2013 German Tischler
    Copyright (C) 2013 Genome Research Limited

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
