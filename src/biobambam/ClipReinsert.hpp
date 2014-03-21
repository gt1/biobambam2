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
#if ! defined(BIOBAMBAM_CLIPREINSERT_HPP)
#define BIOBAMBAM_CLIPREINSERT_HPP

#include <libmaus/bambam/BamAlignment.hpp>

bool clipReinsert(
	libmaus::bambam::BamAlignment & algn,
	libmaus::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > & auxtags,
 	libmaus::bambam::BamAuxFilterVector & bafv,
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> & cigop,
	libmaus::bambam::BamAlignment::D_array_type & Tcigar,
	std::stack < libmaus::bambam::cigar_operation > & hardstack,
 	libmaus::bambam::BamAuxFilterVector const & auxfilterout
);
#endif
