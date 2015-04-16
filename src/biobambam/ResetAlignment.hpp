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
#if ! defined(BIOBAMBAM_RESETALIGNMENT_HPP)
#define BIOBAMBAM_RESETALIGNMENT_HPP

#include <libmaus/bambam/BamAlignment.hpp>

uint64_t resetAlignment(
	uint8_t * const D, uint64_t blocksize, bool const resetaux = true,
	libmaus::bambam::BamAuxFilterVector const * rgfilter = 0
);
bool resetAlignment(libmaus::bambam::BamAlignment & algn, bool const resetaux = true, 
	uint32_t const excludeflags = 
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
	libmaus::bambam::BamAuxFilterVector const * rgfilter = 0
);
#endif
