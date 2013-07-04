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
#include <biobambam/AttachRank.hpp>

bool attachRank(libmaus::bambam::BamAlignment & algn, uint64_t const c, libmaus::bambam::BamAuxFilterVector const & zzbafv)
{
	algn.filterOutAux(zzbafv);

	uint8_t const R[8] = {
		static_cast<uint8_t>((c >> ((8-0-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-1-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-2-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-3-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-4-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-5-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-6-1)*8)) & 0xFF),
		static_cast<uint8_t>((c >> ((8-7-1)*8)) & 0xFF)
	};

	algn.putAuxNumberArray("zz", &R[0], sizeof(R)/sizeof(R[0]));

	return true;
}
