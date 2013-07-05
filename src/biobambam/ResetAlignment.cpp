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
#include <biobambam/ResetAlignment.hpp>

bool resetAlignment(libmaus::bambam::BamAlignment & algn, bool const resetaux)
{
	if ( resetaux )
		algn.eraseAux();
	algn.putRefId(-1);
	algn.putPos(-1);
	algn.putNextRefId(-1);
	algn.putNextPos(-1);
	algn.putTlen(0);
	algn.eraseCigarString();
	
	if ( algn.isReverse() )
		algn.reverseComplementInplace();
	
	if ( algn.isPaired() )
		algn.putFlags(
			(algn.getFlags() |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP)
			&
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPROPER_PAIR))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMREVERSE)))
		);
	else
		algn.putFlags(
			(algn.getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPROPER_PAIR))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMREVERSE)))
		);
	
	
	if ( algn.isSecondary() )
		return false;
		
	return true;
}
