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

uint64_t resetAlignment(uint8_t * const D, uint64_t blocksize, bool const resetaux, 
	libmaus::bambam::BamAuxFilterVector const * rgfilter
)
{
	libmaus::bambam::BamAlignmentEncoderBase::putRefId(D,-1);
	libmaus::bambam::BamAlignmentEncoderBase::putPos(D,-1);
	libmaus::bambam::BamAlignmentEncoderBase::putNextRefId(D,-1);
	libmaus::bambam::BamAlignmentEncoderBase::putNextPos(D,-1);
	libmaus::bambam::BamAlignmentEncoderBase::putTlen(D,0);
	libmaus::bambam::BamAlignmentEncoderBase::putMapQ(D,0);
	
	if ( rgfilter )
		blocksize = libmaus::bambam::BamAlignmentDecoderBase::filterAux(D,blocksize,*rgfilter);
	if ( resetaux )
		blocksize = libmaus::bambam::BamAlignment::eraseAux(D);
		
	blocksize = libmaus::bambam::BamAlignment::eraseCigarString(D,blocksize);
	
	uint32_t const inflags = libmaus::bambam::BamAlignmentDecoderBase::getFlags(D);
	
	if ( libmaus::bambam::BamAlignmentDecoderBase::isReverse(inflags) )
		libmaus::bambam::BamAlignment::reverseComplementInplace(D);

	if ( libmaus::bambam::BamAlignmentDecoderBase::isPaired(inflags) )
	{
		libmaus::bambam::BamAlignmentEncoderBase::putFlags(D,
			(inflags | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPROPER_PAIR))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMREVERSE)))
		
		);	
	}
	else
	{
		libmaus::bambam::BamAlignmentEncoderBase::putFlags(D,
			(inflags | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FDUP))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPROPER_PAIR))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE))) &
			(~(static_cast<uint32_t>(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMREVERSE)))
		);	
	}
	
	return blocksize;
}

bool resetAlignment(
	libmaus::bambam::BamAlignment & algn, bool const resetaux,
	uint32_t const excludeflags,
	libmaus::bambam::BamAuxFilterVector const * rgfilter
)
{
	algn.blocksize = resetAlignment(algn.D.begin(),algn.blocksize,resetaux,rgfilter);
	
	if ( algn.getFlags() & excludeflags )
		return false;
		
	return true;
}
