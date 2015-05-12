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
#include <biobambam2/ResetAlignment.hpp>

uint64_t resetAlignment(uint8_t * const D, uint64_t blocksize, bool const resetaux, 
	libmaus2::bambam::BamAuxFilterVector const * rgfilter
)
{
	libmaus2::bambam::BamAlignmentEncoderBase::putRefId(D,-1);
	libmaus2::bambam::BamAlignmentEncoderBase::putPos(D,-1);
	libmaus2::bambam::BamAlignmentEncoderBase::putNextRefId(D,-1);
	libmaus2::bambam::BamAlignmentEncoderBase::putNextPos(D,-1);
	libmaus2::bambam::BamAlignmentEncoderBase::putTlen(D,0);
	libmaus2::bambam::BamAlignmentEncoderBase::putMapQ(D,0);
	
	if ( rgfilter )
		blocksize = libmaus2::bambam::BamAlignmentDecoderBase::filterAux(D,blocksize,*rgfilter);
	if ( resetaux )
		blocksize = libmaus2::bambam::BamAlignment::eraseAux(D);
		
	blocksize = libmaus2::bambam::BamAlignment::eraseCigarString(D,blocksize);
	
	uint32_t const inflags = libmaus2::bambam::BamAlignmentDecoderBase::getFlags(D);
	
	if ( libmaus2::bambam::BamAlignmentDecoderBase::isReverse(inflags) )
		libmaus2::bambam::BamAlignment::reverseComplementInplace(D);

	if ( libmaus2::bambam::BamAlignmentDecoderBase::isPaired(inflags) )
	{
		libmaus2::bambam::BamAlignmentEncoderBase::putFlags(D,
			(inflags | libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FUNMAP | libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FMUNMAP) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FDUP))) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FPROPER_PAIR))) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE))) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FMREVERSE)))
		
		);	
	}
	else
	{
		libmaus2::bambam::BamAlignmentEncoderBase::putFlags(D,
			(inflags | libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FUNMAP) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FDUP))) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FPROPER_PAIR))) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE))) &
			(~(static_cast<uint32_t>(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FMREVERSE)))
		);	
	}
	
	return blocksize;
}

bool resetAlignment(
	libmaus2::bambam::BamAlignment & algn, bool const resetaux,
	uint32_t const excludeflags,
	libmaus2::bambam::BamAuxFilterVector const * rgfilter
)
{
	algn.blocksize = resetAlignment(algn.D.begin(),algn.blocksize,resetaux,rgfilter);
	
	if ( algn.getFlags() & excludeflags )
		return false;
		
	return true;
}
