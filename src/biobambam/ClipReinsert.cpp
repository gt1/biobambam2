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
#include <biobambam/ClipReinsert.hpp>

bool clipReinsert(
	libmaus::bambam::BamAlignment & algn,
	libmaus::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > & auxtags,
 	libmaus::bambam::BamAuxFilterVector & bafv,
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> & cigop,
	libmaus::bambam::BamAlignment::D_array_type & Tcigar,
	std::stack < libmaus::bambam::cigar_operation > & hardstack,
 	libmaus::bambam::BamAuxFilterVector const & auxfilterout
)
{
	uint64_t const numaux = algn.enumerateAuxTags(auxtags);
	for ( uint64_t i = 0; i < numaux; ++i )
		bafv.set(auxtags[i].first,auxtags[i].second);

	if ( (bafv('a','s') || bafv('a','h')) && bafv('q','s') && bafv('q','q') )
	{
		std::string const qs = algn.getAuxAsString("qs");
		std::string const qq = algn.getAuxAsString("qq");
		assert ( qs.size() == qq.size() );

		std::string const read = algn.getRead();
		std::string const qual = algn.getQual();

		// modify cigar string if read is mapped			
		if ( algn.isMapped() )
		{
			// get current cigar string
			uint32_t numcigop = algn.getCigarOperations(cigop);
		
			// make space for one more
			if ( cigop.size() == numcigop )
				cigop.resize(cigop.size()+1);
		
			// reverse cigar operations vector if read is mapped to reverse strand
			if ( algn.isReverse() )
				std::reverse(cigop.begin(),cigop.begin()+numcigop);
		
			// move hard clip operations to stack	
			while ( numcigop && cigop[numcigop-1].first == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP )
				hardstack.push(cigop[--numcigop]);
				
			// if last operation is soft clip, then add number of bases
			if ( numcigop && cigop[numcigop-1].first == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP )
				cigop[numcigop-1].second += qs.size();
			// otherwise add new operation
			else
				cigop[numcigop++] = libmaus::bambam::cigar_operation(
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP,qs.size()
				);
				
			// reinsert hardclip operations
			while ( hardstack.size() )
			{
				cigop[numcigop++] = hardstack.top();
				hardstack.pop();
			}

			// reverse cigar operations vector if read is mapped to reverse strand
			if ( algn.isReverse() )
				std::reverse(cigop.begin(),cigop.begin()+numcigop);
				
			algn.replaceCigarString(cigop.begin(),numcigop,Tcigar);
		}			

		// straight
		if ( ! algn.isReverse () )
			algn.replaceSequence(read + qs, qual + qq);
		else
			algn.replaceSequence( libmaus::fastx::reverseComplementUnmapped(qs) + read, std::string(qq.rbegin(),qq.rend()) + qual);
			
		algn.filterOutAux(auxfilterout);
	}

	for ( uint64_t i = 0; i < numaux; ++i )
		bafv.clear(auxtags[i].first,auxtags[i].second);

	return true;
}
