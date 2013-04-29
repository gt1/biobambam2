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
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/util/ArgInfo.hpp>

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		uint64_t const minmapped = arginfo.getValue<uint64_t>("minmapped",0);
		uint64_t const maxmapped = arginfo.getValue<uint64_t>("maxmapped",std::numeric_limits<uint64_t>::max());
		uint64_t const minlen = arginfo.getValue<uint64_t>("minlen",0);

		int const level = arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION);
		
		switch ( level )
		{
			case Z_NO_COMPRESSION:
			case Z_BEST_SPEED:
			case Z_BEST_COMPRESSION:
			case Z_DEFAULT_COMPRESSION:
				break;
			default:
			{
				::libmaus::exception::LibMausException se;
				se.getStream()
					<< "Unknown compression level, please use"
					<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
					<< " level=" << Z_BEST_SPEED << " (fast) or"
					<< " level=" << Z_BEST_COMPRESSION << " (best) or"
					<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
				se.finish();
				throw se;
			}
				break;
		}

		::libmaus::bambam::BamDecoder BD(std::cin);
		::libmaus::bambam::BamHeader const & bamheader = BD.bamheader;
		::libmaus::bambam::BamAlignment & alignment = BD.getAlignment();
		::libmaus::bambam::BamWriter writer(std::cout,bamheader,level);
		
		while ( BD.readAlignment() )
		{
			bool const a_1_mapped = !(alignment.getFlags() & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP);
			bool const a_2_mapped = !(alignment.getFlags() & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP);
			bool const proper     =  (alignment.getFlags() & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPROPER_PAIR);

			uint64_t const nummapped = (a_1_mapped?1:0)+(a_2_mapped?1:0)+(proper?1:0);

			if ( 
				nummapped >= minmapped && 
				nummapped <= maxmapped && 
				alignment.getLseq() >= static_cast<int64_t>(minlen)
			)
				alignment.serialise(writer.bgzfos);
		}	
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
