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
#include "config.h"

#include <iostream>
#include <queue>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamAlignmentNameComparator.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 1; }

int bamchecksort(libmaus::util::ArgInfo const & arginfo)
{
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());

	libmaus::bambam::BamDecoder bamdec(std::cin);
	libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	std::string const sortorder = libmaus::bambam::BamHeader::getSortOrderStatic(header.text);
	libmaus::bambam::BamAlignment prevalgn;
	
	if ( bamdec.readAlignment() )
	{
		prevalgn.swap(algn);
		
		if ( sortorder == "coordinate" )
		{
			uint64_t c = 0;
		
			while ( bamdec.readAlignment() )
			{
				bool const ok =
					(static_cast<uint32_t>(    algn.getRefID()) >
					 static_cast<uint32_t>(prevalgn.getRefID())
					)
					||
					(
						(static_cast<uint32_t>(    algn.getRefID()) ==
						 static_cast<uint32_t>(prevalgn.getRefID())
						)
						&&
						(static_cast<uint32_t>(    algn.getPos()) >=
						 static_cast<uint32_t>(prevalgn.getPos())
						)
					);
					
				if ( ! ok )
				{
					libmaus::exception::LibMausException se;
					se.getStream() << "Broken order:";
					se.getStream() << prevalgn.formatAlignment(header) << std::endl;
					se.getStream() <<     algn.formatAlignment(header) << std::endl;
					se.finish();
					throw se;
				}
				
				prevalgn.swap(algn);

				if ( verbose && ( ((++c) & ((1ull<<20)-1)) == 0 ) )
					std::cerr << "[V] " << c << std::endl;
			}
			
			if ( verbose )
				std::cerr << "[V] " << c << std::endl;
			
			std::cerr << "Alignments sorted by coordinate." << std::endl;
		}
		else if ( sortorder == "queryname" )
		{
			uint64_t c = 0;
			
			while ( bamdec.readAlignment() )
			{
				// bool const ok = libmaus::bambam::BamAlignmentNameComparator::compareInt(prevalgn,algn) <= 0;
				bool const ok = 
					!libmaus::bambam::BamAlignmentNameComparator::compare(algn,prevalgn);

				if ( ! ok )
				{
					libmaus::exception::LibMausException se;
					se.getStream() << "Broken order:";
					se.getStream() << prevalgn.formatAlignment(header) << std::endl;
					se.getStream() <<     algn.formatAlignment(header) << std::endl;
					se.getStream() << libmaus::bambam::BamAlignmentNameComparator::compareInt(prevalgn,algn) << std::endl;
					se.finish();
					throw se;
				}

				prevalgn.swap(algn);

				if ( verbose && ( ((++c) & ((1ull<<20)-1)) == 0 ) )
					std::cerr << "[V] " << c << std::endl;
			}

			if ( verbose )
				std::cerr << "[V] " << c << std::endl;

			std::cerr << "Alignments sorted by query name." << std::endl;
		}
		else
		{
			std::cerr << "[V] not checking order for \"" << sortorder << "\"" << std::endl;
		}
	}


	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		
		for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
			if ( 
				arginfo.restargs[i] == "-v"
				||
				arginfo.restargs[i] == "--version"
			)
			{
				std::cerr << ::biobambam::Licensing::license();
				return EXIT_SUCCESS;
			}
			else if ( 
				arginfo.restargs[i] == "-h"
				||
				arginfo.restargs[i] == "--help"
			)
			{
				std::cerr << ::biobambam::Licensing::license();
				std::cerr << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamchecksort(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

