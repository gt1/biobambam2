/**
    biobambam
    Copyright (C) 2009-2014 German Tischler
    Copyright (C) 2011-2014 Genome Research Limited

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

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::bambam::BamDecoder bamdec1(std::cin);
		libmaus2::bambam::BamHeader const & header1(bamdec1.getHeader());
		libmaus2::bambam::BamAlignment const &  al1 = bamdec1.getAlignment();
		
		std::cout << header1.text;
		while ( bamdec1.readAlignment() )
			std::cout << al1.formatAlignment(header1) << "\n";		
		std::cout.flush();
		
		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
