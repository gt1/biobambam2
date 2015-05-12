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
		std::string const fn1 = arginfo.stringRestArg(0);
		std::string const fn2 = arginfo.stringRestArg(1);
		libmaus2::bambam::BamDecoder bamdec1(fn1);
		libmaus2::bambam::BamDecoder bamdec2(fn2);
		libmaus2::bambam::BamHeader const & header1(bamdec1.getHeader());
		libmaus2::bambam::BamHeader const & header2(bamdec2.getHeader());
		libmaus2::bambam::BamAlignment const &  al1 = bamdec1.getAlignment();
		libmaus2::bambam::BamAlignment const &  al2 = bamdec2.getAlignment();
		
		uint64_t lcnt = 0;
		while ( bamdec1.readAlignment() )
		{
			if ( ! bamdec2.readAlignment() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "EOF on second file\n";
				lme.finish();
				throw lme;
			}
			
			std::string const s1 = al1.formatAlignment(header1);
			std::string const s2 = al2.formatAlignment(header2);
			
			if ( s1 != s2 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Difference in line " << lcnt << "\n" << s1 << "\n" << s2 << "\n";
				lme.finish();
				throw lme;			
			}
			
			lcnt += 1;
		}
		
		if ( bamdec2.readAlignment() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "EOF on first file\n";
			lme.finish();
			throw lme;		
		}
		
		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
