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
#include <biobambam/Strip12.hpp>

bool strip12(libmaus::bambam::BamAlignment & algn)
{
	char const * name = algn.getName();
	
	char const * u1 = name;
	
	while ( *u1 && *u1 != '_' )
		++u1;
					
	if ( ! *u1 )
		return true;
	else
	{
		bool ok = true;
		uint64_t ranka = 0;
			
		for ( char const * t1 = name; t1 != u1; ++t1 )
		{	
			ranka *= 10;
			ranka += ((*t1)-'0');
			ok = ok && isdigit(*t1);
		}

		int const read1 = algn.isRead1() ? 1 : 0;
		int const read2 = algn.isRead2() ? 1 : 0;
			
		if ( (read1+read2 != 1) || (!ok) )
		{
			return true;
		}
		else
		{
			std::ostringstream upnamestr;

			upnamestr << (u1+1);

			std::string const upname = upnamestr.str();
				
			algn.replaceName(upname.begin(),upname.size());
			
			return true;
		}
	}

}
