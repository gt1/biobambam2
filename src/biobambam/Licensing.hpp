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
#if ! defined(LICENSEMESSAGE_HPP)
#define LICENSEMESSAGE_HPP

#include <string>
#include <vector>
#include <map>
#include <ostream>

namespace biobambam
{
	struct Licensing
	{
		static std::string license();
		static std::string printLeft(std::string const & s, uint64_t const w, char const fill = ' ');
		static std::ostream & printMap(std::ostream & out, std::vector< std::pair<std::string,std::string> > const & M);
		static std::string formatNumber(int64_t const n);
		static std::string formatFloatingPoint(double const n);
	};
}
#endif
