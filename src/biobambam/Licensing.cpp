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
#include <biobambam/Licensing.hpp>
#include <sstream>
#include <iomanip>
#include "config.h"

std::string biobambam::Licensing::license()
{
	std::ostringstream ostr;
	ostr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
	ostr << PACKAGE_NAME << " is distributed under version 3 of the GNU General Public License." << std::endl;
	return ostr.str();
}

std::string biobambam::Licensing::printLeft(std::string const & s, uint64_t const w, char const fill)
{
	std::ostringstream ostr;
	ostr << std::setiosflags(std::ios::left);
	ostr << std::setw(w);
	ostr << std::setfill(fill);
	ostr << s;
	return ostr.str();
}

std::ostream & biobambam::Licensing::printMap(std::ostream & out, std::vector< std::pair<std::string,std::string> > const & M)
{
	uint64_t maxfield = 0;
	for ( std::vector< std::pair<std::string,std::string> >::const_iterator ita = M.begin(); ita != M.end(); ++ita )
		maxfield = std::max(maxfield,static_cast<uint64_t>(ita->first.size()));

	for ( std::vector< std::pair<std::string,std::string> >::const_iterator ita = M.begin(); ita != M.end(); ++ita )
		out << printLeft(ita->first,maxfield) << " : " << ita->second << std::endl;
	
	return out;
}

std::string biobambam::Licensing::formatNumber(int64_t const n) { std::ostringstream ostr; ostr << n; return ostr.str(); }
std::string biobambam::Licensing::formatFloatingPoint(double const n) { std::ostringstream ostr; ostr << n; return ostr.str(); }
