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
#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <config.h>

#include <libmaus/fastx/StreamFastAReader.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>

static unsigned int getDefaultCols()
{
	return 80;
}

std::string stripName(std::string const & s)
{
	uint64_t i = 0;
	while ( i < s.size() && !isspace(s[i]) )
		++i;
		
	return s.substr(0,i);
}

void normalisefasta(libmaus::util::ArgInfo const & arginfo)
{
	libmaus::fastx::StreamFastAReaderWrapper in(std::cin);
	libmaus::fastx::StreamFastAReaderWrapper::pattern_type pattern;
	unsigned int const cols = arginfo.getValue<unsigned int>("cols",getDefaultCols());
	uint64_t offset = 0;
	
	while ( in.getNextPatternUnlocked(pattern) )
	{
		std::string const name = pattern.getStringId();
		std::string const shortname = stripName(name);
	
		std::cerr <<  shortname << "\t" << pattern.patlen << "\t" 
			<< offset+pattern.getStringId().size()+2 << "\t" << cols << "\t" << cols+1 << std::endl;
			
		pattern.printMultiLine(std::cout,cols,offset);
	}
	
	std::cout << std::flush;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::timing::RealTimeClock rtc; rtc.start();
		
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
				std::cerr << ::biobambam::Licensing::license() << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;

				V.push_back ( std::pair<std::string,std::string> ( std::string("cols=<[")+libmaus::util::NumberSerialisation::formatNumber(getDefaultCols(),0)+"]>", "column width" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
		
			
		normalisefasta(arginfo);		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
