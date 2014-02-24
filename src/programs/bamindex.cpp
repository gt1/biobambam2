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
#include <config.h>

#include <libmaus/bambam/BamIndexGenerator.hpp>
#include <libmaus/lz/BgzfInflate.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>

#include <biobambam/Licensing.hpp>

bool getDefaultVerbose() { return true; }
bool getDefaultDisableValidation() { return false; }

int bamindex(libmaus::util::ArgInfo const & arginfo, std::istream & in, std::ostream & out)
{
	bool const debug = arginfo.getValue<unsigned int>("debug",0);
	unsigned int const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	bool const validate = !(arginfo.getValue<unsigned int>("disablevalidation",getDefaultDisableValidation()));
	std::string const tmpfileprefix = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());

	libmaus::lz::BgzfInflate<std::istream> rec(in);
	
	libmaus::bambam::BamIndexGenerator BIG(tmpfileprefix,verbose,validate,debug);

	libmaus::autoarray::AutoArray<uint8_t> B(libmaus::lz::BgzfConstants::getBgzfMaxBlockSize());
	libmaus::lz::BgzfInflateInfo rinfo;
	while ( ! (rinfo=rec.readAndInfo(reinterpret_cast<char *>(B.begin()),B.size())).streameof )
		BIG.addBlock(B.begin(),rinfo.compressed,rinfo.uncompressed);

	BIG.flush(out);
	
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

				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report (default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable alignment validation (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<["+arginfo.getDefaultTmpFileName()+"]>", "temporary file prefix (default: create in current directory)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamindex(arginfo,std::cin,std::cout);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

