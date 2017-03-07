/**
    biobambam2
    Copyright (C) 2009-2017 German Tischler
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

#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <biobambam2/Licensing.hpp>

static int getDefaultVerbose() { return 1; }

int bamscrapcount(::libmaus2::util::ArgInfo const & arginfo)
{
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());

	char const * tag = "sc";

	std::map<char,uint64_t> entryCount;
	std::map<char,uint64_t> baseCount;
	uint64_t c = 0;

	uint64_t nonscrapentries = 0;
	uint64_t nonscrapbasecount = 0;

	for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
	{
		libmaus2::util::ArgInfo arg_a = arginfo;
		arg_a.replaceKey("I",arginfo.restargs.at(i));

		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper_a(libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arg_a,false /* put rank */));
		::libmaus2::bambam::BamAlignmentDecoder * ppdec_a = &(decwrapper_a->getDecoder());
		::libmaus2::bambam::BamAlignmentDecoder & dec_a = *ppdec_a;
		//::libmaus2::bambam::BamHeader const & header_a = dec_a.getHeader();
		::libmaus2::bambam::BamAlignment const & algn_a = dec_a.getAlignment();

		while ( dec_a.readAlignment() )
		{
			uint8_t const * p = ::libmaus2::bambam::BamAlignmentDecoderBase::getAux(algn_a.D.get(),algn_a.blocksize,tag);

			if (
				p
				&&
				p[0] == tag[0]
				&&
				p[1] == tag[1]
				&&
				p[2] == 'A'
			)
			{
				char const type = p[3];
				entryCount[type] += 1;
				baseCount[type] += algn_a.getLseq();
			}
			else
			{
				nonscrapentries++;
				nonscrapbasecount += algn_a.getLseq();
			}

			if ( (++c % (16*1024)) == 0 && verbose )
				std::cerr << "[V] " << c << std::endl;
		}
	}

	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

	for ( std::map<char,uint64_t>::const_iterator ita = entryCount.begin(); ita != entryCount.end(); ++ita )
	{
		assert ( baseCount.find(ita->first) != baseCount.end() );
		std::cout << "[S]\t" << ita->first << "\t" << ita->second << "\t" << baseCount.find(ita->first)->second << std::endl;
	}
	if ( nonscrapentries )
	{
		std::cout << "[N]\t" << nonscrapentries << "\t" << nonscrapbasecount << std::endl;
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo const arginfo(argc,argv);

		for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
			if (
				arginfo.restargs[i] == "-v"
				||
				arginfo.restargs[i] == "--version"
			)
			{
				std::cerr << ::biobambam2::Licensing::license();
				return EXIT_SUCCESS;
			}
			else if (
				arginfo.restargs[i] == "-h"
				||
				arginfo.restargs[i] == "--help"
			)
			{
				std::cerr << ::biobambam2::Licensing::license();
				std::cerr << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;

				std::vector< std::pair<std::string,std::string> > V;

				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;

				return EXIT_SUCCESS;
			}

		return bamscrapcount(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
