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

#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus/util/ArgInfo.hpp>

#include <biobambam/Licensing.hpp>
#include <biobambam/KmerPoisson.hpp>

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
			
				V.push_back ( std::pair<std::string,std::string> ( "L", "length of reference sequence" ) );
				V.push_back ( std::pair<std::string,std::string> ( "n", "number of occurences" ) );
				V.push_back ( std::pair<std::string,std::string> ( "KA", "number of A bases in the query sequence" ) );
				V.push_back ( std::pair<std::string,std::string> ( "KC", "number of C bases in the query sequence" ) );
				V.push_back ( std::pair<std::string,std::string> ( "KG", "number of G bases in the query sequence" ) );
				V.push_back ( std::pair<std::string,std::string> ( "KT", "number of T bases in the query sequence" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pA", "relative frequency of A base in reference" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pC", "relative frequency of C base in reference" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pG", "relative frequency of G base in reference" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pT", "relative frequency of T base in reference" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
								
				return EXIT_SUCCESS;
			}
			
		if ( !arginfo.hasArg("L") )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Need length of reference sequence (set key L)" << std::endl;
			se.finish();
			throw se;
		}
		if ( !arginfo.hasArg("KA") )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Need number of A bases in query (set key KA)" << std::endl;
			se.finish();
			throw se;
		}
		if ( !arginfo.hasArg("KC") )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Need number of C bases in query (set key KC)" << std::endl;
			se.finish();
			throw se;
		}
		if ( !arginfo.hasArg("KG") )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Need number of G bases in query (set key KG)" << std::endl;
			se.finish();
			throw se;
		}
		if ( !arginfo.hasArg("KT") )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Need number of T bases in query (set key KT)" << std::endl;
			se.finish();
			throw se;
		}
		if ( !arginfo.hasArg("n") )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Need number of occurences (set key n)" << std::endl;
			se.finish();
			throw se;
		}
		
		uint64_t const L = arginfo.getValue<uint64_t>("L",0);
		uint64_t const KA = arginfo.getValue<uint64_t>("KA",0);
		uint64_t const KC = arginfo.getValue<uint64_t>("KC",0);
		uint64_t const KG = arginfo.getValue<uint64_t>("KG",0);
		uint64_t const KT = arginfo.getValue<uint64_t>("KT",0);
		uint64_t const n = arginfo.getValue<uint64_t>("n",0);
		double const pA = arginfo.getValue<double>("pA",0.25);
		double const pC = arginfo.getValue<double>("pC",0.25);
		double const pG = arginfo.getValue<double>("pG",0.25);
		double const pT = arginfo.getValue<double>("pT",0.25);

		double const p = pA + pC + pG + pT;
		double const eps = 1e-4;
		
		if ( p < 1.0 - eps || p > 1.0 + eps )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "pA + pC + pG + pT = " << p << " which is not within eps=" << eps << " of 1" << std::endl;
			se.finish();
			throw se;			
		}

		std::cout << kmerPoisson(L,KA,KC,KG,KT,n,pA,pC,pG,pT) << std::endl;
			
		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

