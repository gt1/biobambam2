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

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }

int bamclipreinsert(::libmaus::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	
	switch ( level )
	{
		case Z_NO_COMPRESSION:
		case Z_BEST_SPEED:
		case Z_BEST_COMPRESSION:
		case Z_DEFAULT_COMPRESSION:
			break;
		default:
		{
			::libmaus::exception::LibMausException se;
			se.getStream()
				<< "Unknown compression level, please use"
				<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
				<< " level=" << Z_BEST_SPEED << " (fast) or"
				<< " level=" << Z_BEST_COMPRESSION << " (best) or"
				<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
			se.finish();
			throw se;
		}
			break;
	}

	::libmaus::bambam::BamDecoder dec(std::cin,false);
	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamclipreinsert", // ID
		"bamclipreinsert", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
		
	// construct new header
	libmaus::bambam::BamHeader const uphead(upheadtext);
 	libmaus::bambam::BamWriter writer(std::cout,uphead,level);
 	libmaus::bambam::BamAuxFilterVector bafv;
 	// bafv.set('z','z');
 	// std::vector<uint8_t> R(8);
 	// std::string const zz("zz");
 	
	libmaus::bambam::BamAlignment & algn = dec.getAlignment();
	uint64_t c = 0;
	libmaus::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > auxtags;

	while ( dec.readAlignment() )
	{
		#if 0
		algn.filterOutAux(bafv);

		for ( uint64_t i = 0; i < R.size(); ++i )
			R[i] = (c >> ((8-i-1)*8)) & 0xFF;

		algn.putAuxNumberArray(zz, R);
		#endif

		uint64_t const numaux = algn.enumerateAuxTags(auxtags);
		for ( uint64_t i = 0; i < numaux; ++i )
			bafv.set(auxtags[i].first,auxtags[i].second);

		if ( (bafv('a','s') || bafv('a','h')) && bafv('q','s') && bafv('q','q') )
		{
			std::string const qs = algn.getAuxAsString("qs");
			std::string const qq = algn.getAuxAsString("qq");
			assert ( qs.size() == qq.size() );

			std::string const read = algn.getRead();
			std::string const qual = algn.getQual();
			
			std::string const cigar = algn.getCigarString();
			std::string const cigaradd = libmaus::util::NumberSerialisation::formatNumber(qs.size(),0) + "S";

			// straight
			if ( ! algn.isReverse () )
			{	
				algn.replaceSequence(read + qs, qual + qq);
				algn.replaceCigarString(cigar + cigaradd);
			}
			else
			{				
				algn.replaceSequence( libmaus::fastx::reverseComplementUnmapped(qs) + read, std::string(qq.rbegin(),qq.rend()) + qual);				
				algn.replaceCigarString(cigaradd + cigar);
			}
		}

		for ( uint64_t i = 0; i < numaux; ++i )
			bafv.clear(auxtags[i].first,auxtags[i].second);
/*
as:i:44 
aa:Z:Nextflex-PCR-2
af:f:1
ah:i:1
a3:i:94
qs:Z:AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCCAGATCTGATATCGTATGCCGTCTTCTGCTTGAACAAAAAAAATCACTACAG
qq:Z:B0BF?FDAA/0FGCCGGEFHGGC/?0?1DGGCGC..<C--<=GHCG0:C00;0;;EH.GHEB9;-..FB0;0C0.00;/...99-.://///9/
*/

		algn.serialise(writer.getStream());

		++c;
		
		if ( verbose && (c & (1024*1024-1)) == 0 )
 			std::cerr << "[V] " << c/(1024*1024) << std::endl;
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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress information" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				
				std::cerr << "The keep and remove keys are mutually exclusive. Tags are given by their two character ids. Multiple ids are separated by commas." << std::endl;
				
				return EXIT_SUCCESS;
			}
			
		return bamclipreinsert(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

