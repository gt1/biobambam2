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

#include <libmaus/aio/CheckedOutputStream.hpp>

#include <libmaus/bambam/AdapterFilter.hpp>
#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/fastx/acgtnMap.hpp>
#include <libmaus/rank/popcnt.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/Histogram.hpp>

#include <biobambam/Licensing.hpp>
#include <biobambam/Split12.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1024*1024; }
static uint64_t getDefaultRankSplit() { return 0; }

int bam12auxmerge(::libmaus::util::ArgInfo const & arginfo)
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
	
	std::string const prefilename = arginfo.getRestArg<std::string>(0);
	libmaus::bambam::BamDecoder bampredec(prefilename);

	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const ranksplit = arginfo.getValue<int>("ranksplit",getDefaultRankSplit());
	uint64_t const mod = arginfo.getValue<int>("mod",getDefaultMod());
	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;
	
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

	libmaus::autoarray::AutoArray<char> Aread;

	::libmaus::bambam::BamDecoder bamdec(std::cin,false);
	::libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	::libmaus::bambam::BamHeader const & preheader = bampredec.getHeader();

	std::string const headertext(header.text);
	std::string const preheadertext(preheader.text);
	
	libmaus::bambam::ProgramHeaderLineSet headerlines(headertext);
	libmaus::bambam::ProgramHeaderLineSet preheaderlines(preheadertext);
	
	std::vector<libmaus::bambam::HeaderLine> allheaderlines = libmaus::bambam::HeaderLine::extractLines(headertext);
	
	std::string const lastid = preheaderlines.getLastIdInChain();
	
	std::stack < std::pair<uint64_t,std::string> > pgtodo;
	for ( uint64_t i = 0; i < headerlines.roots.size(); ++i )
		pgtodo.push(std::pair<uint64_t,std::string>(headerlines.roots[i],lastid));
	
	std::string upheadtext = preheadertext;
	while ( pgtodo.size() )
	{
		uint64_t const hid = pgtodo.top().first;
		std::string const PP = pgtodo.top().second;
		pgtodo.pop();
		libmaus::bambam::HeaderLine const & line = headerlines.lines[hid];
		
		// ID, PP, PN, CL, VN
		std::string       ID = (line.M.find("ID") != line.M.end()) ? line.M.find("ID")->second : "";
		std::string const PN = (line.M.find("PN") != line.M.end()) ? line.M.find("PN")->second : "";
		std::string const CL = (line.M.find("CL") != line.M.end()) ? line.M.find("CL")->second : "";
		std::string const VN = (line.M.find("VN") != line.M.end()) ? line.M.find("VN")->second : "";
		
		upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLineRef(
			upheadtext,
			ID,
			PN,
			CL,
			PP,
			VN
		);
	
		if ( headerlines.edges.find(hid) != headerlines.edges.end() )
		{
			std::vector<uint64_t> const & children = headerlines.edges.find(hid)->second;
			
			for ( uint64_t j = 0; j < children.size(); ++j )
				pgtodo.push(std::pair<uint64_t,std::string>(children[j],ID));
		}
	}
	
	/* copy SQ lines */
	std::ostringstream sqconcstr;
	sqconcstr << upheadtext;
	for ( uint64_t i = 0; i < allheaderlines.size(); ++i )
		if ( allheaderlines[i].type == "SQ" )
			sqconcstr << allheaderlines[i].line << "\n";
	upheadtext = sqconcstr.str();

	::libmaus::bambam::BamHeader uphead(upheadtext);
	uphead.changeSortOrder("unknown");
	
	::libmaus::bambam::BamWriter writer(std::cout,uphead);
	
	::libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	::libmaus::bambam::BamAlignment & prealgn = bampredec.getAlignment();
	int64_t curid = -1;
	
	libmaus::autoarray::AutoArray< std::pair<uint8_t,uint8_t> > auxpre;
	libmaus::autoarray::AutoArray< std::pair<uint8_t,uint8_t> > auxnew;
	
	libmaus::bambam::BamAuxFilterVector auxfilter;
	
	// loop over aligned BAM file
	while ( bamdec.readAlignment() )
	{
		if ( ranksplit )
			split12(algn);
	
		// extract rank
		char const * name = algn.getName();
		char const * u1 = name;
		bool ok = true;
		uint64_t rank = 0;
		while ( *u1 && *u1 != '_' )
		{
			rank *= 10;
			rank += (*u1-'0');
			ok = ok && isdigit(*u1);
			++u1;
		}
		
		// unable to find rank?	write out as is and continue
		if ( ! ok )
		{
			algn.serialise(writer.getStream());
			continue;
		}
		
		// loop over unaligned BAM file
		while ( curid != static_cast<int64_t>(rank) )
		{
			bool const a_ok = bampredec.readAlignment();
			
			if ( ! a_ok )
			{
				libmaus::exception::LibMausException se;
				se.getStream() << "Found unexpected EOF on file " << prefilename << std::endl;
				se.finish();
				throw se;
			}
			assert ( a_ok );
			++curid;
			
			if ( verbose && (! (curid & bmask)) )
				std::cerr << "[V] " << (curid / bmod) << std::endl;
		}
			
		if ( verbose > 1 )
			std::cerr << "Merging:\n" << algn.formatAlignment(header) << "\n" << prealgn.formatAlignment(preheader) << std::endl;
		
		uint64_t pretagnum = prealgn.enumerateAuxTags(auxpre);
		uint64_t newtagnum = algn.enumerateAuxTags(auxnew);
		
		std::sort(auxpre.begin(),auxpre.begin()+pretagnum);
		std::sort(auxnew.begin(),auxnew.begin()+newtagnum);
		
		if ( verbose > 1 )
			std::cerr << "pretagnum=" << pretagnum << " newtagnum=" << newtagnum << std::endl;
		
		std::pair<uint8_t,uint8_t> * prec = auxpre.begin();
		std::pair<uint8_t,uint8_t> * pree = prec + pretagnum;
		std::pair<uint8_t,uint8_t> * preo = prec;

		std::pair<uint8_t,uint8_t> * newc = auxnew.begin();
		std::pair<uint8_t,uint8_t> * newe = newc + newtagnum;
		std::pair<uint8_t,uint8_t> * newo = newc;
		
		while ( prec != pree && newc != newe )
		{
			// pre which is not in new
			if ( *prec < *newc )
			{
				*(preo++) = *(prec++);
			}
			// tag in both, drop pre
			else if ( *prec == *newc )
			{
				*(newo++) = *(newc++);
				prec++;
			}
			// new not in pre
			else
			{
				*(newo++) = *(newc++);					
			}
		}
		
		while ( prec != pree )
			*(preo++) = *(prec++);
		while ( newc != newe )
			*(newo++) = *(newc++);
			
		pretagnum = preo-auxpre.begin();
		newtagnum = newo-auxnew.begin();
		
		for ( uint64_t i = 0; i < pretagnum; ++i )
			auxfilter.set(auxpre[i].first,auxpre[i].second);
		
		algn.copyAuxTags(prealgn, auxfilter);

		for ( uint64_t i = 0; i < pretagnum; ++i )
			auxfilter.clear(auxpre[i].first,auxpre[i].second);

		if ( verbose > 1 )
		{
			std::cerr << "pretagnum=" << pretagnum << " newtagnum=" << newtagnum << std::endl;	
			std::cerr << "result: " << algn.formatAlignment(header) << std::endl;
		}
		
		algn.serialise(writer.getStream());
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
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for every mod'th alignment if verbose" ) );
				V.push_back ( std::pair<std::string,std::string> ( "ranksplit=<["+::biobambam::Licensing::formatNumber(getDefaultRankSplit())+"]>", "split ranks (see bam12split command)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bam12auxmerge(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

