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

#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/fastx/acgtnMap.hpp>

#include <libmaus/rank/popcnt.hpp>

#include <libmaus/util/ArgInfo.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static unsigned int getDefaultSeedLength() { return 12; }

template<typename iterator>
bool compareHamming(
	iterator ita,
	iterator const itae,
	iterator itb,
	double const maxmismatchrate,
	unsigned int const baselength,
	unsigned int const basemismatches
)
{
	unsigned int const maxmismatches = maxmismatchrate * (baselength + (itae-ita));
	unsigned int nummis = basemismatches;
	
	while ( ita != itae && nummis <= maxmismatches )
		if ( *(ita++) != *(itb++) )
			nummis++;
			
	return (ita == itae) && (nummis <= maxmismatches);
}

int bamadapterfind(::libmaus::util::ArgInfo const & arginfo)
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
	uint64_t const seedlength = 
		std::min(
			static_cast<unsigned int>((8*sizeof(uint64_t))/3),
			std::max(arginfo.getValue<unsigned int>("seedlength",getDefaultSeedLength()),1u));
	double const mismatchrate = std::min(100u,arginfo.getValue<unsigned int>("maxmismatchpercent",4))/100.0;
	unsigned int const maxseedmismatches = arginfo.getValue<unsigned int>("maxseedmismatches",static_cast<unsigned int>(std::ceil(seedlength * mismatchrate)));
	
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

	::libmaus::bambam::BamDecoder bamdec(std::cin,false);
	::libmaus::bambam::BamHeader const & header = bamdec.getHeader();

	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamadapterfind", // ID
		"bamadapterfind", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);
	::libmaus::bambam::BamWriter writer(std::cout,uphead);
	
	bool running = true;
	libmaus::bambam::BamAlignment algns[2];
	libmaus::bambam::BamAlignment & inputalgn = bamdec.getAlignment();
	libmaus::autoarray::AutoArray<char> seqs[2];
	
	libmaus::autoarray::AutoArray<uint8_t> S(256,false);
	std::fill(S.begin(),S.end(),4);
	S['a'] = S['A'] = 0;
	S['c'] = S['C'] = 1;
	S['g'] = S['G'] = 2;
	S['t'] = S['T'] = 3;

	libmaus::autoarray::AutoArray<uint8_t> R(256,false);
	std::fill(R.begin(),R.end(),5);
	R['a'] = R['A'] = 0;
	R['c'] = R['C'] = 1;
	R['g'] = R['G'] = 2;
	R['t'] = R['T'] = 3;
	
	static uint64_t const mmask =
		(1ull << 0) |
		(1ull << 3) |
		(1ull << 6) |
		(1ull << 9) |
		(1ull << 12) |
		(1ull << 15) |
		(1ull << 18) |
		(1ull << 21) |
		(1ull << 24) |
		(1ull << 27) |
		(1ull << 30) |
		(1ull << 33) |
		(1ull << 36) |
		(1ull << 39) |
		(1ull << 42) |
		(1ull << 45) |
		(1ull << 48) |
		(1ull << 51) |
		(1ull << 54) |
		(1ull << 57) |
		(1ull << 60) |
		(1ull << 63);
	
	while ( (running = bamdec.readAlignment()) )
	{
		// if this is a single end read, then write it back and try the next one
		if ( ! inputalgn.isPaired() )
		{
			inputalgn.serialise(writer.getStream());
			continue;
		}
	
		// store alignment in algns[0]
		algns[0].swap(inputalgn);
		
		// read next alignment
		bool const okb = bamdec.readAlignment();
		
		if ( ! okb )
		{
			std::cerr << "[D] warning: orphan alignment"  << std::endl;
			algns[0].serialise(writer.getStream());
			break;
		}
		
		// if next is not paired or name does not match
		if ( (! inputalgn.isPaired()) || strcmp(algns[0].getName(), inputalgn.getName()) )
		{
			std::cerr << "[D] warning: orphan alignment" << std::endl;
			algns[0].serialise(writer.getStream());
			bamdec.putback();
			continue;
		}
		
		// sanity checks
		assert ( algns[0].isPaired() );
		assert ( inputalgn.isPaired() );
		assert ( strcmp(algns[0].getName(),inputalgn.getName()) == 0 );
		
		// are the read in the correct order? if not, write them out without touching them
		if ( !(algns[0].isRead1() && inputalgn.isRead2()) )
		{
			std::cerr << "[D] warning: reads are not in the correct order" << std::endl;
			algns[0].serialise(writer.getStream());
			inputalgn.serialise(writer.getStream());
			continue;
		}

		// are the reads both non empty?
		if ( !(algns[0].getLseq() && inputalgn.getLseq()) )
		{
			std::cerr << "[D] warning: empty read" << std::endl;
			algns[0].serialise(writer.getStream());
			inputalgn.serialise(writer.getStream());
			continue;
		}
		
		// put second read in algns[1]
		algns[1].swap(inputalgn);
		
		unsigned int const rev0 = algns[0].isReverse() ? 1 : 0;
		unsigned int const rev1 = algns[1].isReverse() ? 1 : 0;
		unsigned int const map0 = algns[0].isMapped() ? 1 : 0;
		unsigned int const map1 = algns[1].isMapped() ? 1 : 0;
		
		/* 
		 * if both are mapped then both are supposed to be stored relative to
		 * the forward strand;
		 * if no end is mapped, and the second read is marked as one the
		 * reverse strand, then we do the same
		 */
		if ( 
			(map0 + map1 == 2) 
			||
			( ((map0+map1)==0) && (!rev0) && (rev1) )
		)
		{
			algns[0].decodeRead(seqs[0]);
			algns[1].decodeRead(seqs[1]);
		}
		else if ( (map0+map1) == 0 && (!rev0) && (!rev1) )
		{
			algns[0].decodeRead(seqs[0]);
			algns[1].decodeReadRC(seqs[1]);	
		}
		else
		{
			algns[0].serialise(writer.getStream());
			algns[1].serialise(writer.getStream());		
			continue;
		}
			
		uint64_t const l0 = algns[0].getLseq();
		uint64_t const l1 = algns[1].getLseq();
		uint64_t const lm = std::min(l0,l1);
		/*
		 * local seed length; this may be smaller then the
		 * global seed length if one of the too reads
		 * is shorter than the seed length
		 */
		uint64_t const lseedlength = std::min(seedlength,lm);
		uint64_t const lseedmask = libmaus::math::lowbits(3*lseedlength);

		#if 0
		std::cerr << std::string(seqs[0].begin(),seqs[0].begin()+l0) << std::endl;
		std::cerr << std::string(seqs[1].begin(),seqs[1].begin()+l1) << std::endl;
		#endif
		
		// compute seed from last lseedlength bases of second read
		uint64_t seed = 0;
		uint8_t const * p = reinterpret_cast<uint8_t const *>(seqs[1].begin() + l1);
		for ( uint64_t i = 0; i < lseedlength; ++i )
		{
			seed <<= 3;
			seed  |= S[*(--p)];
		}

		// compute mask of length lseedlength-1 from back of first read
		uint64_t query = 0;
		uint8_t const * const qe = reinterpret_cast<uint8_t const *>(seqs[0].begin());
		uint8_t const * q = qe + l0;
		uint64_t matchpos = l0-lseedlength;
		
		for ( uint64_t i = 0; i < lseedlength-1; ++i )
		{
			query <<= 3;
			query |= R[*(--q)];
		}
				
		// try to find seed in first read
		do
		{
			// add next base (step backward from end)
			query <<= 3;			
			query |= R[*(--q)];
			query &= lseedmask;

			// compute number of mismatches			
			uint64_t dif = (query ^ seed);
			dif = (dif | (dif >> 1) | (dif >> 2)) & mmask;
			unsigned int const difcnt = libmaus::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(dif);
			
			#if defined(DIFCNTDEBUG)
			unsigned int debdifcnt = 0;
			for ( unsigned int i = 0; i < lseedlength; ++i )
				if ( S[*(seqs[0].begin()+matchpos+i)] != R[*(seqs[1].begin()+l1-lseedlength + i)] )
					debdifcnt++;
			assert ( debdifcnt == difcnt );
			#endif
			
			// if the seed matches, then look at the rest
			if ( difcnt <= maxseedmismatches )
			{
				uint64_t const end0 = matchpos + lseedlength;
				uint64_t const end1 = l1;
				uint64_t const seedp0 = end0 - lseedlength;
				uint64_t const seedp1 = end1 - lseedlength;
				uint64_t const restoverlap = std::min(seedp0,seedp1);

				uint8_t const * check0  = reinterpret_cast<uint8_t const *>(seqs[0].begin()+seedp0-restoverlap);
				uint8_t const * check0e = check0 + restoverlap;
				uint8_t const * check1  = reinterpret_cast<uint8_t const *>(seqs[1].begin()+seedp1-restoverlap);
				
				uint64_t const maxmis = (restoverlap+lseedlength) * mismatchrate;
				
				uint64_t nummis = difcnt;
				while ( nummis <= maxmis && check0 != check0e  )
					if ( S[*check0++] != R[*check1++] )
						nummis++;
				
				if ( check0 == check0e && nummis <= maxmis )
				{
					// adapter length for read 1
					uint64_t const al0 = (l0-end0);
					// adapter length for read 2
					uint64_t const al1 = end1-(restoverlap+lseedlength);
					// maximum number of bases to be compared
					uint64_t const almax = 12;
					// number of bases compared in this case
					uint64_t const alcmp = std::min(almax,std::min(al0,al1));
					
					char const * ap0 = seqs[0].begin()+end0;
					char const * ap0e = ap0 + alcmp;
					char const * ap1 = seqs[1].begin()+al1;
					unsigned int aldif = 0;
					
					while ( ap0 != ap0e )
						if ( S[*(ap0++)] != R[libmaus::fastx::invertUnmapped(*(--ap1))] )
							aldif++;

					if ( (ap0 == ap0e) && ! aldif  /* && aldif < 5 */ )
					{
						#if 0
						std::cerr
							<< "overlap: " << algns[0].getName() 
							<< " mismatchrate=" << 
								nummis << "/" << (restoverlap+lseedlength) << "=" <<
								static_cast<double>(nummis)/(restoverlap+lseedlength)
							<< " map0=" << map0
							<< " map1=" << map1 
							<< " al0=" << al0 
							<< " al1=" << al1
							<< " aldif=" << aldif
							<< " alcmp=" << alcmp
							<< "\n" <<
						std::string(
							seqs[0].begin()+seedp0-restoverlap,
							seqs[0].begin()+end0) << "\n" <<
						std::string(
							seqs[1].begin()+seedp1-restoverlap,
							seqs[1].begin()+end1) << "\n";
							
						std::cerr << "assumed adapter on read 1 ["<<al0<<"]: "
							<< std::string(
								seqs[0].begin()+end0,
								seqs[0].begin()+l0
							) << std::endl;
						
						std::cerr << "assumed adapter on read 2 ["<<al1<<"]: "
							<< libmaus::fastx::reverseComplementUnmapped(std::string(
								seqs[1].begin(),
								seqs[1].begin()+al1)) << std::endl;
						#endif
						
						algns[0].putAuxNumber("ah",'i',1);
						algns[1].putAuxNumber("ah",'i',1);
						algns[0].putAuxNumber("a3",'i',al0);
						algns[1].putAuxNumber("a3",'i',al1);
						// ah:i:1
						
						break;
					}
				}
			}
			
			--matchpos;
		} while ( q != qe );

		// std::cerr << "pair for " << algns[0].getName() << std::endl;
		
		algns[0].serialise(writer.getStream());		
		algns[1].serialise(writer.getStream());		
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

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamadapterfind(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

