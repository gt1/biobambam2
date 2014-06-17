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
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/fastx/acgtnMap.hpp>
#include <libmaus/rank/popcnt.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/Histogram.hpp>

#include <biobambam/Licensing.hpp>
#include <biobambam/ClipAdapters.hpp>
#include <biobambam/KmerPoisson.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1024*1024; }
static uint64_t getDefaultPCT_MISMATCH() { return 10; }
static unsigned int getDefaultSEED_LENGTH() { return 12; }
static uint64_t getDefaultMIN_OVERLAP() { return 32; }
static uint64_t getDefaultADAPTER_MATCH() { return 12; }
static uint64_t getDefaultMatchMinScore() { return 16; }
static double getDefaultMatchMinFrac() { return 0.75; }
static double getDefaultMatchMinPFrac() { return 0.8; }
static int getDefaultClip() { return 0; }
static uint64_t getDefaultRefLen() { return 3000000000ull; }
static double getDefaultpA() { return 0.25; }
static double getDefaultpC() { return 0.25; }
static double getDefaultpG() { return 0.25; }
static double getDefaultpT() { return 0.25; }

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

void adapterListMatch(
	libmaus::autoarray::AutoArray<char> & Aread,
	libmaus::util::PushBuffer<libmaus::bambam::AdapterOffsetStrand> & AOSPB,
	libmaus::bambam::BamAlignment & algn,
	libmaus::bambam::AdapterFilter & AF,
	int const verbose,
	uint64_t const adpmatchminscore, // = getDefaultMatchMinScore(),
	double const adpmatchminfrac, // = getDefaultMatchMinFrac(),
	double const adpmatchminpfrac, // = getDefaultMatchMinPFrac()
	uint64_t const reflen,
	double const pA,
	double const pC,
	double const pG,
	double const pT
)
{
	uint64_t const len = algn.decodeRead(Aread);
	uint8_t const * const ua = reinterpret_cast<uint8_t const *>(Aread.begin());
			
	bool const matched = AF.searchAdapters(
		ua, len, 
		2 /* max mismatches */, 
		AOSPB,
		adpmatchminscore /* min score */,
		adpmatchminfrac /* minfrac */,
		adpmatchminpfrac /* minpfrac */
	);
			
	if ( matched )
	{
		std::sort(AOSPB.begin(),AOSPB.end(),libmaus::bambam::AdapterOffsetStrandMatchStartComparator());
		uint64_t const clip = len - AOSPB.begin()[0].getMatchStart();

		uint64_t fA, fC, fG, fT;
		AF.getAdapterMatchFreqs(len,AOSPB.begin()[0],fA,fC,fG,fT);

		// confidence that this is not a random event
		double const randconfid = kmerPoisson(reflen,fA,fC,fG,fT,0/*n*/,pA,pC,pG,pT);
		
		if ( verbose > 1 )
		{
			std::cerr << "\n" << std::string(80,'-') << "\n\n";
		
			std::cerr << "read length " << len << " clip " << clip << " randconfig=" << randconfid << " fA=" << fA << " fC=" << fC << " fG=" << fG << " fT=" << fT << std::endl;
				
			for ( libmaus::bambam::AdapterOffsetStrand const * it = AOSPB.begin(); it != AOSPB.end(); ++it )
				AF.printAdapterMatch(ua,len,*it);
		}
		
		algn.putAuxNumber("as",'i',clip);
		algn.putAuxString("aa",AF.getAdapterName(AOSPB.begin()[0]));
		algn.putAuxNumber("af",'f',AOSPB.begin()[0].pfrac);
		algn.putAuxNumber("ar",'f',randconfid);
	}
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

	uint64_t const adpmatchminscore  = arginfo.getValue<uint64_t>("adpmatchminscore",getDefaultMatchMinScore());
	double   const adpmatchminfrac   = arginfo.getValue<double>("adpmatchminfrac",getDefaultMatchMinFrac());
	double   const adpmatchminpfrac  = arginfo.getValue<double>("adpmatchminpfrac",getDefaultMatchMinPFrac());
	uint64_t const reflen = arginfo.getValue<uint64_t>("reflen",getDefaultRefLen());
	double   const pA = arginfo.getValue<double>("pA",getDefaultpA());
	double   const pC = arginfo.getValue<double>("pC",getDefaultpC());
	double   const pG = arginfo.getValue<double>("pG",getDefaultpG());
	double   const pT = arginfo.getValue<double>("pT",getDefaultpT());

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const clip = arginfo.getValue<int>("clip",getDefaultClip());
	uint64_t const mod = arginfo.getValue<int>("mod",getDefaultMod());
	// length of seed
	uint64_t const seedlength = 
		std::min(
			static_cast<unsigned int>((8*sizeof(uint64_t))/3),
			std::max(arginfo.getValue<unsigned int>("SEED_LENGTH",getDefaultSEED_LENGTH()),1u));
	// maximum mismatch rate in overlap
	double const mismatchrate = std::min(100u,arginfo.getValue<unsigned int>("PCT_MISMATCH",getDefaultPCT_MISMATCH()))/100.0;
	// maximum number of mismatches in seed
	unsigned int const maxseedmismatches = arginfo.getValue<unsigned int>("MAX_SEED_MISMATCHES",static_cast<unsigned int>(std::ceil(seedlength * mismatchrate)));
	// minimum length of overlap between reads
	uint64_t const minoverlap = arginfo.getValue<uint64_t>("MIN_OVERLAP",getDefaultMIN_OVERLAP());
	// maximum number of adapter bases to be compared
	uint64_t const almax = arginfo.getValue<uint64_t>("ADAPTER_MATCH",getDefaultADAPTER_MATCH());

	libmaus::bambam::AdapterFilter::unique_ptr_type AF;
	
	if ( arginfo.hasArg("adaptersbam") )
	{
		libmaus::aio::CheckedInputStream adapterCIS(arginfo.getUnparsedValue("adaptersbam","adapters.bam"));
		libmaus::bambam::AdapterFilter::unique_ptr_type tAF(
                                new libmaus::bambam::AdapterFilter(adapterCIS,12 /* seed length */)
                        );
		AF = UNIQUE_PTR_MOVE(tAF);
	}
	else
	{
		std::string const builtinAdapters = libmaus::bambam::BamDefaultAdapters::getDefaultAdapters();
		std::istringstream builtinAdaptersStr(builtinAdapters);
		libmaus::bambam::AdapterFilter::unique_ptr_type tAF(
                                new libmaus::bambam::AdapterFilter(builtinAdaptersStr,12 /* seed length */)
                        );
		AF = UNIQUE_PTR_MOVE(tAF);
	}

	libmaus::autoarray::AutoArray<char> Aread;
	libmaus::util::PushBuffer<libmaus::bambam::AdapterOffsetStrand> AOSPB;

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

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string md5filename;
	std::string indexfilename;

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
			md5filename = arginfo.getUnparsedValue("md5filename","");
		else
			std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

		if ( md5filename.size() )
		{
			::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
			Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
			cbs.push_back(Pmd5cb.get());
		}
	}
	libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	if ( arginfo.getValue<unsigned int>("index",getDefaultIndex()) )
	{
		if ( arginfo.hasArg("indexfilename") &&  arginfo.getUnparsedValue("indexfilename","") != "" )
			indexfilename = arginfo.getUnparsedValue("indexfilename","");
		else
			std::cerr << "[V] no filename for index given, not creating index" << std::endl;

		if ( indexfilename.size() )
		{
			libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Tindex(new libmaus::bambam::BgzfDeflateOutputCallbackBamIndex(tmpfileindex));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
			cbs.push_back(Pindex.get());
		}
	}
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(std::cout,uphead,level,Pcbs));
	
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

	libmaus::bambam::BamAuxFilterVector auxfilter;
	auxfilter.set("a3");
	auxfilter.set("ah");
	
	uint64_t alcnt = 0;
	uint64_t adptcnt = 0;
	uint64_t lastalcnt = std::numeric_limits<uint64_t>::max();
	uint64_t paircnt = 0;
	uint64_t orphcnt = 0;

	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	// uint64_t const bmask = bmod-1;
	uint64_t const bshift = libmaus::math::ilog(bmod);
	
	libmaus::util::Histogram overlaphist;
	libmaus::util::Histogram adapterhist;

	libmaus::autoarray::AutoArray<char> CR;
	libmaus::autoarray::AutoArray<char> CQ;
	libmaus::bambam::BamSeqEncodeTable const seqenc;
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	libmaus::bambam::BamAlignment::D_array_type T;
	
	// std::cerr << "bmask=" << bmask << std::endl;
	
	while ( (running = bamdec.readAlignment()) )
	{
		if ( verbose && ( (alcnt >> bshift) != (lastalcnt >> bshift) ) )
		{
			std::cerr << "[V]\t" << alcnt << "\t" << paircnt << "\t" << adptcnt << std::endl;
			lastalcnt = alcnt;
		}
	
		alcnt++;

		// find adapters in given list
		adapterListMatch(Aread,AOSPB,inputalgn,*AF,verbose,adpmatchminscore,adpmatchminfrac,adpmatchminpfrac,reflen,pA,pC,pG,pT);

		// if this is a single end read, then write it back and try the next one
		if ( ! inputalgn.isPaired() )
		{
			if ( clip )
				clipAdapters(inputalgn,CR,CQ,seqenc,cigop,T);
			inputalgn.serialise(writer->getStream());
			continue;
		}
	
		// store alignment in algns[0]
		algns[0].swap(inputalgn);
		
		// read next alignment
		bool const okb = bamdec.readAlignment();
		alcnt++;
		
		// no next alignment, algns[0] is an orphan
		if ( ! okb )
		{
			if ( clip )
				clipAdapters(algns[0],CR,CQ,seqenc,cigop,T);
		
			++orphcnt;
			// std::cerr << "[D] warning: orphan alignment"  << std::endl;
			algns[0].serialise(writer->getStream());
			break;
		}
		
		// if next is not paired or name does not match
		if ( (! inputalgn.isPaired()) || strcmp(algns[0].getName(), inputalgn.getName()) )
		{
			if ( clip )
				clipAdapters(algns[0],CR,CQ,seqenc,cigop,T);
			++orphcnt;
			//std::cerr << "[D] warning: orphan alignment" << std::endl;
			algns[0].serialise(writer->getStream());
			bamdec.putback();
			alcnt--;
			continue;
		}
		
		// sanity checks
		assert ( algns[0].isPaired() );
		assert ( inputalgn.isPaired() );
		assert ( strcmp(algns[0].getName(),inputalgn.getName()) == 0 );

		// put second read in algns[1]
		algns[1].swap(inputalgn);

		// find adapters in given list
		adapterListMatch(Aread,AOSPB,algns[1],*AF,verbose,adpmatchminscore,adpmatchminfrac,adpmatchminpfrac,reflen,pA,pC,pG,pT);
		
		// are the read in the correct order? if not, write them out without touching them
		if ( !(algns[0].isRead1() && algns[1].isRead2()) )
		{
			std::cerr << "[D] warning: reads are not in the correct order" << std::endl;
			
			if ( clip )
			{
				clipAdapters(algns[0],CR,CQ,seqenc,cigop,T);
				clipAdapters(algns[1],CR,CQ,seqenc,cigop,T);
			}
			
			algns[0].serialise(writer->getStream());
			algns[1].serialise(writer->getStream());
			continue;
		}

		// are the reads both non empty?
		if ( !(algns[0].getLseq() && algns[1].getLseq()) )
		{
			std::cerr << "[D] warning: empty read" << std::endl;

			if ( clip )
			{
				clipAdapters(algns[0],CR,CQ,seqenc,cigop,T);
				clipAdapters(algns[1],CR,CQ,seqenc,cigop,T);
			}

			algns[0].serialise(writer->getStream());
			algns[1].serialise(writer->getStream());
			continue;
		}
		
		paircnt++;
		
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
			if ( clip )
			{
				clipAdapters(algns[0],CR,CQ,seqenc,cigop,T);
				clipAdapters(algns[1],CR,CQ,seqenc,cigop,T);
			}
			
			algns[0].serialise(writer->getStream());
			algns[1].serialise(writer->getStream());		
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
				// end of total match on read 1
				uint64_t const end0 = matchpos + lseedlength;
				// end of total match on read 2
				uint64_t const end1 = l1;
				// position of seed on read 1
				uint64_t const seedp0 = end0 - lseedlength;
				// position of seed on read 2
				uint64_t const seedp1 = end1 - lseedlength;
				// length of match between read 1 and read 2 (no indels allowed)
				uint64_t const restoverlap = std::min(seedp0,seedp1);

				// if length of overlap is sufficient
				if ( lseedlength + restoverlap >= minoverlap )
				{
					// iterator range of non seed overlap on read 1
					uint8_t const * check0  = reinterpret_cast<uint8_t const *>(seqs[0].begin()+seedp0-restoverlap);
					uint8_t const * check0e = check0 + restoverlap;
					// start iterator of non seed overlap on read 2
					uint8_t const * check1  = reinterpret_cast<uint8_t const *>(seqs[1].begin()+seedp1-restoverlap);
					
					// maximum number of mismatches allowed
					uint64_t const maxmis = (restoverlap+lseedlength) * mismatchrate;
					
					// compute number of mismatches
					uint64_t nummis = difcnt;
					while ( nummis <= maxmis && check0 != check0e  )
						if ( S[*check0++] != R[*check1++] )
							nummis++;
					
					// if match is within the maximal number of mismatches
					if ( check0 == check0e && nummis <= maxmis )
					{
						// update overlap histogram
						overlaphist(lseedlength + restoverlap);
					
						// adapter length for read 1
						uint64_t const al0 = (l0-end0);
						// adapter length for read 2
						uint64_t const al1 = end1-(restoverlap+lseedlength);
						// number of bases compared in this case
						uint64_t const alcmp = std::min(almax,std::min(al0,al1));

						// adapter range on read 1
						char const * ap0 = seqs[0].begin()+end0;
						char const * ap0e = ap0 + alcmp;
						// adapter start on read 2 (actually adapter end, we will read it backwards)
						char const * ap1 = seqs[1].begin()+al1;
						// number of mismatches on adapter (up to a length of almax)
						unsigned int aldif = 0;
						
						// calculate mismatches in adapter
						while ( ap0 != ap0e )
							if ( S[*(ap0++)] != R[libmaus::fastx::invertUnmapped(*(--ap1))] )
								aldif++;

						// if number of mismatches is within the tolerated range (0 at this time)
						if ( (ap0 == ap0e) && ! aldif )
						{
							if ( verbose > 1 )
							{
								std::cerr
									<< "[V2] overlap: " << algns[0].getName() 
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
									"[V2] " << std::string(
										seqs[0].begin()+seedp0-restoverlap,
										seqs[0].begin()+end0) << "\n" <<
									"[V2] " << std::string(
										seqs[1].begin()+seedp1-restoverlap,
										seqs[1].begin()+end1) << "\n";
									
								std::cerr << "[V2] assumed adapter on read 1 ["<<al0<<"]: "
									<< std::string(
										seqs[0].begin()+end0,
										seqs[0].begin()+l0
									) << std::endl;
								
								std::cerr << "[V2] assumed adapter on read 2 ["<<al1<<"]: "
									<< libmaus::fastx::reverseComplementUnmapped(std::string(
										seqs[1].begin(),
										seqs[1].begin()+al1)) << std::endl;
							}

							algns[0].filterOutAux(auxfilter);
							algns[1].filterOutAux(auxfilter);
							
							algns[0].putAuxNumber("ah",'i',1);
							algns[1].putAuxNumber("ah",'i',1);
							algns[0].putAuxNumber("a3",'i',al0);
							algns[1].putAuxNumber("a3",'i',al1);
							
							adptcnt += 1;
							adapterhist(lseedlength + restoverlap);
							
							break;
						}
					}
				}
			}
			
			--matchpos;
		} while ( q != qe );

		// std::cerr << "pair for " << algns[0].getName() << std::endl;
		
		if ( clip )
		{
			clipAdapters(algns[0],CR,CQ,seqenc,cigop,T);
			clipAdapters(algns[1],CR,CQ,seqenc,cigop,T);
		}
		
		algns[0].serialise(writer->getStream());		
		algns[1].serialise(writer->getStream());		
	}

	if ( verbose )
		std::cerr << "[V] processed " << alcnt << std::endl;

	if ( verbose )
		std::cerr << "[V] paircnt=" << paircnt << " adptcnt=" << adptcnt << " alcnt=" << alcnt << " orphcnt=" << orphcnt << std::endl;
		
	std::map<uint64_t,uint64_t> const overlaphistmap = overlaphist.get();
	std::map<uint64_t,uint64_t> const adapterhistmap = adapterhist.get();
	
	uint64_t totOverlaps = 0;
	uint64_t totAdapters = 0;
	
	// print histogram
	for ( std::map<uint64_t,uint64_t>::const_iterator ita = overlaphistmap.begin();
		ita != overlaphistmap.end(); ++ita )
	{
		// number of adapters for overlap
		uint64_t const ladpt = (adapterhistmap.find(ita->first) != adapterhistmap.end()) ?
			adapterhistmap.find(ita->first)->second : 0;
			
		std::cerr 
			<< std::setfill('0') << std::setw(3) << ita->first << std::setw(0)
			<< " "
			<< std::setfill('0') << std::setw(7) << ita->second << std::setw(0)
			<< " "
			<< std::setfill('0') << std::setw(7) << ladpt << std::setw(0)
			<< std::endl;
			
		totOverlaps += ita->second;
		totAdapters += ladpt;
	}

	double const pctOverlaps = (static_cast<double>(totOverlaps)/paircnt)*100.0;
	double const pctAdapters = (static_cast<double>(totAdapters)/paircnt)*100.0;
	
	std::cerr << "Processed " << std::setw(7) << std::setfill('0') << paircnt << std::setw(0) << " read pairs." << std::endl;
	std::cerr << "Found     " << std::setw(7) << std::setfill('0') << totOverlaps << std::setw(0) << " overlaps (" << pctOverlaps << ")" << std::endl;
	std::cerr << "Found     " << std::setw(7) << std::setfill('0') << totAdapters << std::setw(0) << " adapters (" << pctAdapters << ")" << std::endl;

	writer.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam::Licensing::formatNumber(getDefaultMod())+"]>", "print progress every mod'th line (if verbose>0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adaptersbam=<[]>", "list of adapters/primers stored in a BAM file (use internal list if not given)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "clip=<["+::biobambam::Licensing::formatNumber(getDefaultClip())+"]>", "clip off adapters (see bamadapterclip program)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "SEED_LENGTH=<["+::biobambam::Licensing::formatNumber(getDefaultSEED_LENGTH())+"]>", "length of seed for matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "PCT_MISMATCH=<["+::biobambam::Licensing::formatNumber(getDefaultPCT_MISMATCH())+"]>", "maximum percentage of mismatches in matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "MAX_SEED_MISMATCHES=<[SEED_LENGTH*PCT_MISMATCH]>", "maximum number of mismatches in seed (up to 2)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "MIN_OVERLAP=<["+::biobambam::Licensing::formatNumber(getDefaultMIN_OVERLAP())+"]>", "minimum overlap between mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "ADAPTER_MATCH=<["+::biobambam::Licensing::formatNumber(getDefaultADAPTER_MATCH())+"]>", "maximum adapter match check" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adpmatchminscore=<["+::biobambam::Licensing::formatNumber(getDefaultMatchMinScore())+"]>", "minimum score for adapter list matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adpmatchminfrac=<["+::biobambam::Licensing::formatFloatingPoint(getDefaultMatchMinFrac())+"]>", "minimum fraction for adapter list matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "adpmatchminpfrac=<["+::biobambam::Licensing::formatFloatingPoint(getDefaultMatchMinPFrac())+"]>", "minimum fraction of overlap for adapter list matching" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reflen=<["+::biobambam::Licensing::formatNumber(getDefaultRefLen())+"]>", "length of reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pA=<["+::biobambam::Licensing::formatFloatingPoint(getDefaultpA())+"]>", "relative frequency of base A in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pC=<["+::biobambam::Licensing::formatFloatingPoint(getDefaultpC())+"]>", "relative frequency of base C in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pG=<["+::biobambam::Licensing::formatFloatingPoint(getDefaultpG())+"]>", "relative frequency of base G in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pT=<["+::biobambam::Licensing::formatFloatingPoint(getDefaultpT())+"]>", "relative frequency of base T in reference sequence/genome" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

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

