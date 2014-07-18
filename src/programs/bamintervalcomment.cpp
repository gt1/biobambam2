/**
    bambam
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
#include "config.h"

#include <iostream>
#include <queue>

#include <libmaus/aio/CheckedOutputStream.hpp>

#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/bambam/GeneFlatFile.hpp>
#include <libmaus/bambam/MdNmRecalculation.hpp>

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/lz/BufferedGzipStream.hpp>

#include <libmaus/regex/PosixRegex.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>


static int getDefaultMD5() { return 0; }
static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static bool getDefaultDisableValidation() { return false; }
static std::string getDefaultInputFormat() { return "bam"; }
static int getDefaultIndex() { return 0; }

struct NamedInterval
{
	uint64_t name;
	uint64_t refseq;
	uint64_t from;
	uint64_t to;
	uint64_t bin;
	
	NamedInterval() : name(0), refseq(0), from(0), to(0), bin(0) {}
	NamedInterval(
		uint64_t const rname,
		uint64_t const rrefseq,
		uint64_t const rfrom,
		uint64_t const rto,
		uint64_t const rbin = 0
	) : name(rname), refseq(rrefseq), from(rfrom), to(rto), bin(rbin) {}
	
	bool operator<(NamedInterval const & O) const
	{
		if ( refseq != O.refseq )
			return refseq < O.refseq;
		else if ( from != O.from )
			return from < O.from;
		else if ( to != O.to )
			return to < O.to;
		else if ( name != O.name )
			return name < O.name;
		else
			return false;
	}
};

struct NamedIntervalLowComparator
{
	bool operator()(NamedInterval const & A, NamedInterval const & B) const
	{
		if ( A.refseq != B.refseq )
			return A.refseq < B.refseq;
		else
			return A.from < B.from;
	}
};

struct NamedIntervalBinComparator
{
	bool operator()(NamedInterval const & A, NamedInterval const & B) const
	{
		if ( A.refseq != B.refseq )
			return A.refseq < B.refseq;
		else
			return A.bin < B.bin;
	}
};

std::ostream & operator<<(std::ostream & out, NamedInterval const & N)
{
	return out << "NamedInterval(" << N.name << "," << N.refseq << "," << N.from << "," << N.to << ")";
}

struct NamedIntervalGeneMeta
{
	std::string genename;
	std::string name;
	
	NamedIntervalGeneMeta()
	: genename(), name()
	{}
	NamedIntervalGeneMeta(std::string const & rgenename, std::string const & rname)
	: genename(rgenename), name(rname)
	{}
};

std::ostream & operator<<(std::ostream & out, NamedIntervalGeneMeta const & N)
{
	if ( N.genename == N.name )
		return out << "(" << N.genename << ")";
	else
		return out << "(" << N.genename << "," << N.name << ")";
}

struct NamedIntervalGeneSet
{
	static unsigned int const maxbinbits = 16;
	int ilog;
	uint64_t blockmask;

	std::vector<NamedInterval> intervals;
	std::vector<NamedIntervalGeneMeta> meta;
	std::vector< std::pair<uint64_t,uint64_t> > refidintervals;

	static uint64_t parseNumber(std::string const & s)
	{
		if ( ! s.size() )
		{
			libmaus::exception::LibMausException lme;
			lme.getStream() << "[E] cannot parse empty string as number." << std::endl;
			lme.finish();
			throw lme;	
		}

		unsigned char const * c = reinterpret_cast<unsigned char const *>(s.c_str());
		uint64_t n = 0;
		for ( uint64_t i = 0; i < s.size(); ++i )
			if ( ! isdigit(c[i]) )
			{
				libmaus::exception::LibMausException lme;
				lme.getStream() << "[E] cannot parse " << s << " as number." << std::endl;
				lme.finish();
				throw lme;
			}
			else
			{
				n *= 10;
				n += c[i]-'0';
			}
			
		return n;
	}
	
	static std::string unifyTranscripts(std::string const & fn, std::string const & refregex)
	{
		libmaus::regex::PosixRegex regex(refregex);
		libmaus::bambam::GeneFlatFile::unique_ptr_type PGFF(libmaus::bambam::GeneFlatFile::construct(fn));
		std::ostringstream ostr;
		libmaus::bambam::GeneFlatFileEntry entry;
		
		// compute map of chromosomes per gene
		std::map< std::string, std::set<std::string> > genechroms;

		for ( uint64_t i = 0; i < PGFF->size(); ++i )
		{
			PGFF->get(i,entry);
			std::string const chrom(entry.chrom.first,entry.chrom.second);

			if ( regex.findFirstMatch(chrom.c_str()) != -1 )
			{
				std::string const geneName(entry.geneName.first,entry.geneName.second);
				genechroms[geneName].insert(chrom);
			}
		}

		std::map< std::string,std::pair<uint64_t,uint64_t> > geneCoord;
		std::map< std::string,std::vector<uint64_t> > geneInst;
		
		for ( uint64_t i = 0; i < PGFF->size(); ++i )
		{
			PGFF->get(i,entry);
			std::string const chrom(entry.chrom.first,entry.chrom.second);

			if ( regex.findFirstMatch(chrom.c_str()) != -1 )
			{
				std::string geneName(entry.geneName.first,entry.geneName.second);

				// append ref seq name if gene appears on multiple ref seqs
				if ( genechroms.find(geneName)->second.size() > 1 )
				{
					uint64_t minsize = std::numeric_limits<uint64_t>::max();
					
					for ( std::set<std::string>::const_iterator ita = genechroms.find(geneName)->second.begin();
						ita != genechroms.find(geneName)->second.end(); ++ita )
						minsize = std::min(
							static_cast<uint64_t>(minsize),
							static_cast<uint64_t>(ita->size())
						);
						
					uint64_t minsizecnt = 0;
					for ( std::set<std::string>::const_iterator ita = genechroms.find(geneName)->second.begin();
						ita != genechroms.find(geneName)->second.end(); ++ita )
						if ( ita->size() == minsize )
							minsizecnt++;
					
					if ( 
						minsizecnt > 1 
						||
						chrom.size() > minsize
					)
						geneName += (std::string("_")+chrom);
				}
				
				std::map< std::string,std::pair<uint64_t,uint64_t> >::iterator ita =
					geneCoord.find(geneName);
			
				if ( ita == geneCoord.end() )
					geneCoord[geneName] = std::pair<uint64_t,uint64_t>(entry.txStart,entry.txEnd);
				else
				{
					ita->second.first  = std::min(ita->second.first ,entry.txStart);
					ita->second.second = std::max(ita->second.second,entry.txEnd);
				}

				geneInst[geneName].push_back(i);
			}
			else
			{
				// std::cerr << "[V] dropping " << chrom << std::endl; 
			}		
		}
		
		uint64_t nonuniquecnt = 0;
		for ( std::map< std::string,std::vector<uint64_t> >::const_iterator ita = geneInst.begin(); ita != geneInst.end(); ++ita )
		{
			std::string const & geneName = ita->first;
			std::vector<uint64_t> const & inst = ita->second;
			
			assert ( inst.size() > 0 );
			
			// single instance?
			if ( inst.size() == 1 )
			{
				PGFF->get(inst[0],entry);
				entry.geneName.first  = geneName.c_str();
				entry.geneName.second = geneName.c_str() + geneName.size();
				ostr << entry << "\n";
			}
			else
			{
				std::pair<uint64_t,uint64_t> maxcoord = geneCoord.find(geneName)->second;
				uint64_t maxcnt = 0, maxid = inst.size();
				
				for ( uint64_t i = 0; i < inst.size(); ++i )
				{
					PGFF->get(inst[i],entry);
					
					if ( entry.txStart == maxcoord.first && entry.txEnd == maxcoord.second )
					{
						maxcnt++;
						maxid = i;
					}
				}
				
				if ( maxcnt )
				{
					PGFF->get(inst[maxid],entry);
					
					entry.geneName.first  = geneName.c_str();
					entry.geneName.second = geneName.c_str() + geneName.size();

					if ( maxcnt > 1 )
					{
						entry.txStart = entry.cdsStart = maxcoord.first;
						entry.txEnd   = entry.cdsEnd   = maxcoord.second;
						entry.name = entry.geneName;
						entry.exons.resize(0);
					}
					ostr << entry << "\n";		
				}
				else
				{
					nonuniquecnt++;

					PGFF->get(inst[0],entry);

					entry.geneName.first  = geneName.c_str();
					entry.geneName.second = geneName.c_str() + geneName.size();

					entry.txStart = entry.cdsStart = maxcoord.first;
					entry.txEnd   = entry.cdsEnd   = maxcoord.second;
					entry.name    = entry.geneName;
					entry.exons.resize(0);

					ostr << entry << "\n";		
				}
			}
		}
		
		// std::cerr << "[V] found " << nonuniquecnt << " non unique transcript genes" << std::endl;
		
		// std::cerr << ostr.str();
		
		return ostr.str();
	}

	NamedIntervalGeneSet(
		std::string const & fn, libmaus::bambam::BamHeader const & header,
		bool const unify = false,
		std::string const chromregex = ".*"
	)
	: ilog(-1), intervals(), meta()
	{
		
		std::istream * Pistr = 0;
		libmaus::util::unique_ptr<std::istringstream>::type Puistr;		
		libmaus::aio::CheckedInputStream CIS(fn);
		libmaus::lz::BufferedGzipStream::unique_ptr_type PBGS;
		std::string unif;

		// std::cerr << "unify " << unify << std::endl;

		if ( unify )
		{
			unif = unifyTranscripts(fn,chromregex);
			libmaus::util::unique_ptr<std::istringstream>::type Tuistr(new std::istringstream(unif));
			Puistr = UNIQUE_PTR_MOVE(Tuistr);
			Pistr = Puistr.get();
		}
		else
		{
			bool const isgz = libmaus::bambam::MdNmRecalculation::isGzip(fn);
			
			if ( isgz )
			{
				libmaus::lz::BufferedGzipStream::unique_ptr_type TBGS(new libmaus::lz::BufferedGzipStream(CIS));
				PBGS = UNIQUE_PTR_MOVE(TBGS);
				Pistr = PBGS.get();
			}
			else
			{
				Pistr = &CIS;
			}
		}
		
		std::istream & istr = *Pistr;

		std::map<std::string,uint64_t> refmap;

		for ( uint64_t i = 0; i < header.getNumRef(); ++i )
		{
			refmap[header.getRefIDName(i)] = i;
		}
		
		refidintervals.resize(header.getNumRef());
		
		while ( istr )
		{
			std::string line;
			std::getline(istr,line);
			
			if ( line.size() )
			{
				std::deque<std::string> const Q = libmaus::util::stringFunctions::tokenize<std::string>(line,std::string("\t"));
				
				// ignore line if it does not have a  sufficient number of columns
				if ( Q.size() < 6 )
					continue;
					
				std::string const & genename = Q.at(0);
				std::string const & name = Q.at(1);
				std::string const & refname = Q.at(2);
				// std::string const & codingstrand = Q.at(3);
				std::string const & sstart = Q.at(4);
				std::string const & send = Q.at(5);
				
				std::map<std::string,uint64_t>::const_iterator it = refmap.find(refname);
				
				if ( it == refmap.end() )
				{
					std::cerr << "[W] cannot find reference " << refname << " in BAM file, not searching for gene " << genename << std::endl;
					continue;
					
					#if 0
					libmaus::exception::LibMausException lme;
					lme.getStream() << "[E] cannot find reference " << refname << " in BAM file" << std::endl;
					lme.finish();
					throw lme;
					#endif
				}

				uint64_t const refseq = it->second;
				uint64_t const start = parseNumber(sstart);
				uint64_t const end = parseNumber(send);
				
				NamedIntervalGeneMeta nmeta(genename,name);
				NamedInterval interval(meta.size(),refseq,start,end);
				
				intervals.push_back(interval);
				meta.push_back(nmeta);
			}
		}
		
		std::sort(intervals.begin(),intervals.end());
		
		uint64_t maxtop = 0;
		for ( uint64_t i = 0; i < intervals.size(); ++i )
			maxtop = std::max(maxtop,intervals[i].to);
			
		// next power of two
		uint64_t const ntp = libmaus::math::nextTwoPow(maxtop+1);
		// number of bits used
		ilog = libmaus::math::ilog(ntp);
		// mask
		int const maxbinbits1 = maxbinbits-1;
		blockmask = ~((1ull << ((ilog>=maxbinbits1) ? (ilog-maxbinbits1) : 0))-1);

		// compute bins for all intervals		
		for ( uint64_t i = 0; i < intervals.size(); ++i )
		{
			uint64_t bin = 0;
			uint64_t binadd = 1;
			uint64_t from = intervals[i].from;
			uint64_t to = intervals[i].to;

			for ( 
				uint64_t m = 1ull << (ilog-1); 
				(m!=0) && ((from & m) == (to & m)); 
				(m = ((m >> 1) & blockmask)), (binadd <<= 1) 
			)
			{
				bin <<= 1;
				bin += 1+((from&m)!=0);
			}
			
			#if 0
			std::cerr << "[" << from << "," << to << "] bin " << bin << std::endl;
			#endif
			
			intervals[i].bin = bin;
		}
		
		// sort by reference sequence, bin
		std::sort(intervals.begin(),intervals.end(),NamedIntervalBinComparator());
		
		// compute interval for each reference sequence
		uint64_t low = 0;
		while ( low != intervals.size() )
		{
			uint64_t high = low;
			while ( high != intervals.size() && intervals[high].refseq == intervals[low].refseq )
				++high;

			refidintervals[intervals[low].refseq] = std::pair<uint64_t,uint64_t>(low,high);
			
			low = high;
		}
	}

	struct FindIntervalsQueueElement
	{
		int s;
		uint64_t bin;
		uint64_t subfrom;
		uint64_t subto;
		
		FindIntervalsQueueElement() : s(0), bin(0), subfrom(0), subto(0) {}
		FindIntervalsQueueElement(int const rs, uint64_t const rbin, uint64_t const rsubfrom, uint64_t const rsubto) 
		: s(rs), bin(rbin), subfrom(rsubfrom), subto(rsubto)
		{}
	};
	
	struct Interval
	{
		int64_t from;
		int64_t to;
		
		Interval() : from(0), to(0) {}
		Interval(int64_t const rfrom, int64_t const rto) : from(rfrom), to(rto) {}
		
		bool empty() const
		{
			return to < from;
		}
		
		static Interval intersection(Interval const & A, Interval const & B)
		{
			if ( A.empty() || B.empty() )
				return Interval(0,-1);
			if ( A.from > B.from )
				return intersection(B,A);
			
			assert ( A.from <= B.from );
			
			if ( A.to < B.from )
				return Interval(0,-1);
				
			return Interval(B.from,std::min(B.to,A.to));
		}
		
		Interval intersection(Interval const & B) const
		{
			return intersection(*this,B);
		}
	};

	void findIntervals(
		libmaus::bambam::BamAlignment const & algn,
		std::vector<uint64_t> & binMatchingIntervals
	)
	const
	{
		binMatchingIntervals.resize(0);
	
		if ( algn.isMapped() )
		{
			uint64_t const refid = algn.getRefID();
			uint64_t const from = algn.getPos();
			uint64_t const to = algn.getAlignmentEnd();
			Interval const A(from,to);

			#if defined(FIND_INTERVALS_DEBUG)
			#if 0
			std::cerr << "from=" << from << " to=" << to << std::endl;
			#endif
			std::vector<uint64_t> matchingIntervals;
						
			for ( uint64_t i = refidintervals[refid].first; i < refidintervals[refid].second; ++i )
			{
				assert ( intervals[i].refseq == refid );

				if ( ! A.intersection(Interval(intervals[i].from,intervals[i].to)).empty() )
				{
					matchingIntervals.push_back(i);
					#if 0
					std::cerr << "matching " << intervals[i] << std::endl;						
					#endif
				}
			}
			#endif
			
			std::deque<FindIntervalsQueueElement> Q;
			Q.push_back(FindIntervalsQueueElement(ilog-1,0,from,to));

			while ( Q.size() )
			{
				FindIntervalsQueueElement C = Q.front();
				Q.pop_front();
				
				#if 0
				std::cerr << "[" << C.subfrom << "," << C.subto << "] bin " << C.bin << std::endl;
				#endif
				
				std::vector<NamedInterval>::const_iterator it = std::lower_bound(
					intervals.begin()+refidintervals[refid].first,
					intervals.begin()+refidintervals[refid].second,
					NamedInterval(0,refid,0,0,C.bin),
					NamedIntervalBinComparator()
				);
				
				for ( ; it != intervals.end() && it->bin == C.bin; ++it )
					if ( ! A.intersection(Interval(it->from,it->to)).empty() )
					{
						binMatchingIntervals.push_back(it-intervals.begin());
						
						#if 0
						std::cerr << "bin-matching " << (*it) << std::endl;									
						#endif
					}
			
				uint64_t const m = (1ull << C.s) & blockmask;
				
				if ( m )
				{
					if ( ( C.subfrom & m ) != (C.subto & m) )
					{
						uint64_t leftfrom = C.subfrom;
						uint64_t leftto = C.subfrom | libmaus::math::lowbits(C.s);
					
						uint64_t rightfrom = C.subto & (~libmaus::math::lowbits(C.s));
						uint64_t rightto = C.subto;

						#if 0					
						std::cerr << "split [" << C.subfrom << "," << C.subto << "] into [" << leftfrom << "," << leftto << "] and ["
							<< rightfrom << "," << rightto << "]\n";
						#endif
							
						Q.push_back(
							FindIntervalsQueueElement(
								C.s-1,
								(C.bin<<1)+1,
								leftfrom,
								leftto
							)
						);
						Q.push_back(
							FindIntervalsQueueElement(
								C.s-1,
								(C.bin<<1)+2,
								rightfrom,
								rightto
							)
						);
					}
					else if ( !(C.subfrom & m) )
					{
						Q.push_back(
							FindIntervalsQueueElement(
								C.s-1,
								(C.bin<<1)+1,
								C.subfrom,
								C.subto
							)
						);
						
						#if 0
						std::cerr << "unsplit left [" << C.subfrom << "," << C.subto << "]" << std::endl;
						#endif
					}
					else
					{
						Q.push_back(
							FindIntervalsQueueElement(
								C.s-1,
								(C.bin<<1)+2,
								C.subfrom,
								C.subto
							)
						);				

						#if 0
						std::cerr << "unsplit right [" << C.subfrom << "," << C.subto << "]" << std::endl;
						#endif
					}
				}
			}
		
			#if defined(FIND_INTERVALS_DEBUG)
			std::sort(matchingIntervals.begin(),matchingIntervals.end());
			std::sort(binMatchingIntervals.begin(),binMatchingIntervals.end());
			
			assert ( matchingIntervals == binMatchingIntervals );
			#endif
		}
	}
};

std::ostream & operator<<(std::ostream & out, NamedIntervalGeneSet const & NIGS)
{
	out << "NamedIntervalGeneSet(\n";
	
	for ( uint64_t i = 0; i < NIGS.intervals.size(); ++i )
		out << "\t" << NIGS.intervals[i] << "\n";
	for ( uint64_t i = 0; i < NIGS.meta.size(); ++i )
		out << "\t" << NIGS.meta[i] << "\n";
	
	out << ")\n";
	
	return out;
}

int bamintervalcomment(::libmaus::util::ArgInfo const & arginfo)
{
	::libmaus::util::TempFileRemovalContainer::setup();
	
	bool const inputisstdin = (!arginfo.hasArg("I")) || (arginfo.getUnparsedValue("I","-") == "-");

	if ( isatty(STDIN_FILENO) && inputisstdin && (arginfo.getValue<std::string>("inputformat","bam") != "sam") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	if ( ! arginfo.hasArg("intervals") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "[E] required intervals key is missing" << std::endl;
		se.finish();
		throw se;
	}

	std::string const inputformat = arginfo.getUnparsedValue("inputformat",getDefaultInputFormat());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const disablevalidation = arginfo.getValue<int>("disablevalidation",getDefaultDisableValidation());

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	if ( disablevalidation )
		dec.disableValidation();
	::libmaus::bambam::BamHeader const & header = dec.getHeader();
	
	std::string const chromregex = arginfo.getUnparsedValue("chromregex",".*");
	bool const unifytranscripts = arginfo.getValue<unsigned int>("unifytranscripts",false);
	
	NamedIntervalGeneSet const NIGS(arginfo.getUnparsedValue("intervals","not set"),header,unifytranscripts,chromregex);
	
	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilenameout = tmpfilenamebase + "_bamintervalcomment";
	std::string const tmpfileindex = tmpfilenamebase + "_bamintervalcomment_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenameout);
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);
	
	::libmaus::bambam::BamHeader::unique_ptr_type genuphead(
		libmaus::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bamintervalcomment",std::string(PACKAGE_VERSION))
	);

	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;

	if ( arginfo.hasArg("md5") && arginfo.hasArg("md5filename") && arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Pmd5 = UNIQUE_PTR_MOVE(Tmd5);
		cbs.push_back(Pmd5.get());
	}

	libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	std::string indexfilename;
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
	
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(
		libmaus::bambam::BamBlockWriterBaseFactory::construct(
			*genuphead,arginfo,
			cbs.size() ? (&cbs) : 0
		)
	);

	libmaus::bambam::BamAlignment & curalgn = dec.getAlignment();
	libmaus::bambam::BamAuxFilterVector COfilter;
	COfilter.set("CO");
	uint64_t c = 0;
	std::vector<uint64_t> matchingIntervals;
	bool const coord = arginfo.getValue<unsigned int>("coord",0);
	
	while ( dec.readAlignment() )
	{
		NIGS.findIntervals(curalgn,matchingIntervals);
		
		if ( matchingIntervals.size() )
		{
			std::ostringstream ostr;
			for ( uint64_t i = 0; i < matchingIntervals.size(); ++i )
			{
				ostr << ((i>0)?";":"");
				
				NamedInterval const & NI = NIGS.intervals[matchingIntervals[i]];
				uint64_t const nameid = NI.name;
				NamedIntervalGeneMeta const & meta = NIGS.meta[nameid];
				
				if ( coord )
				{
					ostr << "(";

					#if 0					
					if ( meta.name != meta.genename )
						ostr << meta.genename << "," << meta.name;
					else
					#endif
						ostr << meta.genename;
					
					ostr << "," << header.getRefIDName(NI.refseq);
					ostr << "," << NI.from;
					ostr << "," << NI.to;
					
					ostr << ")";					
				}
				else
					ostr << NIGS.meta[nameid];
			}
			curalgn.filterOutAux(COfilter);
			curalgn.putAuxString("CO",ostr.str());
		}
		
		Pwriter->writeAlignment(curalgn);

		if ( (++c & (1024*1024-1)) == 0 && verbose )
			std::cerr << "[V] " << c << std::endl;		
	}
	
	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

	Pwriter.reset();

	if ( Pmd5 )
		Pmd5->saveDigestAsFile(arginfo.getUnparsedValue("md5filename","not set"));
	if ( Pindex )
		Pindex->flush(std::string(indexfilename));
	
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
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );				

				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "intervals=<>", "intervals file (plain tab separated or tab separated gzip compressed)" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamintervalcomment(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
