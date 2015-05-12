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
#include <libmaus2/aio/PosixFdOutputStream.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/math/iabs.hpp>
#include <libmaus2/util/GrowingFreeList.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/bambam/BamHeaderUpdate.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/StrCmpNum.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

#include <config.h>

enum link_type_enum { link_type_chain, link_type_cluster };

static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

static int getDefaultVerbose()
{
	return 1;
}

static int getDefaultDisableValidation()
{
	return false;
}

static std::string getDefaultInputFormat()
{
	return "bam";
}

struct AlignmentComparator
{
	bool operator()(::libmaus2::bambam::BamAlignment const * A, ::libmaus2::bambam::BamAlignment const * B)
	{
		// put mapped reads before 
		if ( A->isMapped() != B->isMapped() )
			return A->isMapped();
		// dont care about order of unmapped reads
		if ( ! A->isMapped() )
			return false;
		if ( A->getRefID() != B->getRefID() )
			return A->getRefID() < B->getRefID();
		else
			return A->getPos() < B->getPos();
	}
};

struct MappingRegion
{
	uint64_t readFrom;
	uint64_t readTo;
	uint64_t refId;
	uint64_t refFrom;
	uint64_t refTo;
	libmaus2::bambam::BamAlignment * algn;
	
	MappingRegion() : readFrom(0), readTo(0), refFrom(0), refTo(0), algn(0) {}
	MappingRegion(
		uint64_t const rreadFrom,
		uint64_t const rreadTo,
		uint64_t const rrefId,
		uint64_t const rrefFrom,
		uint64_t const rrefTo,
		libmaus2::bambam::BamAlignment * ralgn
	) : readFrom(rreadFrom), readTo(rreadTo), refId(rrefId), refFrom(rrefFrom), refTo(rrefTo), algn(ralgn) {}
};

std::ostream & operator<<(std::ostream & out, MappingRegion const & MR)
{
	return out << "MappingRegion(" 
		<< MR.readFrom << ","
		<< MR.readTo << ","
		<< MR.refId << ","
		<< MR.refFrom << ","
		<< MR.refTo << ")";
}

struct MappingRegionReadCoordComparator
{
	bool operator()(MappingRegion const & A, MappingRegion const & B) const
	{
		if ( A.readFrom != B.readFrom )
			return A.readFrom < B.readFrom;
		else
			return A.readTo < B.readTo;
	}
};

struct MappingRegionRefCoordComparator
{
	bool operator()(MappingRegion const & A, MappingRegion const & B) const
	{
		if ( A.refId != B.refId )
			return A.refId < B.refId;
		else if ( A.refFrom != B.refFrom )
			return A.refFrom < B.refFrom;
		else 
			return A.refTo < B.refTo;
	}
};

struct ScoredInterval
{
	uint64_t low;
	uint64_t high;
	double score;
	
	ScoredInterval() : low(0), high(0), score(0) {}
	ScoredInterval(
		uint64_t const rlow,
		uint64_t const rhigh,
		double const rscore
	) : low(rlow), high(rhigh), score(rscore) {}
};

struct ScoredIntervalComparator
{
	bool operator()(ScoredInterval const & A, ScoredInterval const & B) const
	{
		return A.score > B.score;
	}
};

struct ScoredAlignment;
std::ostream & operator<<(std::ostream & out, ScoredAlignment const & S);

struct ScoredAlignment
{
	friend std::ostream & operator<<(std::ostream & out, ScoredAlignment const & S);

	private:	
	libmaus2::bambam::BamAlignment * algn;
	int64_t score;
	uint64_t readfrom;
	uint64_t readto;
	uint64_t reffrom;
	uint64_t refto;
	
	public:
	ScoredAlignment() : algn(0), score(0), readfrom(0), readto(0), reffrom(0), refto(0) {}
	ScoredAlignment(
		libmaus2::bambam::BamAlignment * const ralgn,
		int64_t const rscore,
		uint64_t const rreadfrom,
		uint64_t const rreadto,
		uint64_t const rreffrom,
		uint64_t const rrefto
	)  : algn(ralgn), score(rscore), readfrom(rreadfrom), readto(rreadto), reffrom(rreffrom), refto(rrefto) {}
	
	double getNormalisedScore() const
	{
		return score / static_cast<double>(readto-readfrom);
	}
	
	bool operator<(ScoredAlignment const & O) const
	{
		// return score < O.score;
		return getNormalisedScore() < O.getNormalisedScore();
	}
	
	bool overlaps(ScoredAlignment const & O) const
	{
		typedef libmaus2::math::IntegerInterval<uint64_t> int_type;
		return ! int_type::intersection(int_type(readfrom,readto-1),int_type(O.readfrom,O.readto-1)).isEmpty();
	}
	
	bool overlaps(std::vector<ScoredAlignment> const & V) const
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( overlaps(V[i]) )
				return true;
		return false;
	}
	
	libmaus2::bambam::BamAlignment * getAlignment()
	{
		return algn;
	}
};

std::ostream & operator<<(std::ostream & out, ScoredAlignment const & S)
{
	return
		out << "ScoredAlignment(" << S.score << "," << S.readfrom << "," << S.readto 
			<< "," << S.reffrom
			<< "," << S.refto
			<< "," << S.getNormalisedScore() << ")";
}

void handleVector(
	std::vector< ::libmaus2::bambam::BamAlignment * > & samename,
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> & cigop,
	libmaus2::bambam::BamHeader const & header,
	uint64_t & readbases,
	uint64_t & mappedbases,
	libmaus2::bambam::BamBlockWriterBase & writer,
	link_type_enum const link_type,
	double const erate = 0.3
)
{
	std::sort(samename.begin(),samename.end());

	typedef std::pair<uint64_t,uint64_t> up;
	std::vector<MappingRegion> M;
	int64_t readlen = -1;
	
	typedef ScoredAlignment scored_alignment_type;
	std::priority_queue < scored_alignment_type > PQ;
	
	for ( uint64_t i = 0; i < samename.size() && samename[i]->isMapped(); )
	{
		uint64_t j = i;
		while ( 
			j < samename.size() && 
			samename[j]->isMapped() && 
			samename[j]->getRefID() == samename[i]->getRefID()
		)
			++j;
		
		for ( uint64_t z = i ; z < j; ++z )
		{
			::libmaus2::bambam::BamAlignment * palgn = samename[z];
			::libmaus2::bambam::BamAlignment & algn = *palgn;

			uint32_t const numcigop = algn.getCigarOperations(cigop);
			uint64_t const seqlen = algn.getLseq();
			uint64_t readpos = 0;
			uint64_t refpos = algn.getPos();
			uint64_t refid = algn.getRefID();

			uint64_t hleft = 0;
			uint64_t hright = 0;
			uint64_t sleft = 0;
			uint64_t sright = 0;
			uint64_t cl = 0;
			uint64_t cr = numcigop;
				
			while ( cl < cr && cigop[cl].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP )
				hleft += cigop[cl++].second;
			while ( cr > cl && cigop[cr-1].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP )
				hright += cigop[--cr].second;

			while ( cl < cr && cigop[cl].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP )
				sleft += cigop[cl++].second;
			while ( cr > cl && cigop[cr-1].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP )
				sright += cigop[--cr].second;
				
			readpos += sleft;

			for ( uint64_t ci = cl; ci < cr; ++ci )
				if ( 
					cigop[ci].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP
					||
					cigop[ci].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP
				)
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "Malformed cigar string in read " << algn.getName() << "\n";
					lme.getStream() << algn.formatAlignment(header) << "\n";
					lme.finish();
					throw lme;
				}
				
			uint64_t const mpre = readpos;
			uint64_t const rpre = refpos;
			
			for ( uint64_t ci = cl; ci < cr; ++ci )
			{
				uint64_t const ciglen = cigop[ci].second;
				
				switch ( cigop[ci].first )
				{
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH:
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
					{
						readpos += ciglen;
						refpos += ciglen;
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
					{
						readpos += ciglen;
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
					{
						// deleting bases from the reference
						refpos += ciglen;
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CREF_SKIP:
					{
						// skip bases on reference
						refpos += ciglen;
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP:
					{
						// skip bases on read
						readpos += ciglen;
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP:
					{
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CPAD:
					{
						break;
					}
				}
			}
			
			uint64_t const mpost = readpos;
			uint64_t const rpost = refpos;

			M.push_back(MappingRegion(hleft+mpre,hleft+mpost,refid,rpre,rpost,&algn));

			int64_t const score = palgn->getAuxAsNumber<int64_t>("AS");
			PQ.push(scored_alignment_type(samename[i],score,hleft+mpre,hleft+mpost,rpre,rpost));
			
			readpos += sright;

			if ( readpos != seqlen )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Malformed cigar string in read " << algn.getName() << "\n";
				lme.getStream() << algn.formatAlignment(header) << "\n";
				lme.getStream() << "cl=" << cl << " cr=" << cr << std::endl;
				lme.finish();
				throw lme;
			}

			assert ( readpos == seqlen );
			
			if ( readlen < 0 )
				readlen = seqlen + hleft + hright;
			else
				assert ( static_cast<uint64_t>(readlen) == seqlen + hleft + hright );
		}
		
		i = j;
	}

	std::vector<scored_alignment_type> scoredout;
	while ( ! PQ.empty() )
	{
		scored_alignment_type const SA = PQ.top(); PQ.pop();
		
		if ( ! SA.overlaps(scoredout) )
			scoredout.push_back(SA);
	}
	
	std::sort(M.begin(),M.end(),MappingRegionReadCoordComparator());

	std::vector < libmaus2::math::IntegerInterval<uint64_t> > PIIV;
	for ( uint64_t i = 0; i < M.size(); ++i )
		if ( M[i].readFrom < M[i].readTo )
			PIIV.push_back(libmaus2::math::IntegerInterval<uint64_t>(M[i].readFrom,M[i].readTo-1));
	
	std::vector < libmaus2::math::IntegerInterval<uint64_t> > const IIV = libmaus2::math::IntegerInterval<uint64_t>::mergeOverlapping(PIIV);

	readbases += readlen;

	// count bases used in mapped regions of read
	for ( uint64_t i = 0; i < IIV.size(); ++i )
		mappedbases += (IIV[i].to - IIV[i].from + 1);
	
	// sort by ref coordinates	
	std::sort(M.begin(),M.end(),MappingRegionRefCoordComparator());
	// maximum error rate for linking fragments
	// linked fragments
	std::vector<ScoredInterval> alintervals;
	
	for ( uint64_t refidlow = 0; refidlow != M.size(); )
	{
		uint64_t refidhigh = refidlow;
		while ( refidhigh != M.size() && M[refidhigh].refId == M[refidlow].refId )
			++refidhigh;
			
		for ( uint64_t readlow = refidlow; readlow != refidhigh; )
		{
			uint64_t readhigh = readlow+1;
					
			while ( 
				readhigh != refidhigh &&
				// increasing read coordinates
				M[readhigh].readFrom >= M[readhigh-1].readTo &&
				// increasing pos coordinates
				M[readhigh].refFrom >= M[readhigh-1].refTo &&
				// check error rate implied by distance of fragments on read and reference
				(
					(
						libmaus2::math::iabs(
							static_cast<int64_t>(M[readhigh].readFrom - M[readhigh-1].readTo)
							-
							static_cast<int64_t>(M[readhigh].refFrom - M[readhigh-1].refTo)
						)
					)
					<=
					(
						(1+erate) *
							std::min(
								M[readhigh].readFrom - M[readhigh-1].readTo,
								M[readhigh].refFrom - M[readhigh-1].refTo
							)
					)
				)
			)
			{
				++readhigh;
			}
			
			#if 0
			if ( readhigh - readlow > 1 )
			{
				std::cerr << std::string(80,'-') << std::endl;
				for ( uint64_t i = readlow; i != readhigh; ++i )
				{
					std::cerr << M[i];
				}
				std::cerr << std::endl;
			}
			#endif
			double score = 0;
			for ( uint64_t i = readlow; i < readhigh; ++i )
				score += static_cast<double>(M[i].algn->getAuxAsNumber<int64_t>("AS"));
			
			alintervals.push_back(ScoredInterval(readlow,readhigh,score));
			
			readlow = readhigh;
		}
			
		refidlow = refidhigh;
	}
	
	// sort by descending score
	std::sort(alintervals.begin(),alintervals.end(),ScoredIntervalComparator());
	
	std::vector< ::libmaus2::bambam::BamAlignment * > outputvec;
	
	if ( link_type == link_type_chain )
	{
		if ( alintervals.size() )
		{
			ScoredInterval const SI = alintervals.front();
			
			#if 0
			std::cerr << std::string(80,'-') << std::endl;
			for ( uint64_t i = SI.low; i != SI.high; ++i )
				std::cerr << M[i];
			std::cerr << std::endl;
			#endif
			
			for ( uint64_t i = SI.low; i < SI.high; ++i )
				outputvec.push_back(M[i].algn);			
		}
	}
	else if ( link_type == link_type_cluster )
	{	
		for ( uint64_t i = 0; i < scoredout.size(); ++i )
			outputvec.push_back(scoredout[i].getAlignment());
	}
	
	if ( outputvec.size() == 1 )
	{
		writer.writeAlignment(*(outputvec[0]));
	}
	else if ( outputvec.size() > 1 )
	{
		for ( uint64_t i = 1; i < outputvec.size(); ++i )
		{
			outputvec[i]->putFlags(outputvec[i]->getFlags() | libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSUPPLEMENTARY);
		}
		for ( uint64_t i = 0; i < outputvec.size(); ++i )
		{
			uint64_t const next = (i+1) % outputvec.size();
			outputvec[i]->putNextPos(outputvec[next]->getPos());
			outputvec[i]->putNextRefId(outputvec[next]->getRefID());
		}
		for ( uint64_t i = 0; i < outputvec.size(); ++i )
			writer.writeAlignment(*(outputvec[i]));
	}
}

int bamlastfilter(libmaus2::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue("verbose",getDefaultVerbose());
	std::string const reference = arginfo.getUnparsedValue("reference",std::string());
	std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpfile",arginfo.getDefaultTmpFileName());	
	link_type_enum link_type = link_type_chain;
	
	if ( arginfo.hasArg("linktype") )
	{
		std::string const ltype = arginfo.getUnparsedValue("linktype","");
		
		if ( ltype == "chain" )
			link_type = link_type_chain;
		else if ( ltype == "cluster" )
			link_type = link_type_cluster;
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Unknown linktype " << ltype << "\n";
			lme.getStream() << "Supported options are chain and cluster\n";
			lme.finish();
			throw lme;		
		}
	}
	
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();	
	::libmaus2::bambam::BamAlignment & algn = dec.getAlignment();
	::libmaus2::bambam::BamAlignment prevalgn;
	bool haveprevalgn = false;
	uint64_t alcnt = 0;
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string md5filename;
	std::string indexfilename;

	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
			md5filename = arginfo.getUnparsedValue("md5filename","");
		else
			std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

		if ( md5filename.size() )
		{
			::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus2::lz::BgzfDeflateOutputCallbackMD5);
			Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
			cbs.push_back(Pmd5cb.get());
		}
	}
	libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	if ( arginfo.getValue<unsigned int>("index",getDefaultIndex()) )
	{
		if ( arginfo.hasArg("indexfilename") &&  arginfo.getUnparsedValue("indexfilename","") != "" )
			indexfilename = arginfo.getUnparsedValue("indexfilename","");
		else
			std::cerr << "[V] no filename for index given, not creating index" << std::endl;

		if ( indexfilename.size() )
		{
			libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Tindex(new libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex(tmpfileindex));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
			cbs.push_back(Pindex.get());
		}
	}
	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	::libmaus2::bambam::BamHeader::unique_ptr_type genuphead(
		libmaus2::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bamlastfilter",std::string(PACKAGE_VERSION))
	);
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(libmaus2::bambam::BamBlockWriterBaseFactory::construct(*genuphead,arginfo,Pcbs));
	libmaus2::bambam::BamBlockWriterBase & wr = *Pwriter;

	libmaus2::util::GrowingFreeList< ::libmaus2::bambam::BamAlignment > alfl;
	std::vector< ::libmaus2::bambam::BamAlignment * > samename;
	uint64_t readbases = 0;
	uint64_t mappedbases = 0;
	
	while ( dec.readAlignment() )
	{
		if ( haveprevalgn )
		{
			int const r = libmaus2::bambam::StrCmpNum::strcmpnum(
				prevalgn.getName(),
				algn.getName()
			);
			
			if ( r > 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Input file is not sorted by query name\n";
				lme.finish();
				throw lme;
			}
		}
	
		// new name
		if ( samename.size() && strcmp(samename.back()->getName(),algn.getName()) )
		{
			handleVector(samename,cigop,header,readbases,mappedbases,wr,link_type);
		
			while ( samename.size() )
			{
				alfl.put(samename.back());
				samename.pop_back();
			}
		}

		::libmaus2::bambam::BamAlignment * calgn = alfl.get();
		calgn->copyFrom(algn);
		samename.push_back(calgn);
		
		algn.swap(prevalgn);
		haveprevalgn = true;
		
		if ( verbose && ((++alcnt % (1024*1024)) == 0) )
			std::cerr << "[V] " << alcnt << std::endl;
	}

	if ( samename.size() )
	{
		handleVector(samename,cigop,header,readbases,mappedbases,wr,link_type);

		while ( samename.size() )
		{
			alfl.put(samename.back());
			samename.pop_back();
		}
	}
	
	std::cerr << "[V]\treadbases=" << readbases << "\tmappedbases=" << mappedbases << std::endl;

	// reset BAM writer
	Pwriter.reset();

	if ( Pmd5cb )
		Pmd5cb->saveDigestAsFile(md5filename);
	if ( Pindex )
		Pindex->flush(std::string(indexfilename));


	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam2::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "linktype=<>", "method used for selecting mapped fragments (chain or cluster)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return bamlastfilter(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
