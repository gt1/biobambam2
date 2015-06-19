/*
    biobambam2
    Copyright (C) 2015 German Tischler

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
*/
#include <iostream>
#include <cassert>

#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/dazzler/align/AlignmentFile.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/lcs/EditDistance.hpp>
#include <libmaus2/lcs/ND.hpp>
#include <libmaus2/lcs/NDextend.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/util/ArgInfo.hpp>

#include <libmaus2/fastx/FastAIndex.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/parallel/SynchronousQueue.hpp>
#include <libmaus2/parallel/TerminatableSynchronousQueue.hpp>

#include <config.h>

int getDefaultLevel()
{
	return -1;
}

int getDefaultVerbose()
{
	return 0;
}

int getDefaultMD5()
{
	return 0;
}

int getDefaultCalMdNm()
{
	return 0;
}

struct RgInfo
{
	std::string ID;
	std::string CN;
	std::string DS;
	std::string DT;
	std::string FO;
	std::string KS;
	std::string LB;
	std::string PG; // = fastqtobam
	std::string PI;
	std::string PL;
	std::string PM;
	std::string PU;
	std::string SM;
	
	RgInfo() {}
	RgInfo(libmaus2::util::ArgInfo const & arginfo)
	:
		ID(arginfo.getUnparsedValue("RGID","")),
		CN(arginfo.getUnparsedValue("RGCN","")),
		DS(arginfo.getUnparsedValue("RGDS","")),
		DT(arginfo.getUnparsedValue("RGDT","")),
		FO(arginfo.getUnparsedValue("RGFO","")),
		KS(arginfo.getUnparsedValue("RGKS","")),
		LB(arginfo.getUnparsedValue("RGLB","")),
		PG(arginfo.getUnparsedValue("RGPG","lasToBAM")),
		PI(arginfo.getUnparsedValue("RGPI","")),
		PL(arginfo.getUnparsedValue("RGPL","")),
		PM(arginfo.getUnparsedValue("RGPM","")),
		PU(arginfo.getUnparsedValue("RGPU","")),
		SM(arginfo.getUnparsedValue("RGSM",""))
	{
		
	}
	
	std::string toString() const
	{
		std::ostringstream ostr;
		
		if ( ID.size() )
		{
			ostr << "@RG\tID:" << ID;
			
			if ( CN.size() ) ostr << "\tCN:" << CN;
			if ( DS.size() ) ostr << "\tDS:" << DS;
			if ( DT.size() ) ostr << "\tDT:" << DT;
			if ( FO.size() ) ostr << "\tFO:" << FO;
			if ( KS.size() ) ostr << "\tKS:" << KS;
			if ( LB.size() ) ostr << "\tLB:" << LB;
			if ( PG.size() ) ostr << "\tPG:" << PG;
			if ( PI.size() ) ostr << "\tPI:" << PI;
			if ( PL.size() ) ostr << "\tPL:" << PL;
			if ( PM.size() ) ostr << "\tPM:" << PM;
			if ( PU.size() ) ostr << "\tPU:" << PU;
			if ( SM.size() ) ostr << "\tSM:" << SM;
			
			ostr << "\n";
		}
		
		return ostr.str();
	}
};

std::string getUsage(libmaus2::util::ArgInfo const & arginfo)
{
	std::ostringstream ostr;
	ostr << "usage: " << arginfo.progname << " [key=value pairs] <reads.db> <alignments.las> or <reads.db> <reads.db> <alignments.las>";
	return ostr.str();
}

struct LASToBAMConverter
{
	typedef LASToBAMConverter this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::dazzler::db::DatabaseFile const & DB2;

	std::vector<libmaus2::dazzler::db::Read> const & readsMeta1;
	std::vector<libmaus2::dazzler::db::Read> const & readsMeta2;

	libmaus2::lcs::EditDistanceTraceContainer ATC;		
	libmaus2::lcs::ND ND;
	::libmaus2::fastx::UCharBuffer ubuffer;
	::libmaus2::fastx::UCharBuffer wbuffer;
	::libmaus2::fastx::UCharBuffer sbuffer;
	::libmaus2::bambam::MdStringComputationContext context;
	::libmaus2::bambam::BamSeqEncodeTable const seqenc;

	char const * aptr;
	
	libmaus2::aio::InputStream::unique_ptr_type PbaseStreamB;
	
	bool const loadallb;
	libmaus2::autoarray::AutoArray<char> const * const readsB;
	std::vector<uint64_t> const * const readsOffB;

	libmaus2::autoarray::AutoArray<char> Bspace;
	libmaus2::autoarray::AutoArray<char> Binvspace;
	int64_t bid;
	char const * bbaseptr;
	char const * bptr;
	bool Binvspacevalid;
	
	uint64_t const tspace;

	bool const calmdnm;

	libmaus2::autoarray::AutoArray<std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Acigop;
	libmaus2::autoarray::AutoArray<uint8_t> ASQ;
	
	enum supplementary_seq_strategy_t {
		supplementary_seq_strategy_soft,
		supplementary_seq_strategy_hard,
		supplementary_seq_strategy_none
	};
	
	supplementary_seq_strategy_t const supplementary_seq_strategy;

	LASToBAMConverter(
		libmaus2::dazzler::db::DatabaseFile const & rDB2,
		std::vector<libmaus2::dazzler::db::Read> const & rreadsMeta1,
		std::vector<libmaus2::dazzler::db::Read> const & rreadsMeta2,
		libmaus2::aio::InputStream::unique_ptr_type & rPbaseStreamB,
		bool const rloadallb,
		libmaus2::autoarray::AutoArray<char> const * rreadsB,
		std::vector<uint64_t> const * rreadsOffB,
		uint64_t const rtspace, /* algn.tspace */
		bool const rcalmdnm,
		supplementary_seq_strategy_t const rsupplementaryStrategy
	)
	: 
	  DB2(rDB2),
	  readsMeta1(rreadsMeta1), readsMeta2(rreadsMeta2),
	  aptr(NULL),
	  PbaseStreamB(UNIQUE_PTR_MOVE(rPbaseStreamB)),
	  loadallb(rloadallb),
	  readsB(rreadsB),
	  readsOffB(rreadsOffB),
	  Bspace(),
	  Binvspace(),
	  bid(-1),
	  bbaseptr(NULL),
	  bptr(NULL),
	  Binvspacevalid(false),
	  tspace(rtspace),
	  calmdnm(rcalmdnm),
	  supplementary_seq_strategy(rsupplementaryStrategy)
	{
	
	}
	
	void setaptr(char const * aptr)
	{
		this->aptr = aptr;
	}
	
	std::pair<uint64_t,uint64_t> operator()(
		libmaus2::dazzler::align::Overlap const & OVL,
		::libmaus2::fastx::UCharBuffer & tbuffer,
		bool const primary
	)
	{
		#if 0
		// check path
		int32_t p = OVL.path.bbpos;
		for ( size_t i = 0; i < OVL.path.path.size(); ++i )
			p += OVL.path.path[i].second;
		assert ( p == OVL.path.bepos );
		#endif
		
		if ( OVL.bread != bid )
		{
			libmaus2::dazzler::db::Read const & R = readsMeta2.at(OVL.bread);
			
			if ( loadallb )
			{
				#if 0
				if ( R.rlen > static_cast<int64_t>(Bspace.size()) )
					Bspace.resize(R.rlen);
				std::copy(readsA.begin()+readsOffA[OVL.bread],readsA.begin()+readsOffA[OVL.bread+1],Bspace.begin());
				#endif
				bbaseptr = readsB->begin()+(*readsOffB)[OVL.bread];
			}
			else
			{
				PbaseStreamB->clear();
				PbaseStreamB->seekg(R.boff);
				libmaus2::dazzler::db::DatabaseFile::decodeRead(*PbaseStreamB,Bspace,R.rlen);
				bbaseptr = Bspace.begin();
			}
			
			bid = OVL.bread;
			Binvspacevalid = false;						
		}
		
		if ( OVL.isInverse() )
		{
			if ( ! Binvspacevalid )
			{
				libmaus2::dazzler::db::Read const & R = readsMeta2.at(OVL.bread);
			
				if ( R.rlen > static_cast<int64_t>(Binvspace.size()) )
					Binvspace.resize(R.rlen);
			
				char const * pin =  bbaseptr;
				char * pout = Binvspace.begin() + R.rlen;
			
				for ( int64_t i = 0; i < R.rlen; ++i )
					*(--pout) = libmaus2::fastx::invertUnmapped(*(pin++));

				Binvspacevalid = true;
			}

			bptr = Binvspace.begin();
		}
		else
		{
			bptr = bbaseptr;
		}
		
		#if 0
		int32_t const bstartforw = OVL.isInverse() ? (readsMeta2->at(OVL.bread).rlen - OVL.path.bepos) : OVL.path.bbpos;
		#endif
		int32_t const blen   = (OVL.path.bepos - OVL.path.bbpos);
		int32_t const bclipleft = OVL.path.bbpos;
		int32_t const bclipright = readsMeta2.at(OVL.bread).rlen - blen - bclipleft;

		// current point on A
		int32_t a_i = ( OVL.path.abpos / tspace ) * tspace;
		// current point on B
		int32_t b_i = ( OVL.path.bbpos );
		
		// reset trace container
		ATC.reset();
		
		for ( size_t i = 0; i < OVL.path.path.size(); ++i )
		{
			// block end point on A
			int32_t const a_i_1 = std::min ( static_cast<int32_t>(a_i + tspace), static_cast<int32_t>(OVL.path.aepos) );
			// block end point on B
			int32_t const b_i_1 = b_i + OVL.path.path[i].second;

			// block on A
			char const * asubsub_b = aptr + std::max(a_i,OVL.path.abpos);
			char const * asubsub_e = asubsub_b + a_i_1-std::max(a_i,OVL.path.abpos);
			
			// block on B
			char const * bsubsub_b = bptr + b_i;
			char const * bsubsub_e = bsubsub_b + (b_i_1-b_i);

			bool const ok = ND.process(asubsub_b,(asubsub_e-asubsub_b),bsubsub_b,bsubsub_e-bsubsub_b);
			assert ( ok );

			#if 0
			ND.printAlignmentLines(std::cout,asubsub_b,asubsub_e-asubsub_b,bsubsub_b,bsubsub_e-bsubsub_b,cols);
			#endif
			
			#if 0
			std::cerr << ED_EDR << "\t" << OVL.path.path[i].first << std::endl;
			ED.printAlignmentLines(std::cout,asubsub,bsubsub,cols);
			#endif
			
			// add trace to full alignment
			ATC.push(ND);
			
			// update start points
			b_i = b_i_1;
			a_i = a_i_1;
		}
		
		size_t const nopblocks = ATC.getOpBlocks(Aopblocks);
		size_t const cigopblocks = nopblocks + (bclipleft?1:0) + (bclipright?1:0);
		
		if ( cigopblocks > Acigop.size() )
			Acigop.resize(cigopblocks);
		uint64_t cigp = 0;
		
		int64_t as = 0;
		if ( bclipleft )
		{
			if ( primary )
				Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipleft);
			else
			{
				switch ( supplementary_seq_strategy )
				{
					case supplementary_seq_strategy_soft:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipleft);
						break;
					case supplementary_seq_strategy_hard:
					case supplementary_seq_strategy_none:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP,bclipleft);
						break;
				}
			}
		}
		for ( uint64_t i = 0; i < nopblocks; ++i )
		{
			switch ( Aopblocks[i].first )
			{
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL,Aopblocks[i].second);
					as += Aopblocks[i].second;
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF,Aopblocks[i].second);
					as -= static_cast<int64_t>(Aopblocks[i].second);
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS,Aopblocks[i].second);
					as -= static_cast<int64_t>(Aopblocks[i].second);
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL,Aopblocks[i].second);
					as -= static_cast<int64_t>(Aopblocks[i].second);
					break;
			}
		}
		if ( bclipright )
		{
			if ( primary )
				Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipright);
			else
			{
				switch ( supplementary_seq_strategy )
				{
					case supplementary_seq_strategy_soft:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipright);
						break;
					case supplementary_seq_strategy_hard:
					case supplementary_seq_strategy_none:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP,bclipright);
						break;
				}
			}
		}
		
		uint64_t const offa = tbuffer.end()-tbuffer.begin();

		std::string const readname = DB2.getReadName(OVL.bread,readsMeta2);
		uint64_t const readlen = readsMeta2.at(OVL.bread).rlen;
		
		if ( readlen > ASQ.size() )
		{
			ASQ.resize(readlen);
			std::fill(ASQ.begin(),ASQ.end(),255);
		}
		
		uint64_t storeseqlen;
		
		if ( primary )
		{
			storeseqlen = readlen;
		}
		else
		{
			switch ( supplementary_seq_strategy )
			{
				case supplementary_seq_strategy_hard:
					storeseqlen = readlen - (bclipleft+bclipright);
					break;
				case supplementary_seq_strategy_none:
					storeseqlen = 0;
					break;
				case supplementary_seq_strategy_soft:
				default:
					storeseqlen = readlen;
					break;
			}
		}

		libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment(
			tbuffer,
			seqenc,
			readname.begin(),
			readname.size(),
			OVL.aread, // ref id
			OVL.path.abpos, // pos
			255, // mapq (none given)
			(OVL.isInverse() ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0)
			|
			(primary ? 0 : libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSUPPLEMENTARY)
			,
			Acigop.begin(),
			cigp,
			-1, // next ref id
			-1, // next pos
			0, // template length (not given)
			primary ? bptr : (bptr + ((supplementary_seq_strategy==supplementary_seq_strategy_hard) ? bclipleft : 0)),
			storeseqlen,
			ASQ.begin(), /* quality (not given) */
			0, /* quality offset */
			false /* reset buffer */
		);

		if ( calmdnm )
		{		
			if ( primary || (supplementary_seq_strategy != supplementary_seq_strategy_none) )
			{
				libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(
					tbuffer.begin()+offa,tbuffer.end()-(tbuffer.begin()+offa),
					context,aptr + OVL.path.abpos,false
				);
			}
			else
			{
				sbuffer.reset();
				libmaus2::bambam::BamAlignmentEncoderBase::encodeSeq(sbuffer,seqenc,bptr+bclipleft,readlen-(bclipleft+bclipright));

				libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(
					tbuffer.begin()+offa,tbuffer.end()-(tbuffer.begin()+offa),
					context,aptr + OVL.path.abpos,
					sbuffer.begin(),
					readlen-(bclipleft+bclipright),
					false
				);
			}
			libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(tbuffer,"MD",context.md.get());
			libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(tbuffer,"NM",'i',context.nm);
			libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(tbuffer,"AS",'i',as);
		}

		uint64_t const offe = tbuffer.end()-tbuffer.begin();
		
		#if defined(LASTOBAM_ALIGNMENT_PRINT_DEBUG)
		// print alignment if requested
		if ( printAlignments )
		{
			libmaus2::lcs::AlignmentStatistics const AS = ATC.getAlignmentStatistics();
			
			if ( AS.getErrorRate() <= eratelimit )
			{
				std::cout << OVL << std::endl;
				std::cout << AS << std::endl;
				ATC.printAlignmentLines(std::cout,aptr + OVL.path.abpos,OVL.path.aepos-OVL.path.abpos,bptr + OVL.path.bbpos,OVL.path.bepos-OVL.path.bbpos,80);
			}
		}
		#endif
		
		return std::pair<uint64_t,uint64_t>(offa,offe);
	}
};

struct OverlapLoader : public libmaus2::parallel::PosixThread
{
	libmaus2::dazzler::align::AlignmentFile & algn;
	libmaus2::aio::InputStream & algnfile;
	
	typedef libmaus2::dazzler::align::Overlap o_type;
	typedef std::vector<o_type> v_type;
	typedef libmaus2::util::shared_ptr<v_type>::type p_type;
	
	libmaus2::parallel::SynchronousQueue             < p_type > freelist;
	libmaus2::parallel::TerminatableSynchronousQueue < p_type > filledlist;
	
	std::stack < libmaus2::dazzler::align::Overlap > stallstack;
	
	uint64_t const thres;
	
	OverlapLoader(
		libmaus2::dazzler::align::AlignmentFile & ralgn, libmaus2::aio::InputStream & ralgnfile, uint64_t const rthres = 32*1024
	)
	: algn(ralgn), algnfile(ralgnfile), thres(rthres)
	{
		for ( unsigned int i = 0; i < 4; ++i )
		{
			p_type p(new v_type);
			freelist.enque(p);
		}
	}
	
	virtual void * run()
	{
		libmaus2::dazzler::align::Overlap OVL;
		bool running = true;

		while ( running )
		{
			p_type V = freelist.deque();
			V->clear();
			libmaus2::dazzler::align::Overlap OVL;
			
			while ( 
				algn.peekNextOverlap(algnfile,OVL) &&
				(
					(!V->size()) 
					||
					(
						(V->size() < thres) &&
						(OVL.aread == V->front().aread)
					)
				)
			)
			{
				bool const ok = algn.getNextOverlap(algnfile,OVL);
				assert ( ok );
				V->push_back(OVL);
			}
			
			while ( 
				V->size() && algn.peekNextOverlap(algnfile,OVL) && (OVL.aread == V->front().aread) &&
				(OVL.bread == V->back().bread)	
			)
			{			
				bool const ok = algn.getNextOverlap(algnfile,OVL);
				assert ( ok );
				V->push_back(OVL);
			}
			
			if ( V->size() )
			{
				filledlist.enque(V);
			}
			else
			{
				freelist.enque(V);
				filledlist.terminate();
				running = false;
			}
		}
	
		return NULL;
	}
};

std::string getDefaultSupStoreStrat()
{
	return "none";
}

int lasToBAM(libmaus2::util::ArgInfo const & arginfo)
{
	#if defined(LASTOBAM_ALIGNMENT_PRINT_DEBUG)
	bool const printAlignments = arginfo.getValue<int>("print",1);
	#endif
	bool const loadall  = arginfo.getValue<int>("loadalla",false);
	bool loadalla = arginfo.getValue<int>("loadalla",loadall);
	bool loadallb = arginfo.getValue<int>("loadallb",loadall);
	#if 0
	double const eratelimit = arginfo.getValue<double>("eratelimit",1.0);
	#endif
	
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const calmdnm = arginfo.getValue<int>("calmdnm",getDefaultCalMdNm());

	uint64_t const numprocs = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	uint64_t const numthreads = std::max ( static_cast<uint64_t>(1), arginfo.getValueUnsignedNumeric<uint64_t>("threads",numprocs));
	
	std::string const supstorestrat_s = arginfo.getUnparsedValue("supstrat",getDefaultSupStoreStrat());
	
	LASToBAMConverter::supplementary_seq_strategy_t supstorestrat;
	
	if ( supstorestrat_s == "soft" )
		supstorestrat = LASToBAMConverter::supplementary_seq_strategy_soft;
	else if ( supstorestrat_s == "hard" )
		supstorestrat = LASToBAMConverter::supplementary_seq_strategy_hard;
	else if ( supstorestrat_s == "none" )
		supstorestrat = LASToBAMConverter::supplementary_seq_strategy_none;
	else
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "unknown strategy for storing supplementary alignments: " << supstorestrat_s << " (available options are soft, hard and none)" << std::endl;
		lme.finish();
		throw lme;	
	}
		
	std::string const referencekey = "reference";
	if ( ! arginfo.hasArg(referencekey) )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "reference key is missing (set reference=ref.fasta on the command line)" << std::endl;
		lme.finish();
		throw lme;
	}
	
	std::string const reference = arginfo.getUnparsedValue(referencekey,"");

	libmaus2::fastx::FastAIndex::unique_ptr_type Prefindex(libmaus2::fastx::FastAIndex::load(reference+".fai"));
	libmaus2::fastx::FastAIndex const & refindex = *Prefindex;
	
	std::string const curdir = libmaus2::util::ArgInfo::getCurDir();
	std::string const absreference =  (reference.size() && reference[0] == '/') ? reference : (curdir+'/'+reference);
	std::ostringstream sqstream;

	for ( uint64_t i = 0; i < refindex.size(); ++i )
	{
		libmaus2::fastx::FastAIndexEntry const & entry = refindex[i];
		sqstream << "@SQ\t" << "SN:" << entry.name << "\tLN:" << entry.length << "\tUR:file://" << absreference << std::endl;
	}
	
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB1;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB2;
	libmaus2::dazzler::db::DatabaseFile * DB1 = 0;
	libmaus2::dazzler::db::DatabaseFile * DB2 = 0;
	std::string aligns;
	std::vector<libmaus2::dazzler::db::Read> VreadsMeta1;
	std::vector<libmaus2::dazzler::db::Read> VreadsMeta2;
	std::vector<libmaus2::dazzler::db::Read> const * readsMeta1 = 0;
	std::vector<libmaus2::dazzler::db::Read> const * readsMeta2 = 0;

	libmaus2::autoarray::AutoArray<char> AreadsA;
	std::vector<uint64_t> AreadsOffA;	
	libmaus2::autoarray::AutoArray<char> AreadsB;
	std::vector<uint64_t> AreadsOffB;	

	libmaus2::autoarray::AutoArray<char> const * readsA = 0;
	std::vector<uint64_t> const * readsOffA = 0;	
	libmaus2::autoarray::AutoArray<char> const * readsB = 0;
	std::vector<uint64_t> const * readsOffB = 0;	
	
	if ( arginfo.restargs.size() == 2 )
	{
		std::string const dbfn = arginfo.restargs.at(0);
		aligns = arginfo.restargs.at(1);
	
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbfn));
		PDB->computeTrimVector();
		PDB1 = UNIQUE_PTR_MOVE(PDB);	
		DB1 = PDB1.get();
		DB2 = PDB1.get();
		
		DB1->getAllReads(VreadsMeta1);
		readsMeta1 = &VreadsMeta1;
		readsMeta2 = &VreadsMeta1;
		
		if ( loadalla || loadallb )
		{
			DB1->decodeAllReads(AreadsA,AreadsOffA);
			loadalla = true;
			loadallb = true;
			
			readsA = &AreadsA;
			readsOffA = &AreadsOffA;
			readsB = &AreadsA;
			readsOffB = &AreadsOffA;
		}
	}
	else if ( arginfo.restargs.size() == 3 )
	{
		std::string const db1fn = arginfo.restargs.at(0);
		std::string const db2fn = arginfo.restargs.at(1);
		aligns = arginfo.restargs.at(2);
	
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type TPDB1(new libmaus2::dazzler::db::DatabaseFile(db1fn));
		TPDB1->computeTrimVector();
		PDB1 = UNIQUE_PTR_MOVE(TPDB1);

		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type TPDB2(new libmaus2::dazzler::db::DatabaseFile(db2fn));
		TPDB2->computeTrimVector();
		PDB2 = UNIQUE_PTR_MOVE(TPDB2);
			
		DB1 = PDB1.get();
		DB2 = PDB2.get();
		
		DB1->getAllReads(VreadsMeta1);
		readsMeta1 = &VreadsMeta1;

		DB2->getAllReads(VreadsMeta2);
		readsMeta2 = &VreadsMeta2;

		if ( loadalla )
		{
			DB1->decodeAllReads(AreadsA,AreadsOffA);
			readsA = &AreadsA;
			readsOffA = &AreadsOffA;
		}
		if ( loadallb )
		{
			DB2->decodeAllReads(AreadsB,AreadsOffB);
			readsB = &AreadsB;
			readsOffB = &AreadsOffB;
		}
	}
	else
	{
		std::cerr << getUsage(arginfo) << std::endl;
		return EXIT_FAILURE;
	}
				
	libmaus2::aio::InputStream::unique_ptr_type Palgnfile(libmaus2::aio::InputStreamFactoryContainer::constructUnique(aligns));
	libmaus2::dazzler::align::AlignmentFile algn(*Palgnfile);

	libmaus2::lcs::EditDistanceTraceContainer ATC;		
	// libmaus2::lcs::NDextendDNA ND;
	libmaus2::lcs::ND ND;
	//libmaus2::dazzler::align::Overlap OVL;
	
	// number of alignments processed
	libmaus2::autoarray::AutoArray<char> Aspace;
	int64_t aid = -1;
	char const * aptr = NULL;
	
	#if 0
	libmaus2::autoarray::AutoArray<char> Bspace;
	libmaus2::autoarray::AutoArray<char> Binvspace;
	char const * bbaseptr = NULL;
	char const * bptr = NULL;
	bool Binvspacevalid = false;
	#endif

	libmaus2::aio::InputStream::unique_ptr_type PbaseStreamA(DB1->openBaseStream());
	libmaus2::aio::InputStream::unique_ptr_type PbaseStreamB(DB2->openBaseStream());

	// construct new header
	RgInfo const rginfo(arginfo);
	std::string const rgid = rginfo.ID; // ;arginfo.getUnparsedValue("RG","");

	std::ostringstream headerostr;
	headerostr << "@HD\tVN:1.4\tSO:unknown\n";
	headerostr 
		<< "@PG"<< "\t" 
		<< "ID:" << "lasToBAM" << "\t" 
		<< "PN:" << "lasToBAM" << "\t"
		<< "CL:" << arginfo.commandline << "\t"
		<< "VN:" << std::string(PACKAGE_VERSION)
		<< std::endl;
	headerostr << rginfo.toString();
	headerostr << sqstream.str();
	::libmaus2::bambam::BamHeader bamheader(headerostr.str());
	
	/*
	 * start md5 callbacks
	 */
	std::string md5filename;

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
	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	// construct writer
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pout ( 
		libmaus2::bambam::BamBlockWriterBaseFactory::construct(bamheader, arginfo, Pcbs) 
	);

	libmaus2::bambam::BamBlockWriterBase & bamwr = *Pout;

	::libmaus2::fastx::UCharBuffer ubuffer;
	::libmaus2::bambam::MdStringComputationContext context;
	::libmaus2::bambam::BamSeqEncodeTable const seqenc;
	::libmaus2::timing::RealTimeClock rtc; rtc.start();

	libmaus2::autoarray::AutoArray < LASToBAMConverter::unique_ptr_type > Aconverters(numthreads);
	
	for ( uint64_t i = 0; i < Aconverters.size(); ++i )
	{	
		libmaus2::aio::InputStream::unique_ptr_type TPbaseStreamB(DB2->openBaseStream());
		LASToBAMConverter::unique_ptr_type conv(new LASToBAMConverter(
			*DB2,
			*readsMeta1,
			*readsMeta2,
			TPbaseStreamB,
			loadallb,
			readsB,
			readsOffB,
			algn.tspace,
			calmdnm,
			supstorestrat
		));
		Aconverters[i] = UNIQUE_PTR_MOVE(conv);
	}

	
	uint64_t processed = 0;
	
	OverlapLoader OL(algn,*Palgnfile,32*1024);
	OL.start();
	bool reading = true;
	libmaus2::bitio::BitVector BV(readsMeta2->size());
	
	while ( reading )
	{
		libmaus2::util::shared_ptr < std::vector<libmaus2::dazzler::align::Overlap> >::type PVOVL;
		
		try
		{
			PVOVL = OL.filledlist.deque();
			std::vector<libmaus2::dazzler::align::Overlap> & VOVL = *PVOVL;

			if ( VOVL.size() && VOVL.front().aread != aid )
			{
				libmaus2::dazzler::db::Read const & R = readsMeta1->at(VOVL.front().aread);
				
				if ( loadalla )
				{
					#if 0
					if ( R.rlen > static_cast<int64_t>(Aspace.size()) )
						Aspace.resize(R.rlen);
					std::copy(readsA.begin()+readsOffA[VOFL.front().aread],readsA.begin()+readsOffA[VOFL.front().aread+1],Aspace.begin());
					#endif
					
					aptr = readsA->begin()+(*readsOffA)[VOVL.front().aread];
				}
				else
				{
					PbaseStreamA->clear();
					PbaseStreamA->seekg(R.boff);
					libmaus2::dazzler::db::DatabaseFile::decodeRead(*PbaseStreamA,Aspace,R.rlen);
					aptr = Aspace.begin();
				}
				
				aid = VOVL.front().aread;
			}

			for ( uint64_t i = 0; i < Aconverters.size(); ++i )
			{
				Aconverters[i]->ubuffer.reset();
				Aconverters[i]->setaptr(aptr);
				Aconverters[i]->wbuffer.reset();
			}
				
			std::vector < std::pair < uint64_t, std::pair<uint64_t,uint64_t> > > off(VOVL.size());
				
			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
			#endif
			for ( uint64_t i = 0; i < VOVL.size(); ++i )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 1;
				#endif

				LASToBAMConverter & converter = *(Aconverters[tid]);
				bool const primary = 
					(!BV.get(VOVL[i].bread))
					&&
					(
						(i == 0) 
						||
						(VOVL[i].bread != VOVL[i-1].bread)
					)
				;
				std::pair<uint64_t,uint64_t> const loff = converter(VOVL[i],converter.ubuffer,primary);	
				off[i] = std::pair < uint64_t, std::pair<uint64_t,uint64_t> >(tid,loff);
				if ( primary )
					BV.setSync(VOVL[i].bread);
			}

			#if 0
			uint64_t const onumthreads = libmaus2::parallel::OMPNumThreadsScope::getMaxThreads();
			uint64_t const operthread = (off.size() + onumthreads - 1) / onumthreads;
			
			for ( uint64_t t = 0; t < onumthreads; ++t )
			{
				uint64_t const low = std::min(t * operthread, static_cast<uint64_t>(off.size()));
				uint64_t const high = std::min(low + operthread, static_cast<uint64_t>(off.size()));
				::libmaus2::fastx::UCharBuffer const & wbuffer = *(Aconverters[t]).wbuffer;
				
				for ( uint64_t i = low; i < high; ++i )
				{
					uint64_t const tid = off[i].first;
					std::pair<uint64_t,uint64_t> const & loff = off[i].second;
					LASToBAMConverter & converter = *(Aconverters[tid]);
					::libmaus2::fastx::UCharBuffer const & ubuffer = converter.ubuffer;

					wbuffer.put(ubuffer.begin() + loff.first,loff.second-loff.first);
				}
			}
			#endif

			for ( uint64_t i = 0; i < off.size(); ++i )
			{
				uint64_t const tid = off[i].first;
				std::pair<uint64_t,uint64_t> const & loff = off[i].second;
				LASToBAMConverter & converter = *(Aconverters[tid]);
				::libmaus2::fastx::UCharBuffer const & ubuffer = converter.ubuffer;
				bamwr.writeBamBlock(ubuffer.begin() + loff.first,loff.second-loff.first);
			}
			
			processed += VOVL.size();

			if ( verbose )
				std::cerr << "[V] " << processed << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " " << processed/rtc.getElapsedSeconds() << std::endl;

			VOVL.clear();
			OL.freelist.enque(PVOVL);
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
			reading = false;
		}
	}                                
		
	Pout.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}

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

				std::cerr << getUsage(arginfo) << std::endl << std::endl;

				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus2::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus2::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("calmdnm=<[")+::biobambam2::Licensing::formatNumber(getDefaultCalMdNm())+"]>", "calculate MD and NM aux fields (for coordinate sorted output only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("supstrat=<[")+getDefaultSupStoreStrat()+"]>", "query string storage strategy for supplementary reads (soft, hard or none, default: none)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return lasToBAM(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
