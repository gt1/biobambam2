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
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/aio/OutputStreamInstance.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

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

static char const padsym = '*';

struct ConsensusAccuracy;
ConsensusAccuracy operator+(ConsensusAccuracy const & A, ConsensusAccuracy const & B);

struct ConsensusAccuracy
{
	uint64_t matches;
	uint64_t mismatches;
	uint64_t insertions;
	uint64_t deletions;
	uint64_t refbasesprocessed;
	uint64_t refbasesexpected;
	libmaus2::util::Histogram depthhistogram;
	
	ConsensusAccuracy() : matches(0), mismatches(0), insertions(0), deletions(0), refbasesprocessed(0), refbasesexpected(0) {}
	ConsensusAccuracy(
		ConsensusAccuracy const & o
	) : matches(o.matches), mismatches(o.mismatches), insertions(o.insertions), deletions(o.deletions), refbasesprocessed(o.refbasesprocessed),
	    refbasesexpected(o.refbasesexpected)
	{
	
	}
	ConsensusAccuracy(uint64_t const & rrefbasesexpected) 
	: matches(0), mismatches(0), insertions(0), deletions(0), refbasesprocessed(0), refbasesexpected(rrefbasesexpected) {}
	
	ConsensusAccuracy & operator+=(ConsensusAccuracy const & O)
	{
		matches += O.matches;
		mismatches += O.mismatches;
		insertions += O.insertions;
		deletions += O.deletions;
		refbasesprocessed += O.refbasesprocessed;
		refbasesexpected += O.refbasesexpected;
		depthhistogram.merge(O.depthhistogram);
		return *this;
	}
};

ConsensusAccuracy operator+(ConsensusAccuracy const & A, ConsensusAccuracy const & B)
{
	ConsensusAccuracy C(A);
	C += B;
	return C;
}

std::ostream & operator<<(std::ostream & out, ConsensusAccuracy const & C)
{
	return out << "ConsensusAccuracy("
		<< "matches=" << C.matches << ","
		<< "mismatches=" << C.mismatches << ","
		<< "insertions=" << C.insertions << ","
		<< "deletions=" << C.deletions << ","
		<< "refbasesprocessed=" << C.refbasesprocessed << ","
		<< "refbasesexpected=" << C.refbasesexpected << ","
		<< "refmissing=" << static_cast<int64_t>(C.refbasesexpected)-static_cast<int64_t>(C.refbasesprocessed) << ","
		<< "errorrate=" << 
			(static_cast<double>(C.mismatches + C.insertions + C.deletions + static_cast<int64_t>(C.refbasesexpected)-static_cast<int64_t>(C.refbasesprocessed))
			/
			(C.matches + C.mismatches + C.insertions + C.deletions + static_cast<int64_t>(C.refbasesexpected)-static_cast<int64_t>(C.refbasesprocessed)))
		<< ")";
}

struct ConsensusAux
{
	libmaus2::autoarray::AutoArray<uint64_t> M;
	libmaus2::autoarray::AutoArray<uint64_t> C;

	ConsensusAux() : M(256), C(256)
	{
		std::fill(M.begin(),M.end(),1);
		std::fill(C.begin(),C.end(),0);
	}
};

struct HeapEntry
{
	// bases inserted before reference base
	std::vector< std::vector< std::pair<char,uint8_t> > > I;
	// bases aligned to reference base
	std::vector< std::pair<char,uint8_t> > V;
	uint64_t iadd;
	
	HeapEntry() : I(), V(), iadd(0)
	{
	}
	
	static uint8_t getConsensusBase(ConsensusAux const & caux) 
	{
		uint64_t maxval = 0;
		uint64_t maxindex = 0;
		uint64_t maxcnt = 0;
		
		for ( uint64_t i = 0; i < caux.C.size(); ++i )
			if ( caux.C[i] > maxval )
			{
				maxval = caux.C[i];
				maxindex = i;
				maxcnt = 1;
			}
			else if ( caux.C[i] == maxval )
			{
				maxcnt++;
			}
			
		if ( maxcnt == 1 )
			return maxindex;
		else
			return 'N';
	}
	
	std::ostream & toStream(
		std::ostream & out, 
		uint64_t refpos, 
		std::string const & refidname, 
		int const refbase,
		ConsensusAux & caux,
		ConsensusAccuracy * consacc = 0,
		std::ostream * consostr = 0
	)
	{
		// insertions
		if ( I.size() )
		{
			uint64_t const cov = V.size() + iadd;
			
			uint64_t maxi = 0;
			for ( uint64_t i = 0; i < I.size(); ++i )
				maxi = std::max(maxi,static_cast<uint64_t>(I[i].size()));
			
			// vector of active indices	
			std::vector<uint64_t> active(I.size());
			for ( uint64_t i = 0; i < active.size(); ++i )
				active[i] = i;
			
			// vector for symbol sorting	
			std::vector< std::pair<char,uint8_t> > sortvec(active.size());
			for ( uint64_t j = 0; active.size(); ++j )
			{
				// number still active in next round
				uint64_t o = 0;
				sortvec.resize(active.size());
				for ( uint64_t i = 0; i < active.size(); ++i )
				{
					sortvec[i] = I[active[i]][j];
					
					if ( j+1 < I[active[i]].size() )
						active[o++] = active[i];
				}
				active.resize(o);
					
				std::sort(sortvec.begin(),sortvec.end());

				out << refidname << '\t' << refpos << '\t' << (-static_cast<int64_t>(maxi))+static_cast<int64_t>(j) << '\t';

				if ( cov > sortvec.size() )
					for ( uint64_t i = 0; i < cov-sortvec.size(); ++i )
					{
						out.put(padsym);
						caux.C[padsym] += caux.M[padsym];
					}
				for ( uint64_t i = 0; i < sortvec.size(); ++i )
				{
					out.put(sortvec[i].first);
					caux.C[sortvec[i].first] += caux.M[sortvec[i].first];
				}
				
				out.put('\t');
				
				uint8_t const consbase = getConsensusBase(caux);
				
				out.put(consbase);
				
				if ( consostr && consbase != padsym )
					consostr->put(consbase);
				
				if ( consacc && consbase != padsym )
					consacc->insertions++;

				out.put('\n');

				if ( cov > sortvec.size() )
					for ( uint64_t i = 0; i < cov-sortvec.size(); ++i )
						caux.C[padsym] -= caux.M[padsym];
				for ( uint64_t i = 0; i < sortvec.size(); ++i )
					caux.C[sortvec[i].first] -= caux.M[sortvec[i].first];
			}
		}
		
		std::sort(V.begin(),V.end());

		out << refidname << '\t' << refpos << '\t' << 0 << '\t';
		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			out.put(V[i].first);
			caux.C[V[i].first] += caux.M[V[i].first];
		}
			
		out << "\t";
		
		if ( refbase != -1 )
			out.put(refbase);

		out.put('\t');

		uint8_t const consbase = getConsensusBase(caux);
		
		if ( consacc && (consbase == padsym) )
		{
			consacc->deletions++;
			consacc->depthhistogram(0);
		}
		if ( refbase != -1 && consacc && (consbase != padsym) )
		{
			if ( consbase == refbase 
				&& 
				(
					consbase == 'a' || consbase == 'A' ||
					consbase == 'c' || consbase == 'C' ||
					consbase == 'g' || consbase == 'G' ||
					consbase == 't' || consbase == 'T'
				)
			)
				consacc->matches++;
			else
				consacc->mismatches++;

			consacc->depthhistogram(V.size());
		}
		
		out.put(consbase);

		if ( consostr && consbase != padsym )
			consostr->put(consbase);
			
		out.put('\n');

		for ( uint64_t i = 0; i < V.size(); ++i )
			caux.C[V[i].first] -= caux.M[V[i].first];
		
		for ( uint64_t i = 0; i < caux.C.size(); ++i )
			assert ( caux.C[i] == 0 );
			
		if ( consacc )
			consacc->refbasesprocessed++;
			
		return out;
	}
};

int bamheap2(libmaus2::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue("verbose",getDefaultVerbose());
	std::string const reference = arginfo.getUnparsedValue("reference",std::string());
	std::string const outputprefix = arginfo.getUnparsedValue("outputprefix",std::string());
	
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus2::bambam::BamHeader const & header = dec.getHeader();	
	::libmaus2::bambam::BamAlignment const & algn = dec.getAlignment();
	
	double const damult = arginfo.getValue<double>("amult",1);
	double const dcmult = arginfo.getValue<double>("cmult",1);
	double const dgmult = arginfo.getValue<double>("gmult",1);
	double const dtmult = arginfo.getValue<double>("tmult",1);
	double const dpadmult = arginfo.getValue<double>("padmult",1);
	
	double maxmult = 0;
	maxmult = std::max(damult,maxmult);
	maxmult = std::max(dcmult,maxmult);
	maxmult = std::max(dgmult,maxmult);
	maxmult = std::max(dtmult,maxmult);
	maxmult = std::max(dpadmult,maxmult);
	
	uint64_t const amult = std::floor((damult / maxmult) * (1ull<<16) + 0.5);
	uint64_t const cmult = std::floor((dcmult / maxmult) * (1ull<<16) + 0.5);
	uint64_t const gmult = std::floor((dgmult / maxmult) * (1ull<<16) + 0.5);
	uint64_t const tmult = std::floor((dtmult / maxmult) * (1ull<<16) + 0.5);
	uint64_t const padmult = std::floor((dpadmult / maxmult) * (1ull<<16) + 0.5);
	
	libmaus2::fastx::FastAIndex::unique_ptr_type Pindex;
	libmaus2::aio::InputStreamInstance::unique_ptr_type PCIS;
	if ( reference.size() )
	{
		libmaus2::fastx::FastAIndex::unique_ptr_type Tindex(
			libmaus2::fastx::FastAIndex::load(reference+".fai")
		);
		Pindex = UNIQUE_PTR_MOVE(Tindex);
		
		libmaus2::aio::InputStreamInstance::unique_ptr_type TCIS(new libmaus2::aio::InputStreamInstance(reference));
		PCIS = UNIQUE_PTR_MOVE(TCIS);
	}

	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;
	libmaus2::autoarray::AutoArray<char> bases;
	
	int64_t prevrefid = -1;
	std::string refidname = "*";
	
	std::map< uint64_t, HeapEntry > M;
	uint64_t alcnt = 0;
	std::vector< std::pair<char,uint8_t> > pendinginserts;
	int64_t loadedRefId = -1;
	int64_t streamRefId = -1;
	libmaus2::autoarray::AutoArray<char> refseqbases;
	ConsensusAccuracy * consacc = 0;
	std::map<uint64_t,ConsensusAccuracy> Mconsacc;
	typedef libmaus2::util::shared_ptr<std::ostringstream>::type stream_ptr_type;
	stream_ptr_type Pstream;
	ConsensusAux Caux;
	
	Caux.M['a'] = Caux.M['A'] = amult;
	Caux.M['c'] = Caux.M['C'] = cmult;
	Caux.M['g'] = Caux.M['G'] = gmult;
	Caux.M['t'] = Caux.M['T'] = tmult;
	Caux.M[padsym] = padmult;
	
	while ( dec.readAlignment() )
	{
		if ( algn.isMapped() && (!algn.isQCFail()) )
		{
			assert ( ! pendinginserts.size() );
		
			uint32_t const numcigop = algn.getCigarOperations(cigop);
			uint64_t readpos = 0;
			uint64_t refpos = algn.getPos();
			uint64_t const seqlen = algn.decodeRead(bases);
			uint8_t const * qual = libmaus2::bambam::BamAlignmentDecoderBase::getQual(algn.D.begin());
			
			// handle finished columns
			if ( algn.getRefID() != prevrefid )
			{
				while ( M.size() )
				{
					HeapEntry & H = M.begin()->second;
					
					if ( outputprefix.size() && (streamRefId != prevrefid) )
					{
						if ( Pstream )
						{
							std::ostringstream fnostr;
							fnostr << outputprefix << "_" << header.getRefIDName(streamRefId);
							libmaus2::aio::OutputStreamInstance PFOS(fnostr.str());
							PFOS << ">" << header.getRefIDName(streamRefId) << '\n';
							PFOS << Pstream->str() << '\n';
							
							Pstream.reset();
						}
						
						stream_ptr_type Tstream(new std::ostringstream);
						Pstream = Tstream;
						streamRefId = prevrefid;
					}
					
					if ( Pindex && (loadedRefId != prevrefid) )
					{
						refseqbases = Pindex->readSequence(*PCIS, Pindex->getSequenceIdByName(refidname));
						loadedRefId = prevrefid;
						
						if ( Mconsacc.find(loadedRefId) == Mconsacc.end() )
							Mconsacc[loadedRefId] = ConsensusAccuracy(refseqbases.size());
						
						consacc = &(Mconsacc[loadedRefId]);
					}
					
					H.toStream(std::cout,M.begin()->first,refidname,(M.begin()->first < refseqbases.size()) ? static_cast<int>(refseqbases[M.begin()->first]) : -1,Caux,consacc,Pstream.get());
					
					M.erase(M.begin());
				}
			
				prevrefid = algn.getRefID();
				refidname = header.getRefIDName(prevrefid);
			}
			else
			{
				while ( M.size() && M.begin()->first < refpos )
				{
					HeapEntry & H = M.begin()->second;

					if ( outputprefix.size() && (streamRefId != prevrefid) )
					{
						if ( Pstream )
						{
							std::ostringstream fnostr;
							fnostr << outputprefix << "_" << header.getRefIDName(streamRefId);
							libmaus2::aio::OutputStreamInstance PFOS(fnostr.str());
							PFOS << ">" << header.getRefIDName(streamRefId) << '\n';
							PFOS << Pstream->str() << '\n';

							Pstream.reset();
						}
						
						stream_ptr_type Tstream(new std::ostringstream);
						Pstream = Tstream;
						streamRefId = prevrefid;
					}

					if ( Pindex && (loadedRefId != prevrefid) )
					{
						refseqbases = Pindex->readSequence(*PCIS, Pindex->getSequenceIdByName(refidname));
						loadedRefId = prevrefid;

						if ( Mconsacc.find(loadedRefId) == Mconsacc.end() )
							Mconsacc[loadedRefId] = ConsensusAccuracy(refseqbases.size());

						consacc = &(Mconsacc[loadedRefId]);
					}
					
					H.toStream(std::cout,M.begin()->first,refidname,(M.begin()->first < refseqbases.size()) ? static_cast<int>(refseqbases[M.begin()->first]) : -1,Caux,consacc,Pstream.get());
					
					M.erase(M.begin());				
				}
			}
			
			for ( uint64_t ci = 0; ci < numcigop; ++ci )
			{
				uint64_t const ciglen = cigop[ci].second;
				
				switch ( cigop[ci].first )
				{
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH:
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
					{
						if ( pendinginserts.size() )
						{
							M[refpos].I.push_back(pendinginserts);
							pendinginserts.resize(0);
						}
					
						for ( uint64_t i = 0; i < ciglen; ++i )
						{
							M[refpos].V.push_back(std::make_pair(bases[readpos],qual[readpos]));
							readpos++;
							refpos++;
						}
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
					{
						for ( uint64_t i = 0; i < ciglen; ++i, ++readpos )
							pendinginserts.push_back(std::make_pair(bases[readpos],qual[readpos]));
						break;
					}
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
						// handle pending inserts
						if ( pendinginserts.size() )
						{
							M[refpos].I.push_back(pendinginserts);
							pendinginserts.resize(0);
						}
						
						// deleting bases from the reference
						for ( uint64_t i = 0; i < ciglen; ++i, ++refpos )
							M[refpos].V.push_back(std::make_pair(padsym,0));
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CREF_SKIP:
						// handle pending inserts
						if ( pendinginserts.size() )
						{
							M[refpos].I.push_back(pendinginserts);
							pendinginserts.resize(0);
						}

						// skip bases on reference
						for ( uint64_t i = 0; i < ciglen; ++i )
						{
							refpos++;
						}
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP:
						// skip bases on read
						for ( uint64_t i = 0; i < ciglen; ++i )
						{
							readpos++;
						}
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP:
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CPAD:
					{
						for ( uint64_t i = 0; i < ciglen; ++i, ++readpos )
							pendinginserts.push_back(std::make_pair(padsym,0));
						break;
					}
				}
			}

			if ( pendinginserts.size() )
			{
				M[refpos].I.push_back(pendinginserts);
				M[refpos].iadd++;
				pendinginserts.resize(0);
			}

			assert ( readpos == seqlen );
		}
		
		if ( verbose && ((++alcnt % (1024*1024)) == 0) )
			std::cerr << "[V] " << alcnt << std::endl;
	}

	while ( M.size() )
	{
		HeapEntry & H = M.begin()->second;

		if ( outputprefix.size() && (streamRefId != prevrefid) )
		{
			if ( Pstream )
			{
				std::ostringstream fnostr;
				fnostr << outputprefix << "_" << header.getRefIDName(streamRefId);
				libmaus2::aio::OutputStreamInstance PFOS(fnostr.str());
				PFOS << ">" << header.getRefIDName(streamRefId) << '\n';
				PFOS << Pstream->str() << '\n';

				Pstream.reset();
			}
			
			stream_ptr_type Tstream(new std::ostringstream);
			Pstream = Tstream;
			streamRefId = prevrefid;
		}

		if ( Pindex && (loadedRefId != prevrefid) )
		{
			refseqbases = Pindex->readSequence(*PCIS, Pindex->getSequenceIdByName(refidname));
			loadedRefId = prevrefid;

			if ( Mconsacc.find(loadedRefId) == Mconsacc.end() )
				Mconsacc[loadedRefId] = ConsensusAccuracy(refseqbases.size());

			consacc = &(Mconsacc[loadedRefId]);
		}
			
		H.toStream(std::cout,M.begin()->first,refidname,(M.begin()->first < refseqbases.size()) ? static_cast<int>(refseqbases[M.begin()->first]) : -1,Caux,consacc,Pstream.get());
		
		M.erase(M.begin());
	}
	
	if ( Pstream )
	{
		std::ostringstream fnostr;
		fnostr << outputprefix << "_" << header.getRefIDName(streamRefId);
		libmaus2::aio::OutputStreamInstance PFOS(fnostr.str());
		PFOS << ">" << header.getRefIDName(streamRefId) << '\n';
		PFOS << Pstream->str() << '\n';

		Pstream.reset();
	}
	
	ConsensusAccuracy constotal;
	for ( std::map<uint64_t,ConsensusAccuracy>::const_iterator ita = Mconsacc.begin(); ita != Mconsacc.end(); ++ita )
	{
		std::cerr << header.getRefIDName(ita->first) << "\t" << ita->second << std::endl;

		std::map<uint64_t,uint64_t> const M = ita->second.depthhistogram.get();
		uint64_t total = 0;
		uint64_t preavg = 0;
		for ( std::map<uint64_t,uint64_t>::const_iterator aita = M.begin(); aita != M.end(); ++aita )
		{
			total += aita->second;
			preavg += aita->first * aita->second;
		}

		uint64_t acc = 0;		
		for ( std::map<uint64_t,uint64_t>::const_iterator aita = M.begin(); aita != M.end(); ++aita )
		{
			acc += aita->second;
			std::cerr << "H[" << header.getRefIDName(ita->first) << "," << aita->first << ",+]"
				<< "\t" << aita->second << "\t" << static_cast<double>(aita->second)/total
				<< "\t" << acc << "\t" << static_cast<double>(acc)/total << std::endl;
		}
		acc = 0;
		for ( std::map<uint64_t,uint64_t>::const_reverse_iterator aita = M.rbegin(); aita != M.rend(); ++aita )
		{
			acc += aita->second;
			std::cerr << "H[" << header.getRefIDName(ita->first) << "," << aita->first << ",-]"
				<< "\t" << aita->second << "\t" << static_cast<double>(aita->second)/total
				<< "\t" << acc << "\t" << static_cast<double>(acc)/total << std::endl;
		}
		
		std::cerr << "H[" << header.getRefIDName(ita->first) << ",avg]\t" << 
			static_cast<double>(preavg)/total << std::endl;
		
		constotal += ita->second;
	}
	if ( Mconsacc.size() )
	{
		std::cerr << "all\t" << constotal << std::endl;

		std::map<uint64_t,uint64_t> const M = constotal.depthhistogram.get();
		uint64_t total = 0;
		uint64_t preavg = 0;
		for ( std::map<uint64_t,uint64_t>::const_iterator aita = M.begin(); aita != M.end(); ++aita )
		{
			total += aita->second;
			preavg += aita->first * aita->second;
		}

		uint64_t acc = 0;		
		for ( std::map<uint64_t,uint64_t>::const_iterator aita = M.begin(); aita != M.end(); ++aita )
		{
			acc += aita->second;
			std::cerr << "H[" << "all" << "," << aita->first << ",+]"
				<< "\t" << aita->second << "\t" << static_cast<double>(aita->second)/total
				<< "\t" << acc << "\t" << static_cast<double>(acc)/total << std::endl;
		}
		acc = 0;
		for ( std::map<uint64_t,uint64_t>::const_reverse_iterator aita = M.rbegin(); aita != M.rend(); ++aita )
		{
			acc += aita->second;
			std::cerr << "H[" << "all" << "," << aita->first << ",-]"
				<< "\t" << aita->second << "\t" << static_cast<double>(aita->second)/total
				<< "\t" << acc << "\t" << static_cast<double>(acc)/total << std::endl;
		}
		
		std::cerr << "H[all,avg]\t" << static_cast<double>(preavg) / total << std::endl;
		
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
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam2::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return bamheap2(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
