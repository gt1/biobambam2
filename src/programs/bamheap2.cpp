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
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

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
	
	std::ostream & toStream(std::ostream & out, uint64_t refpos, std::string const & refidname)
	{
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

				out << refidname << "\t" << refpos << "\t" << (-static_cast<int64_t>(maxi))+static_cast<int64_t>(j) << "\t";
				if ( cov > sortvec.size() )
					for ( uint64_t i = 0; i < cov-sortvec.size(); ++i )
						out << padsym;
				for ( uint64_t i = 0; i < sortvec.size(); ++i )
					out << sortvec[i].first;
				out << '\n';
			}
		}
		
		std::sort(V.begin(),V.end());
		
		out << refidname << "\t" << refpos << "\t" << 0 << "\t";
		for ( uint64_t i = 0; i < V.size(); ++i )
			out << V[i].first;
			
		out << '\n';
			
		return out;
	}
};

int bamheap2(libmaus::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue("verbose",getDefaultVerbose());
	
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamHeader const & header = dec.getHeader();	
	::libmaus::bambam::BamAlignment const & algn = dec.getAlignment();

	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	libmaus::autoarray::AutoArray<char> bases;
	
	int64_t prevrefid = -1;
	std::string refidname = "*";
	
	std::map< uint64_t, HeapEntry > M;
	uint64_t alcnt = 0;
	std::vector< std::pair<char,uint8_t> > pendinginserts;
	
	while ( dec.readAlignment() )
	{
		if ( algn.isMapped() && (!algn.isQCFail()) )
		{
			assert ( ! pendinginserts.size() );
		
			uint32_t const numcigop = algn.getCigarOperations(cigop);
			uint64_t readpos = 0;
			uint64_t refpos = algn.getPos();
			uint64_t const seqlen = algn.decodeRead(bases);
			uint8_t const * qual = libmaus::bambam::BamAlignmentDecoderBase::getQual(algn.D.begin());
			
			// handle finished columns
			if ( algn.getRefID() != prevrefid )
			{
				while ( M.size() )
				{
					HeapEntry & H = M.begin()->second;
					
					H.toStream(std::cout,M.begin()->first,refidname);
					
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
					
					H.toStream(std::cout,M.begin()->first,refidname);
					
					M.erase(M.begin());				
				}
			}
			
			for ( uint64_t ci = 0; ci < numcigop; ++ci )
			{
				uint64_t const ciglen = cigop[ci].second;
				
				switch ( cigop[ci].first )
				{
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
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
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
					{
						for ( uint64_t i = 0; i < ciglen; ++i, ++readpos )
							pendinginserts.push_back(std::make_pair(bases[readpos],qual[readpos]));
						break;
					}
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL:
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
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP:
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
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP:
						// skip bases on read
						for ( uint64_t i = 0; i < ciglen; ++i )
						{
							readpos++;
						}
						break;
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP:
						break;
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CPAD:
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
			
		H.toStream(std::cout,M.begin()->first,refidname);
		
		M.erase(M.begin());
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::util::ArgInfo const arginfo(argc,argv);

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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

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
