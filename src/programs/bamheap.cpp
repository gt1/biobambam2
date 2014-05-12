/*
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
*/
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/util/SimpleCountingHash.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/fastx/Phred.hpp>
#include <libmaus/math/binom.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 1; }
static std::string getDefaultInputFormat() { return "bam"; }

enum bam_heap_entry_base
{
	bam_heap_entry_base_A,
	bam_heap_entry_base_C,
	bam_heap_entry_base_G,
	bam_heap_entry_base_T,
	bam_heap_entry_base_N,
	bam_heap_entry_base_DEL
	// bam_heap_entry_base_IGN
};

std::ostream & operator<<(std::ostream & out, bam_heap_entry_base const base)
{
	switch ( base )
	{
		case bam_heap_entry_base_A: out << 'A'; break;
		case bam_heap_entry_base_C: out << 'C'; break;
		case bam_heap_entry_base_G: out << 'G'; break;
		case bam_heap_entry_base_T: out << 'T'; break;
		case bam_heap_entry_base_N: out << 'N'; break;
		case bam_heap_entry_base_DEL: out << '-'; break;
		// case bam_heap_entry_base_IGN: 
		default: break;
	}
	
	return out;
}

uint64_t strlen(std::pair<uint8_t,uint8_t> const * const rp)
{
	std::pair<uint8_t,uint8_t> const * p = rp;
	while ( p->first )
		++p;
	return p-rp;
}

struct BamHeapEntry
{
	typedef BamHeapEntry this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t refpos;
	std::vector< std::pair<bam_heap_entry_base,uint8_t> > baseheap;
	std::vector<std::pair<uint8_t,uint8_t> *> insbases;
	
	void reset(std::vector<std::vector<std::pair<uint8_t,uint8_t> *> > & pairstringfreelist)
	{
		baseheap.resize(0);
		
		for ( uint64_t i = 0; i < insbases.size(); ++i )
			pairstringfreelist[strlen(insbases[i])+1].push_back(insbases[i]);
		
		insbases.resize(0);
	}
};


struct EntriesContainer
{
	std::vector< BamHeapEntry::shared_ptr_type > entrypool;
	std::vector< BamHeapEntry * > entryfreelist;
	std::deque < BamHeapEntry * > entries;
	uint64_t entryptr;

	#if 0
	std::vector<std::vector<libmaus::autoarray::AutoArray<uint8_t>::shared_ptr_type > > stringpool;
	std::vector<std::vector<uint8_t *> > stringfreelist;
	#endif
	
	std::vector<std::vector<libmaus::autoarray::AutoArray<std::pair<uint8_t,uint8_t> >::shared_ptr_type > > pairstringpool;
	std::vector<std::vector<std::pair<uint8_t,uint8_t> *> > pairstringfreelist;
	
	EntriesContainer()
	{
	
	}
	
	#if 0
	uint8_t * getString(uint64_t const len)
	{
		uint64_t const ind = len+1;
		
		while ( ! ( ind < stringpool.size() ) )
		{
			stringpool.push_back(std::vector<libmaus::autoarray::AutoArray<uint8_t>::shared_ptr_type >());
			stringfreelist.push_back(std::vector<uint8_t *>());
		}
		
		assert ( ind < stringpool.size() );
		assert ( ind < stringfreelist.size() );
	
		if ( ! stringfreelist[ind].size() )
		{
			libmaus::autoarray::AutoArray<uint8_t>::shared_ptr_type tptr(
				new libmaus::autoarray::AutoArray<uint8_t>(ind,false)
			);
			stringpool[ind].push_back(tptr);
			stringfreelist[ind].push_back(stringpool[ind].back()->get());
		}
		
		uint8_t * str = stringfreelist[ind].back();
		stringfreelist[ind].pop_back();
		str[len] = 0;
				
		return str;
	}
	#endif

	std::pair<uint8_t,uint8_t> * getPairString(uint64_t const len)
	{
		uint64_t const ind = len+1;
		
		while ( ! ( ind < pairstringpool.size() ) )
		{
			pairstringpool.push_back(std::vector<libmaus::autoarray::AutoArray< std::pair<uint8_t,uint8_t> >::shared_ptr_type >());
			pairstringfreelist.push_back(std::vector<std::pair<uint8_t,uint8_t> *>());
		}
		
		assert ( ind < pairstringpool.size() );
		assert ( ind < pairstringfreelist.size() );
	
		if ( ! pairstringfreelist[ind].size() )
		{
			libmaus::autoarray::AutoArray<std::pair<uint8_t,uint8_t> >::shared_ptr_type tptr(
				new libmaus::autoarray::AutoArray< std::pair<uint8_t,uint8_t> >(ind,false)
			);
			pairstringpool[ind].push_back(tptr);
			pairstringfreelist[ind].push_back(pairstringpool[ind].back()->get());
		}
		
		std::pair<uint8_t,uint8_t> * pairstr = pairstringfreelist[ind].back();
		pairstringfreelist[ind].pop_back();
		pairstr[len] = std::pair<uint8_t,uint8_t>(0,0);
				
		return pairstr;
	}
	
	void resetPointer()
	{
		entryptr = 0;
	}
	
	BamHeapEntry * getEntry()
	{
		if ( ! entryfreelist.size() )
		{
			BamHeapEntry::shared_ptr_type tptr(new BamHeapEntry);
			entrypool.push_back(tptr);
			entryfreelist.push_back(entrypool.back().get());
		}
		
		assert ( entryfreelist.size() );

		BamHeapEntry * entry = entryfreelist.back();
		entryfreelist.pop_back();
		
		entry->reset(pairstringfreelist);
		
		return entry;
	}
	
	BamHeapEntry * getEntry(uint64_t refpos)
	{
		for ( uint64_t i = 1; i < entries.size(); ++i )
			assert ( entries[i-1]->refpos < entries[i]->refpos );
	
		while ( entryptr != entries.size() && entries[entryptr]->refpos < refpos )
			entryptr++;
			
		// entry is not there and is to be inserted at the end of the list
		if ( entryptr == entries.size() )
		{
			BamHeapEntry * entry = getEntry();
			entry->refpos = refpos;
			entries.push_back(entry);
			return entry;
		}
		// entry is already present
		else if ( entries[entryptr]->refpos == refpos )
		{
			return entries[entryptr];
		}
		// entry is not present and to be inserted before an existing element
		else
		{
			assert ( entryptr != entries.size() );
			assert ( entries[entryptr]->refpos > refpos );
			uint64_t const prevrefpos = entries[entryptr]->refpos;
			
			uint64_t const elementstomove = entries.size() - entryptr + 1;
			entries.push_back(reinterpret_cast<BamHeapEntry *>(0));
			for ( uint64_t i = 0; i < elementstomove; ++i )
				entries[entries.size()-i-1] = entries[entries.size()-i-2];
				
			assert ( entries[entryptr+1]->refpos == prevrefpos );
			
			entries[entryptr] = getEntry();
			entries[entryptr]->refpos = refpos;
			
			return entries[entryptr];
		}
	}
	
	void handleFinished(
		libmaus::bambam::BamHeader const & header,
		uint64_t refid, 
		uint64_t refpos
	)
	{
		#define INSPRINT
		#define REGPRINT
	
		while ( entries.size() && entries.front()->refpos < refpos )
		{
			BamHeapEntry * entry = entries.front();
			entries.pop_front();
			
			std::sort(entry->baseheap.begin(),entry->baseheap.end());
			
			if ( entry->insbases.size() )
			{
				std::vector<uint64_t> active(entry->insbases.size());
				uint64_t maxlen = 0;
				for ( uint64_t i = 0; i < active.size(); ++i )
				{
					active[i] = i;
					maxlen = std::max(maxlen,strlen(entry->insbases[i]));
				}
				
				for ( uint64_t j = 0; j < maxlen; ++j )
				{
					#if defined(INSPRINT)
					std::cerr << header.getRefIDName(refid) << "," << entry->refpos << ",-" << maxlen-j << " ";
					#endif
					
					uint64_t o = 0;
					for ( uint64_t i = 0; i < active.size(); ++i )
					{
						uint64_t const idx = active[i];
						if ( entry->insbases[idx][j].first )
						{
							#if defined(INSPRINT)
							std::cerr << entry->insbases[idx][j].first;
							#endif
							// std::cerr << "(" << static_cast<int>(entry->insqual[idx][j]) << ")";
							active[o++] = idx;
						}
					}
					active.resize(o);
					
					// std::cerr << "o=" << o << std::endl;
					
					#if defined(INSPRINT)
					for ( uint64_t k = o; k < entry->baseheap.size(); ++k )
						std::cerr.put('-');
					
					std::cerr.put('\n');
					#endif
				}
				
				#if 0
				for ( uint64_t i = 0; i < entry->insbases.size(); ++i )
					std::cerr << "(" << entry->insbases[i] << ")";
				std::cerr << std::endl;
				#endif
			}
		
			#if defined(REGPRINT)	
			// if ( dif /* && del */ )
			{
				double T[6] = {0,0,0,0,0,0};
				unsigned int N[6] = {0,0,0,0,0,0};
				
				std::cerr << header.getRefIDName(refid) << "," << entry->refpos << ",0 ";
				for ( uint64_t i = 0; i < entry->baseheap.size(); ++i )
				{
					double const p = libmaus::fastx::Phred::probCorrect(entry->baseheap[i].second);
					double const q = (1.0-p)/5.0;
				
					switch ( entry->baseheap[i].first )
					{
						case bam_heap_entry_base_A:
							T[0] += p;
							T[1] += q;
							T[2] += q;
							T[3] += q;
							T[4] += q;
							T[5] += q;
							N[0] ++;
							break;
						case bam_heap_entry_base_C:
							T[1] += p;
							T[0] += q;
							T[2] += q;
							T[3] += q;
							T[4] += q;
							T[5] += q;
							N[1] ++;
							break;
						case bam_heap_entry_base_G:
							T[2] += p;
							T[0] += q;
							T[1] += q;
							T[3] += q;
							T[4] += q;
							T[5] += q;
							N[2] ++;
							break;
						case bam_heap_entry_base_T:
							T[3] += p;
							T[0] += q;
							T[1] += q;
							T[2] += q;
							T[4] += q;
							T[5] += q;
							N[3] ++;
							break;
						case bam_heap_entry_base_N:
							T[3] += p;
							T[0] += q;
							T[1] += q;
							T[2] += q;
							T[4] += q;
							T[5] += q;
							N[4] ++;
							break;
						case bam_heap_entry_base_DEL:
							T[5] += p;
							T[0] += q;
							T[1] += q;
							T[2] += q;
							T[3] += q;
							T[4] += q;
							N[5] ++;
							break;
					}
					
					std::cerr << entry->baseheap[i].first;
					std::cerr << "(" << p << ")";
				}
				
				double const maxT = std::max(std::max(std::max(T[0],T[1]),std::max(T[2],T[3])),std::max(T[4],T[5]));
				int maxind[6] = {-1,-1,-1,-1,-1,-1};
				unsigned int nummax = 0;
				for ( unsigned int i = 0; i < sizeof(T)/sizeof(T[0]); ++i )
					if ( T[i] == maxT )
					{
						maxind[nummax++] = i;
					}
				assert ( nummax );
				
				std::cerr << " ";
				for ( unsigned int i = 0; i < nummax; ++i )
				{
					double const p = T[maxind[i]]/entry->baseheap.size();
					unsigned int const n = N[maxind[i]];

					std::cerr << static_cast<bam_heap_entry_base>(maxind[i]);
					std::cerr << "(" << p << "/" << n << ")";
				}

				std::cerr << "\n";
			}
			#endif
						
			entryfreelist.push_back(entry);
		}
	}
	
	void processPendingInsert(uint64_t const refpos, std::vector< std::pair<uint8_t,uint8_t> * > & pendinginserts)
	{
		if ( pendinginserts.size() )
		{
			BamHeapEntry * entry = getEntry(refpos);

			if ( pendinginserts.size() == 1 )
			{
				entry->insbases.push_back(pendinginserts.front());
			}
			else
			{
				uint64_t acclen = 0;
				for ( uint64_t i = 0; i < pendinginserts.size(); ++i )
					acclen += strlen(pendinginserts[i]);
				
				std::pair<uint8_t,uint8_t> * conc = getPairString(acclen);

				acclen = 0;
				for ( uint64_t i = 0; i < pendinginserts.size(); ++i )
				{
					uint64_t const llen = strlen(pendinginserts[i]);
					std::copy(pendinginserts[i], pendinginserts[i]+llen,conc+acclen);
					pairstringfreelist[llen+1].push_back(pendinginserts[i]);
					acclen += llen;
				}
				assert ( 
					(conc[acclen] == std::pair<uint8_t,uint8_t>(0,0))
				);

				entry->insbases.push_back(conc);
			}
			
			pendinginserts.resize(0);
		}	
	}
};

int bamheap(libmaus::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamHeader const & header = dec.getHeader();	
	::libmaus::bambam::BamAlignment const & algn = dec.getAlignment();

	std::vector< ::libmaus::bambam::BamAlignment::shared_ptr_type > algnpool;
	std::vector< ::libmaus::bambam::BamAlignment * > algnfreelist;
	
	
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	libmaus::autoarray::AutoArray<char> bases;
	
	EntriesContainer entcnt;

	uint64_t c = 0;	
	int64_t prevrefid = -1;
	std::vector< std::pair<uint8_t,uint8_t> * > pendinginserts;
	unsigned int const delqual = 40;
	
	while ( dec.readAlignment() )
	{
		if ( algn.isMapped() && (!algn.isQCFail()) )
		{
			assert ( pendinginserts.size() == 0 );
		
			uint32_t const numcigop = algn.getCigarOperations(cigop);
			uint64_t readpos = 0;
			uint64_t refpos = algn.getPos();
			uint64_t const seqlen = algn.decodeRead(bases);
			uint8_t const * qual = libmaus::bambam::BamAlignmentDecoderBase::getQual(algn.D.begin());
			
			if ( algn.getRefID() != prevrefid )
			{
				prevrefid = algn.getRefID();
				entcnt.handleFinished(header,prevrefid,std::numeric_limits<uint64_t>::max());
			}
			else
			{
				entcnt.handleFinished(header,prevrefid,refpos);
			}
			
			entcnt.resetPointer();
									
			for ( uint64_t ci = 0; ci < numcigop; ++ci )
			{
				uint64_t const ciglen = cigop[ci].second;
				
				switch ( cigop[ci].first )
				{
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
					{
						entcnt.processPendingInsert(refpos,pendinginserts);
						
						for ( uint64_t i = 0; i < ciglen; ++i )
						{
							BamHeapEntry * entry = entcnt.getEntry(refpos);
							
							switch ( bases[readpos] )
							{
								case 'A':
									entry->baseheap.push_back(std::pair<bam_heap_entry_base,uint8_t>(bam_heap_entry_base_A,qual[readpos]));
									break;
								case 'C':
									entry->baseheap.push_back(std::pair<bam_heap_entry_base,uint8_t>(bam_heap_entry_base_C,qual[readpos]));
									break;
								case 'G':
									entry->baseheap.push_back(std::pair<bam_heap_entry_base,uint8_t>(bam_heap_entry_base_G,qual[readpos]));
									break;
								case 'T':
									entry->baseheap.push_back(std::pair<bam_heap_entry_base,uint8_t>(bam_heap_entry_base_T,qual[readpos]));
									break;
								default:
									entry->baseheap.push_back(std::pair<bam_heap_entry_base,uint8_t>(bam_heap_entry_base_N,qual[readpos]));
									break;
							}
							
							readpos++;
							refpos++;
						}
						break;
					}
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
					{
						std::pair<uint8_t,uint8_t> * p = entcnt.getPairString(ciglen);
						for ( uint64_t i = 0; i < ciglen; ++i )
							p[i] = std::pair<uint8_t,uint8_t>(bases[readpos+i],qual[readpos+i]);
						p[ciglen] = std::pair<uint8_t,uint8_t>(0,0);
						
						pendinginserts.push_back(p);
												
						readpos += ciglen;
						break;
					}
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL:
						// deleting bases from the reference
						for ( uint64_t i = 0; i < ciglen; ++i )
						{
							BamHeapEntry * entry = entcnt.getEntry(refpos);
							entry->baseheap.push_back(
								std::pair<bam_heap_entry_base,uint8_t>
								(
									bam_heap_entry_base_DEL,delqual)
								);
						
							refpos++;
						}
						break;
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP:
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
						std::pair<uint8_t,uint8_t> * p = entcnt.getPairString(ciglen);
						std::fill(p,p+ciglen,std::pair<uint8_t,uint8_t>('*',0));
						p[ciglen] = std::pair<uint8_t,uint8_t>(0,0);
						pendinginserts.push_back(p);
						break;
					}
				}
			}

			entcnt.processPendingInsert(refpos,pendinginserts);
			
			assert ( readpos == seqlen );
		}
		
		c += 1;
		if ( verbose && (c & (1024*1024-1)) == 0 )
			std::cerr << "[V] " << c/(1024*1024) << std::endl;			
	}

	entcnt.handleFinished(header,prevrefid,std::numeric_limits<uint64_t>::max());

	#if 0
	assert ( entcnt.pairstringpool.size() == entcnt.pairstringfreelist.size() );
	for ( uint64_t i = 0; i < entcnt.pairstringpool.size(); ++i )
		assert (
			entcnt.pairstringpool[i].size()
			==
			entcnt.pairstringfreelist[i].size()
		);
	#endif
		
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
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamheap(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;	
	}
}
