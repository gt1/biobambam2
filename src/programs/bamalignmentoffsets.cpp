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
#include <config.h>
#include <cstdlib>
#include <iostream>
#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/aio/PosixFdOutputStream.hpp>
#include <libmaus/fastx/FastAReader.hpp>
#include <libmaus/fastx/FastABgzfIndex.hpp>
#include <libmaus/fastx/FastAIndex.hpp>
#include <libmaus/fastx/StreamFastAReader.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/bambam/MdNmRecalculation.hpp>
#include <libmaus/lcs/SuffixPrefix.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/lz/RAZFIndex.hpp>
#include <libmaus/lz/RAZFDecoder.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GetFileSize.hpp>
#include <libmaus/util/OutputFileNameTools.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

static uint64_t getDefaultIOBlockSize() { return 128*1024; }

std::string getClippedRead(libmaus::bambam::BamAlignment const & algn)
{
	uint64_t const f = algn.getFrontSoftClipping();
	uint64_t const b = algn.getBackSoftClipping();
	
	std::string r = algn.getRead();
	r = r.substr(0,r.size()-b); 
	r = r.substr(f);
	
	return r;
}

std::vector<libmaus::bambam::BamFlagBase::bam_cigar_ops> getCigarOperations(
	libmaus::bambam::BamAlignment const & algn,
	bool removeFrontSoftClip = false,
	bool removeBackSoftClip = false
)
{
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	algn.getCigarOperations(cigop);
	std::vector<libmaus::bambam::BamFlagBase::bam_cigar_ops> V;
	for ( uint64_t i = 0; i < cigop.size(); ++i )
		for ( int64_t j = 0; j < cigop[i].second; ++j )
			V.push_back(static_cast<libmaus::bambam::BamFlagBase::bam_cigar_ops>(cigop[i].first));
			
	if ( removeFrontSoftClip )
	{
		uint64_t f = algn.getFrontSoftClipping();
		std::reverse(V.begin(),V.end());
		while ( f-- )
			V.pop_back();
		std::reverse(V.begin(),V.end());
	}
	
	if ( removeBackSoftClip )
	{
		uint64_t b = algn.getBackSoftClipping();
		while ( b-- )
			V.pop_back();	
	}
			
	return V;
}

std::string encodeCigarString(std::vector<libmaus::bambam::BamFlagBase::bam_cigar_ops> const & V)
{
	uint64_t low = 0;
	std::ostringstream ostr;
	while ( low != V.size() )
	{
		uint64_t high = low;
		while ( high != V.size() && V[high] == V[low] )
			++high;

		ostr << (high-low);
		switch ( V[low] )
		{
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
				ostr.put('M');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
				ostr.put('I');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL:
				ostr.put('D');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP:
				ostr.put('N');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP:
				ostr.put('S');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP:
				ostr.put('H');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CPAD:
				ostr.put('P');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
				ostr.put('=');
				break;
			case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
				ostr.put('X');
				break;
		}
			
		low = high;
	}
	
	return ostr.str();
}

std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> handleChain(std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> const & chain)
{
	int64_t previd = -1;
	int64_t prevend = -1;
	
	std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> outchain;

	for ( uint64_t i = 0; i < chain.size(); ++i )
	{
		libmaus::bambam::BamAlignment const & algn = *(chain[i]);

		int64_t const thisid = algn.getRefID();
		int64_t const thispos = algn.getPos();
		int64_t const thisend = algn.getAlignmentEnd();

		std::cerr << algn.getName() << "\t" << "[" << thispos << "," << thisend << "]" << "\t" << (thisend-thispos+1);

		int64_t const offset = ( thisid == previd && thisid >= 0 ) ? (thispos-prevend) : std::numeric_limits<int64_t>::min();
		if ( thisid == previd && thisid >= 0 )
		{
			std::cerr << "\t" << offset;
			
		}
		std::cerr << std::endl;

		if ( offset < 0 && offset != std::numeric_limits<int64_t>::min() )
		{	
			std::string prevread = getClippedRead(*(outchain.back()));
			uint64_t const prevreadkeep = std::min(
				static_cast<uint64_t>(prevread.size()),static_cast<uint64_t>(5*-offset)
			);
			prevread = prevread.substr(prevread.size()-prevreadkeep);

			std::string thisread = getClippedRead(algn);
			uint64_t const thisreadkeep = std::min(
				static_cast<uint64_t>(thisread.size()),static_cast<uint64_t>(5*-offset)
			);
			thisread = thisread.substr(0,thisreadkeep);
				::libmaus::lcs::SuffixPrefix SP(prevread.size(),thisread.size());
			::libmaus::lcs::SuffixPrefixResult SSPR = SP.process(prevread.begin(),thisread.begin());

			//if ( !ok && !SP.startsWithDeletion() && !SP.endsWithInsertion() && SSPR.nummat >= 50 )
			std::cerr << SSPR << std::endl;
			SP.printAlignmentLines(std::cerr, prevread, thisread, SSPR, 120);
			
			std::cerr << "back clip a? ";
			int64_t backclipa = -1;
			std::cin >> backclipa;

			if ( backclipa > 0 )
			{
				uint64_t soft = 0, hard = 0;
				
				std::vector<libmaus::bambam::BamFlagBase::bam_cigar_ops> cigops = getCigarOperations(
					*(outchain.back()),false,false);
				
				while ( cigops.size() && 
					cigops.back() == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP
				)
				{
					hard++;
					cigops.pop_back();
				}
				while ( cigops.size() && 
					cigops.back() == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP
				)
				{
					soft++;
					cigops.pop_back();
				}
				
				uint64_t tbackclipa = backclipa;
				while ( tbackclipa > 0 )
				{
					assert ( cigops.size() );
					
					switch ( cigops.back() )
					{
						case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
						case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
						case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
						case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
							tbackclipa -= 1;
							break;
						default:
							break;
					}
					
					cigops.pop_back();
				}
				
				while ( cigops.size() && 
					(
						cigops.back() == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL 
						||
						cigops.back() == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP
						||
						cigops.back() == libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CPAD
					)
				)
				{
					cigops.pop_back();
				}
				
				for ( uint64_t i = 0; i < backclipa+soft; ++i )
					cigops.push_back(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP);
				for ( uint64_t i = 0; i < hard; ++i )
					cigops.push_back(libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP);
					
				std::string const newcig = encodeCigarString(cigops);

				std::cerr << outchain.back()->getCigarString() << std::endl;
				
				outchain.back()->replaceCigarString(newcig);
			
				std::cerr << newcig << std::endl;
			}
		}

		previd = thisid;
		prevend = algn.getAlignmentEnd();
		
		outchain.push_back(algn.sclone());
	}

	return outchain;
}

void printChain(std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> const & chain)
{
	int64_t previd = -1;
	int64_t prevend = -1;
	
	for ( uint64_t i = 0; i < chain.size(); ++i )
	{
		libmaus::bambam::BamAlignment const & algn = *(chain[i]);

		int64_t const thisid = algn.getRefID();
		int64_t const thispos = algn.getPos();
		int64_t const thisend = algn.getAlignmentEnd();

		std::cerr << algn.getName() 
			<< "\t" << algn.getFrontSoftClipping()
			<< "\t" << algn.getBackSoftClipping()
			<< "\t" << "[" << thispos << "," << thisend << "]" << "\t" 
			<< (thisend-thispos+1) << "\t" << algn.getAuxAsString("ZJ");

		int64_t const offset = ( thisid == previd && thisid >= 0 ) ? (thispos-prevend) : std::numeric_limits<int64_t>::min();
		if ( thisid == previd && thisid >= 0 )
		{
			std::cerr << "\t" << offset;
			
		}
		std::cerr << std::endl;

		previd = thisid;
		prevend = algn.getAlignmentEnd();
	}
}

std::pair<bool,std::string> chainToContig(std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> const & chain)
{
	int64_t previd = -1;
	int64_t prevend = -1;
	
	bool ok = true;
	std::ostringstream ostr;
	
	for ( uint64_t i = 0; i < chain.size(); ++i )
	{
		libmaus::bambam::BamAlignment const & algn = *(chain[i]);

		int64_t const thisid = algn.getRefID();
		int64_t const thispos = algn.getPos();
		int64_t const thisend = algn.getAlignmentEnd();

		std::cerr << algn.getName() << "\t" << "[" << thispos << "," << thisend << "]" << "\t" << (thisend-thispos+1);

		int64_t const offset = ( thisid == previd && thisid >= 0 ) ? (thispos-prevend) : std::numeric_limits<int64_t>::min();
		if ( thisid == previd && thisid >= 0 )
			std::cerr << "\t" << offset;
		std::cerr << std::endl;
		
		if ( i == 0 )
			ostr << getClippedRead(algn);		
		else if ( offset >= 0 )
			ostr << std::string(offset,'N') << getClippedRead(algn);
		else
			ok = false;

		previd = thisid;
		prevend = algn.getAlignmentEnd();		
	}

	return std::pair<bool,std::string>(ok,ostr.str());
}

static int bamalignmentoffsets(libmaus::util::ArgInfo const & arginfo)
{
	uint64_t const ioblocksize = arginfo.getValueUnsignedNumeric<uint64_t>("ioblocksize",getDefaultIOBlockSize());
	int64_t const contigsplit = arginfo.getValueUnsignedNumeric<uint64_t>("contigsplit",20000);
	bool const modify = arginfo.getValue<unsigned int>("modify",0);
	std::string const fn = arginfo.restargs.at(0);
	::libmaus::aio::PosixFdInputStream PFIS(fn,ioblocksize);
	libmaus::bambam::BamDecoder bamdec(PFIS);
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	
	std::vector< std::vector<libmaus::bambam::BamAlignment::shared_ptr_type> > chains;
	
	int64_t previd = -1;
	int64_t prevend = -1;

	while ( bamdec.readAlignment() )
	{
		if ( algn.isMapped() )
		{
			int64_t const thisid = algn.getRefID();
			int64_t const thispos = algn.getPos();
			int64_t const thisend = algn.getAlignmentEnd();

			std::cerr << algn.getName() << "\t" << "[" << thispos << "," << thisend << "]" << "\t" << (thisend-thispos+1);

			int64_t const offset = ( thisid == previd && thisid >= 0 ) ? (thispos-prevend) : std::numeric_limits<int64_t>::min();
			if ( thisid == previd && thisid >= 0 )
				std::cerr << "\t" << offset;
				
			if ( offset == std::numeric_limits<int64_t>::min() || offset >= contigsplit )
				chains.push_back(std::vector<libmaus::bambam::BamAlignment::shared_ptr_type>());
				
			chains.back().push_back(algn.sclone());
			
			previd = thisid;
			prevend = algn.getAlignmentEnd();
			
			std::cerr << std::endl;
		}
	}
	
	int64_t maxchainid = -1;
	int64_t maxchainlength = -1;
	int64_t maxchainsize = 0;
	for ( uint64_t i = 0; i < chains.size(); ++i )
	{
		int64_t chainlength = 0;
		for ( uint64_t j = 0; j < chains[i].size(); ++j )
			chainlength += (chains[i][j]->getAlignmentEnd() - chains[i][j]->getPos())+1;
			
		if ( chainlength > maxchainlength )
		{
			maxchainid = i;
			maxchainlength = chainlength;
			maxchainsize = chains[i].size();
		}
	}
	
	if ( modify )
	{
		libmaus::bambam::BamWriter bw(std::cout,header);
		
		uint64_t contigid = 0;
		libmaus::aio::CheckedOutputStream faout(fn+".fa");
		
		if ( maxchainid >= 0 )
		{
			std::cerr << "[I] maximum chain length " << maxchainlength << std::endl;
			chains[maxchainid] = handleChain(chains[maxchainid]);
			for ( uint64_t i = 0; i < chains[maxchainid].size(); ++i )
				chains[maxchainid][i]->serialise(bw.getStream());

			std::pair<bool,std::string> contig = chainToContig(chains[maxchainid]);
			if ( contig.first )
			{
				faout << ">contig_" << contigid++ << " " << contig.second.size() << "\n";
				faout << contig.second << "\n";
			}
			
			for ( uint64_t j = 0; j < chains.size(); ++j )
				if ( static_cast<int64_t>(j) != maxchainid && static_cast<int64_t>(chains[j].size()) == maxchainsize )
				{
					int64_t chainlength = 0;
					for ( uint64_t i = 0; i < chains[j].size(); ++i )
						chainlength += (chains[j][i]->getAlignmentEnd() - chains[j][i]->getPos())+1;
						
					std::cerr << "[I] equivalent chain length " << chainlength << std::endl;	
					chains[j] = handleChain(chains[j]);

					for ( uint64_t i = 0; i < chains[j].size(); ++i )
						chains[j][i]->serialise(bw.getStream());

					std::pair<bool,std::string> contig = chainToContig(chains[maxchainid]);
					if ( contig.first )
					{
						faout << ">contig_" << contigid++ << " " << contig.second.size() << "\n";
						faout << contig.second << "\n";
					}
				}
		}
		
		faout.flush();
		faout.close();
		
		if ( ! contigid )
			remove((fn+".fa").c_str());
	}
	else
	{
		if ( maxchainid >= 0 )
		{
			std::cerr << "[I] maximum chain length " << maxchainlength << std::endl;
			printChain(chains[maxchainid]);

			for ( uint64_t j = 0; j < chains.size(); ++j )
				if ( static_cast<int64_t>(j) != maxchainid && static_cast<int64_t>(chains[j].size()) == maxchainsize )
				{
					int64_t chainlength = 0;
					for ( uint64_t i = 0; i < chains[j].size(); ++i )
						chainlength += (chains[j][i]->getAlignmentEnd() - chains[j][i]->getPos())+1;
						
					std::cerr << "[I] equivalent chain length " << chainlength << std::endl;	
					printChain(chains[j]);
				}

			for ( uint64_t j = 0; j < chains.size(); ++j )
				if ( static_cast<int64_t>(j) != maxchainid && static_cast<int64_t>(chains[j].size()) != maxchainsize )
				{
					int64_t chainlength = 0;
					for ( uint64_t i = 0; i < chains[j].size(); ++i )
						chainlength += (chains[j][i]->getAlignmentEnd() - chains[j][i]->getPos())+1;
						
					std::cerr << "[I] shorter chain length " << chainlength << std::endl;	
					printChain(chains[j]);
				}
		}
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
			
				V.push_back ( std::pair<std::string,std::string> ( "ioblocksize=<["+::biobambam::Licensing::formatNumber(getDefaultIOBlockSize())+"]>", "block size for I/O operations" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}


		return bamalignmentoffsets(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
