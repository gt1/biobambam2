/**
    biobambam
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

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <iomanip>

#include <config.h>

#include <libmaus/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus/bambam/BamToFastqOutputFileSet.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/util/Histogram.hpp>
#include <libmaus/bambam/CollatingBamDecoderAlignmentInputCallback.hpp>

struct DepthHist : public libmaus::bambam::CollatingBamDecoderAlignmentInputCallback
{
	static uint64_t const skipthres = 1024;

	libmaus::util::Histogram hist;
	std::pair<int64_t,int64_t> prev;
	int64_t actchr;
	
	std::deque < uint64_t > Q;
	std::vector < uint32_t > D;
	
	libmaus::bambam::BamHeader const * bamheader;
	
	DepthHist()
	: prev(-1,-1), actchr(-1), bamheader(0)
	{
	
	}
	
	void handleD()
	{
		std::cerr << "[V] handling D array." << std::endl;
		std::string name;
		
		if ( bamheader && actchr >= 0 )
		{
			std::cerr << "[V] getting chromosome data for actchr=" 
				<< actchr 
				<< " header is " << bamheader
				<< " number of chromosomes is " << bamheader->chromosomes.size()
				<< std::endl;
			
			uint64_t const exlen = bamheader->chromosomes.at(actchr).len;
			
			std::cerr << "got exlen=" << exlen << std::endl;
			
			name = bamheader->chromosomes.at(actchr).name;
			
			std::cerr << "[V] got length " << exlen << " name " << name << std::endl;
			
			while ( D.size() < exlen )
				D.push_back(0);
				
			std::cerr << "[V] added missing values." << std::endl;
		}
		
		std::string const fn = "depth_chr_" + libmaus::util::NumberSerialisation::formatNumber(actchr,4) + "_" + name + ".gpl";
		libmaus::aio::CheckedOutputStream COS(fn);
		uint64_t const skip = 1000;

		std::cerr << "[V] opened file " << fn << " for writing." << std::endl;

		int64_t const ws = 500;
		for ( int64_t i = 0; i < static_cast<int64_t>(D.size()); i += skip )
		{
			int64_t const low = std::max(static_cast<int64_t>(0),i-ws);
			int64_t const high = std::min(i+ws,static_cast<int64_t>(D.size()-1));
			int64_t const diam = high-low+1;
			
			int64_t s = 0;
			for ( int64_t j = low; j < high; ++j )
				s += D[j];
			double const ds = static_cast<double>(s)/(diam);

			if ( s )
				COS << i << "\t" << ds << std::endl;
		}
		
		COS.flush();

		std::cerr << "Handle actchr=" << actchr << " size " << D.size() << std::endl;
	}

	void flush()
	{
		while ( Q.size() )
		{
			hist ( Q.front() );
			D.push_back(Q.front());
			Q.pop_front();
		}

		if ( actchr != -1 )
			handleD();
	}

	void operator()(::libmaus::bambam::BamAlignment const & A)
	{
		// consider only mapped pairs
		if (  (!A.isUnmap()) && (!A.isMateUnmap()) )
		{
			// get coordinates
			int64_t const chr = A.getRefID();
			int64_t const pos = A.getPos();
			
			std::pair<int64_t,int64_t> const next(chr,pos);
			
			// check for order
			if ( next < prev )
			{
				libmaus::exception::LibMausException se;
				se.getStream() << "File is not sorted: "
					<< " prev=(" << prev.first << "," << prev.second << ")"
					<< " next=(" << chr << "," << pos << ")"
					<< std::endl;
				se.finish();
				throw se;
			}
			
			// new sequence
			if ( chr != prev.first )
			{
				// handle entries which are still in the queue
				while ( Q.size() )
				{
					assert ( Q.front() );
					hist(Q.front());
					D.push_back(Q.front());
					Q.pop_front();
				}
			
				if ( D.size() )
					handleD();	
				assert ( Q.size() == 0 );
				
				// handle D vector
				
				D.resize(0);
				
				// depth 0 before first mapping
				for ( uint64_t i = 0; i < static_cast<uint64_t>(pos); ++i )
					D.push_back(0);

				actchr = chr;
			}
			else
			{
				// pop old values now out of range
				while ( Q.size() && prev.second != pos )
				{
					hist(Q.front());
					D.push_back(Q.front());
					Q.pop_front();
					prev.second++;
				}
				
				bool const addSkipToHist = (pos-prev.second) <= static_cast<int64_t>(skipthres);
				
				while ( prev.second != pos )
				{
					if ( addSkipToHist )
						hist(0);
					D.push_back(0);
					prev.second++;
				}
			}
			
			uint64_t const reflen = A.getReferenceLength();

			for ( uint64_t i = 0; i < reflen; ++i )
				if ( i < Q.size() )
					Q[i]++;
				else
					Q.push_back(1);
			
			prev.first = chr;
			prev.second = pos;
		}
	}
};

struct BamDistHistInputFileStream
{
	std::string const fn;
	libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	std::istream & in;
	
	static libmaus::aio::CheckedInputStream::unique_ptr_type openFile(std::string const & fn)
	{
		return UNIQUE_PTR_MOVE(libmaus::aio::CheckedInputStream::unique_ptr_type(new libmaus::aio::CheckedInputStream(fn)));
	}
	
	BamDistHistInputFileStream(libmaus::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  CIS(
		(fn != "-") ? UNIQUE_PTR_MOVE(openFile(fn)) : UNIQUE_PTR_MOVE(libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}

	BamDistHistInputFileStream(std::string const & rfn)
	: fn(rfn), CIS(
		(fn != "-") ? UNIQUE_PTR_MOVE(openFile(fn)) : UNIQUE_PTR_MOVE(libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}
};

void bamdisthist(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	libmaus::bambam::BamToFastqOutputFileSet OFS(arginfo);

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	libmaus::util::Histogram hist;
	libmaus::util::Histogram tlenhist;
	libmaus::util::Histogram tlenproperhist;
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( 
			ob->fpair &&
			(!libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Da))) &&
			(!libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Db)))
		)
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getAuxRank(ob->Da, ob->blocksizea);
			uint64_t const rankb = libmaus::bambam::BamAlignmentDecoderBase::getAuxRank(ob->Db, ob->blocksizeb);
			
			if ( ranka > rankb )
				hist(ranka-rankb);
			else
				hist(rankb-ranka);
				
			int32_t const tlen = std::abs(libmaus::bambam::BamAlignmentDecoderBase::getTlen(ob->Da));
			
			if ( libmaus::bambam::BamAlignmentDecoderBase::isProper(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Da)))
				tlenproperhist(tlen);
			tlenhist(tlen);

			cnt += 2;
		}
		
		if ( precnt >> verbshift != cnt >> verbshift )
		{
			std::cerr 
				<< "[V] " << (cnt >> 20) 
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() 
				<< std::endl;
		}
	}
	std::cerr 
		<< "[V] " << (cnt >> 20) 
		<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() 
		<< std::endl;	
	
	libmaus::aio::CheckedOutputStream disthiststr("disthist.gpl");
	hist.print(disthiststr);
	disthiststr.flush();
	disthiststr.close();

	libmaus::aio::CheckedOutputStream tlenhiststr("tlenhist.gpl");
	tlenhist.print(tlenhiststr);
	tlenhiststr.flush();
	tlenhiststr.close();

	libmaus::aio::CheckedOutputStream tlenhistproperstr("tlenhistproper.gpl");
	tlenproperhist.print(tlenhistproperstr);
	tlenhistproperstr.flush();
	tlenhistproperstr.close();
}

void bamdisthist(libmaus::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY"));
	libmaus::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat","bam");
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	DepthHist dh;

	if ( inputformat == "bam" )
	{
		BamDistHistInputFileStream bamin(inputfilename);
		libmaus::bambam::BamCircularHashCollatingBamDecoder CHCBD(bamin.in,
			tmpfilename,excludeflags,
			true /* put rank */
		);
		dh.bamheader = &(CHCBD.getHeader());
		CHCBD.setInputCallback(&dh);
		bamdisthist(arginfo,CHCBD);
		dh.flush();
	}
	#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"r","",
			tmpfilename,excludeflags,
			true /* put rank */
		);
		dh.bamheader = &(CHCBD.getHeader());
		CHCBD.setInputCallback(&dh);
		bamdisthist(arginfo,CHCBD);
		dh.flush();
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			true /* put rank */
		);
		dh.bamheader = &(CHCBD.getHeader());
		CHCBD.setInputCallback(&dh);
		bamdisthist(arginfo,CHCBD);
		dh.flush();
	}
	#endif
	else
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "unknown input format " << inputformat << std::endl;
		se.finish();
		throw se;
	}
	
	libmaus::aio::CheckedOutputStream depthhiststr("depthhist.gpl");
	dh.hist.print(depthhiststr);
	depthhiststr.flush();
	depthhiststr.close();
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
				std::cerr << ::biobambam::Licensing::license() << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
				
				V.push_back ( std::pair<std::string,std::string> ( "filename=<[stdin]>", "input filename (default: read file from standard input)" ) );
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format, cram, bam or sam" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format, bam" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "exclude=<[SECONDARY]>", "exclude alignments matching any of the given flags" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("T=<[") + arginfo.getDefaultTmpFileName() + "]>" , "temporary file name" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				std::cerr << "Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP" << std::endl;

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		bamdisthist(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
