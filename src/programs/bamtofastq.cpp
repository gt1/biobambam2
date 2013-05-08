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
#include <libmaus/util/MemUsage.hpp>

struct BamToFastQInputFileStream
{
	std::string const fn;
	libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	std::istream & in;
	
	static libmaus::aio::CheckedInputStream::unique_ptr_type openFile(std::string const & fn)
	{
		return UNIQUE_PTR_MOVE(libmaus::aio::CheckedInputStream::unique_ptr_type(new libmaus::aio::CheckedInputStream(fn)));
	}
	
	BamToFastQInputFileStream(libmaus::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  CIS(
		(fn != "-") ? UNIQUE_PTR_MOVE(openFile(fn)) : UNIQUE_PTR_MOVE(libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}

	BamToFastQInputFileStream(std::string const & rfn)
	: fn(rfn), CIS(
		(fn != "-") ? UNIQUE_PTR_MOVE(openFile(fn)) : UNIQUE_PTR_MOVE(libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}
};

void bamtofastqNonCollating(libmaus::util::ArgInfo const & arginfo, libmaus::bambam::BamAlignmentDecoder & bamdec)
{
	libmaus::timing::RealTimeClock rtc; rtc.start();
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY"));
	libmaus::bambam::BamAlignment const & algn = bamdec.getAlignment();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	uint64_t cnt = 0;
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
		
	while ( bamdec.readAlignment() )
	{
		uint64_t const precnt = cnt++;
		
		if ( ! (algn.getFlags() & excludeflags) )
		{
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQ(algn.D.begin(),T);
			std::cout.write(reinterpret_cast<char const *>(T.begin()),la);
			bcnt += la;
		}

		if ( precnt >> verbshift != cnt >> verbshift )
		{
			std::cerr 
				<< (cnt >> 20) 
				<< "\t"
				<< (static_cast<double>(bcnt)/(1024.0*1024.0))/rtc.getElapsedSeconds() << "MB/s"
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() << std::endl;
		}
	}

}

void bamtofastqNonCollating(libmaus::util::ArgInfo const & arginfo)
{
	std::string const inputformat = arginfo.getValue<std::string>("inputformat","bam");
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	
	if ( inputformat == "bam" )
	{
		BamToFastQInputFileStream bamin(inputfilename);
		libmaus::bambam::BamDecoder bamdec(bamin.in);
		bamtofastqNonCollating(arginfo,bamdec);
	}
	#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus::bambam::ScramDecoder bamdec(inputfilename,"r","");
		bamtofastqNonCollating(arginfo,bamdec);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramDecoder bamdec(inputfilename,"rc",reference);
		bamtofastqNonCollating(arginfo,bamdec);
	}
	#endif
	else
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "unknown input format " << inputformat << std::endl;
		se.finish();
		throw se;
	}
		
	std::cout.flush();
}

void bamtofastqCollating(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	libmaus::bambam::BamToFastqOutputFileSet OFS(arginfo);

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair )
		{
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
			OFS.Fout.write(reinterpret_cast<char const *>(T.begin()),la);
			uint64_t lb = libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Db,T);
			OFS.F2out.write(reinterpret_cast<char const *>(T.begin()),lb);

			cnt += 2;
			bcnt += (la+lb);
		}
		else if ( ob->fsingle )
		{
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
			OFS.Sout.write(reinterpret_cast<char const *>(T.begin()),la);

			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan1 )
		{
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
			OFS.Oout.write(reinterpret_cast<char const *>(T.begin()),la);

			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan2 )
		{
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
			OFS.O2out.write(reinterpret_cast<char const *>(T.begin()),la);

			cnt += 1;
			bcnt += (la);
		}
		
		if ( precnt >> verbshift != cnt >> verbshift )
		{
			std::cerr 
				<< "[V] "
				<< (cnt >> 20) 
				<< "\t"
				<< (static_cast<double>(bcnt)/(1024.0*1024.0))/rtc.getElapsedSeconds() << "MB/s"
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() << std::endl;
		}
	}
	
	std::cerr << "[V] " << cnt << std::endl;
}

void bamtofastqCollating(libmaus::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY"));
	libmaus::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat","bam");
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");

	if ( inputformat == "bam" )
	{
		BamToFastQInputFileStream bamin(inputfilename);
		libmaus::bambam::BamCircularHashCollatingBamDecoder CHCBD(
			bamin.in,
			tmpfilename,excludeflags,
			false /* put rank */
		);
		bamtofastqCollating(arginfo,CHCBD);
	}
	#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"r","",
			tmpfilename,excludeflags,
			false /* put rank */
		);
		bamtofastqCollating(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			false /* put rank */
		);
		bamtofastqCollating(arginfo,CHCBD);
	}
	#endif
	else
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "unknown input format " << inputformat << std::endl;
		se.finish();
		throw se;
	}
	
	std::cout.flush();
}

void bamtofastq(libmaus::util::ArgInfo const & arginfo)
{
	if ( arginfo.getValue<uint64_t>("collate",1) )
		bamtofastqCollating(arginfo);
	else
		bamtofastqNonCollating(arginfo);
		
	std::cerr << "[V] " << libmaus::util::MemUsage() << std::endl;		
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
				
				V.push_back ( std::pair<std::string,std::string> ( "F=<[stdout]>", "matched pairs first mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "F2=<[stdout]>", "matched pairs second mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "S=<[stdout]>", "single end" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "unmatched pairs first mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O2=<[stdout]>", "unmatched pairs second mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collate=<[1]>", "collate pairs" ) );
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
			
		bamtofastq(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
