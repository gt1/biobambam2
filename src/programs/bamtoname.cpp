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

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

#include <iomanip>

#include <config.h>

#include <libmaus2/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus2/bambam/BamToFastqOutputFileSet.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

struct BamToFastQInputFileStream
{
	std::string const fn;
	libmaus2::aio::InputStreamInstance::unique_ptr_type CIS;
	std::istream & in;
	
	static libmaus2::aio::InputStreamInstance::unique_ptr_type openFile(std::string const & fn)
	{
		libmaus2::aio::InputStreamInstance::unique_ptr_type ptr(new libmaus2::aio::InputStreamInstance(fn));
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	BamToFastQInputFileStream(libmaus2::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  CIS(
		(fn != "-") ? (openFile(fn)) : (libmaus2::aio::InputStreamInstance::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}

	BamToFastQInputFileStream(std::string const & rfn)
	: fn(rfn), CIS(
		(fn != "-") ? (openFile(fn)) : (libmaus2::aio::InputStreamInstance::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}
};

void bamtonameNonCollating(libmaus2::util::ArgInfo const & arginfo, libmaus2::bambam::BamAlignmentDecoder & bamdec)
{
	libmaus2::timing::RealTimeClock rtc; rtc.start();
	uint32_t const excludeflags = libmaus2::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
	::libmaus2::autoarray::AutoArray<uint8_t> T;
	uint64_t cnt = 0;
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
		
	while ( bamdec.readAlignment() )
	{
		uint64_t const precnt = cnt++;
		
		if ( ! (algn.getFlags() & excludeflags) )
		{
			char const * name = libmaus2::bambam::BamAlignmentDecoderBase::getReadName(algn.D.begin());
			uint64_t const la = libmaus2::bambam::BamAlignmentDecoderBase::getLReadName(algn.D.begin())-1;
			std::cout.write(name,la);
			std::cout.put('\n');
			bcnt += la+1;
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

void bamtonameNonCollating(libmaus2::util::ArgInfo const & arginfo)
{
	std::string const inputformat = arginfo.getValue<std::string>("inputformat","bam");
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	
	if ( inputformat == "bam" )
	{
		BamToFastQInputFileStream bamin(inputfilename);
		libmaus2::bambam::BamDecoder bamdec(bamin.in);
		bamtonameNonCollating(arginfo,bamdec);
	}
	#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus2::bambam::ScramDecoder bamdec(inputfilename,"r","");
		bamtonameNonCollating(arginfo,bamdec);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus2::bambam::ScramDecoder bamdec(inputfilename,"r",reference);
		bamtonameNonCollating(arginfo,bamdec);
	}
	#endif
	else
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "unknown input format " << inputformat << std::endl;
		se.finish();
		throw se;
	}
		
	std::cout.flush();
}

void bamtonameCollating(
	libmaus2::util::ArgInfo const & arginfo,
	libmaus2::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	libmaus2::bambam::BamToFastqOutputFileSet OFS(arginfo);

	libmaus2::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus2::timing::RealTimeClock rtc; rtc.start();
	::libmaus2::autoarray::AutoArray<uint8_t> T;
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair )
		{
			char const * namea = libmaus2::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			uint64_t const la = libmaus2::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1;
			OFS.Fout.write(namea,la);
			OFS.Fout.put('\n');

			char const * nameb = libmaus2::bambam::BamAlignmentDecoderBase::getReadName(ob->Db);
			uint64_t const lb = libmaus2::bambam::BamAlignmentDecoderBase::getLReadName(ob->Db)-1;
			OFS.F2out.write(nameb,lb);
			OFS.F2out.put('\n');

			cnt += 2;
			bcnt += (la+lb);
		}
		else if ( ob->fsingle )
		{
			char const * namea = libmaus2::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			uint64_t const la = libmaus2::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1;
			OFS.Sout.write(namea,la);
			OFS.Sout.put('\n');

			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan1 )
		{
			char const * namea = libmaus2::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			uint64_t const la = libmaus2::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1;
			OFS.Oout.write(namea,la);
			OFS.Oout.put('\n');

			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan2 )
		{
			char const * namea = libmaus2::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			uint64_t const la = libmaus2::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1;
			OFS.O2out.write(namea,la);
			OFS.O2out.put('\n');

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
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() 
				<< "\t" << cnt
				<< "\t" << CHCBD.getRank()
				<< std::endl;
		}
	}
	
	std::cerr << "[V] " << cnt << "\t" << CHCBD.getRank() << std::endl;
}

void bamtonameCollating(libmaus2::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus2::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus2::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat","bam");
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");

	if ( inputformat == "bam" )
	{
		BamToFastQInputFileStream bamin(inputfilename);
		libmaus2::bambam::BamCircularHashCollatingBamDecoder CHCBD(bamin.in,
			tmpfilename,excludeflags,
			true /* put rank */
			);
		bamtonameCollating(arginfo,CHCBD);
	}
	#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus2::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"r","",
			tmpfilename,excludeflags,
			true /* put rank */
		);
		bamtonameCollating(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus2::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			true /* put rank */
		);
		bamtonameCollating(arginfo,CHCBD);
	}
	#endif
	else
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "unknown input format " << inputformat << std::endl;
		se.finish();
		throw se;
	}
	
	std::cout.flush();
}

void bamtoname(libmaus2::util::ArgInfo const & arginfo)
{
	if ( arginfo.getValue<uint64_t>("collate",1) )
		bamtonameCollating(arginfo);
	else
		bamtonameNonCollating(arginfo);
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo const arginfo(argc,argv);
		
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
				std::cerr << ::biobambam2::Licensing::license() << std::endl;
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
				#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format, cram, bam or sam" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format, bam" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "exclude=<[SECONDARY,SUPPLEMENTARY]>", "exclude alignments matching any of the given flags" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("T=<[") + arginfo.getDefaultTmpFileName() + "]>" , "temporary file name" ) );
				
				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				std::cerr << "Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" << std::endl;

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		bamtoname(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
