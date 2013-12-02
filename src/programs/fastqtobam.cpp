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
#include <libmaus/fastx/StreamFastQReader.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/util/MemUsage.hpp>

#include <libmaus/lz/BufferedGzipStream.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <iomanip>

#include <config.h>

int getDefaultMD5()
{
	return 0;
}

static int getDefaultLevel() 
{
	return Z_DEFAULT_COMPRESSION;
}
static int getDefaultVerbose() 
{
	return 0;
}
static int getDefaultGz() 
{
	return 0;
}

static int getLevel(libmaus::util::ArgInfo const & arginfo)
{
	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	
	switch ( level )
	{
		case Z_NO_COMPRESSION:
		case Z_BEST_SPEED:
		case Z_BEST_COMPRESSION:
		case Z_DEFAULT_COMPRESSION:
			break;
		default:
		{
			::libmaus::exception::LibMausException se;
			se.getStream()
				<< "Unknown compression level, please use"
				<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
				<< " level=" << Z_BEST_SPEED << " (fast) or"
				<< " level=" << Z_BEST_COMPRESSION << " (best) or"
				<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
			se.finish();
			throw se;
		}
			break;
	}
	
	return level;
}

struct NameInfo
{
	bool ispair;
	bool isfirst;
	uint64_t gl;
	uint64_t gr;
	
	NameInfo() : ispair(false), isfirst(false), gl(0), gr(0) {}
	NameInfo(std::string const & name, libmaus::fastx::SpaceTable const & ST)
	: ispair(false), isfirst(true), gl(0), gr(name.size())
	{
		uint64_t l = 0;
		// skip space at start
		while ( l != name.size() && ST.spacetable[name[l]] )
			++l;
		while ( l != name.size() )
		{
			// skip non space symbols
			uint64_t r = l;
			while ( (r != name.size()) && (!ST.spacetable[name[r]]) )
				++r;
			
			// check whether this part of the name ends on /1 or /2
			if ( r-l >= 2 && name[r-2] == '/' && (name[r-1] == '1' || name[r-1] == '2') )
			{
				if ( name[r-1] == '2' )
					isfirst = false;
					
				gl = l;
				gr = r;
				ispair = true;
			}
			
			l = r;
			// skip spaces
			while ( l != name.size() && ST.spacetable[name[l]] )
				++l;
		}	
	}
};

template<typename writer_type>
void fastqtobamSingle(std::istream & in, writer_type & bamwr, int const verbose = 0)
{
	typedef ::libmaus::fastx::StreamFastQReaderWrapper reader_type;
	typedef reader_type::pattern_type pattern_type;
	::libmaus::fastx::StreamFastQReaderWrapper fqin(std::cin);	
	pattern_type element;
	libmaus::fastx::SpaceTable ST;
	uint64_t proc = 0;
	
	while ( fqin.getNextPatternUnlocked(element) )
	{
		std::string const & name = element.sid;
		NameInfo const NI(name,ST);

		if ( NI.ispair )
		{
			bamwr.encodeAlignment(
				name.substr(NI.gl,NI.gr-NI.gl-2),
				-1,
				-1,
				0,
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
				(NI.isfirst ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
				std::string(),
				-1,
				-1,
				0,
				element.spattern,
				element.quality
			);
			bamwr.commit();
		}
		else
		{
			bamwr.encodeAlignment(
				name,
				-1,
				-1,
				0,
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP,
				std::string(),
				-1,
				-1,
				0,
				element.spattern,
				element.quality
			);
			bamwr.commit();
		
		}

		proc += 1;
		if ( verbose && ((proc & (1024*1024-1)) == 0) )
			std::cerr << "[V] " << proc << std::endl;
	}
	if ( verbose )
		std::cerr << "[V] " << proc << std::endl;
}

template<typename writer_type>
void fastqtobamPair(std::istream & in_1, std::istream & in_2, writer_type & bamwr, int const verbose = 0)
{
	typedef ::libmaus::fastx::StreamFastQReaderWrapper reader_type;
	typedef reader_type::pattern_type pattern_type;
	::libmaus::fastx::StreamFastQReaderWrapper fqin_1(in_1);	
	::libmaus::fastx::StreamFastQReaderWrapper fqin_2(in_2);
	pattern_type element_1, element_2;
	libmaus::fastx::SpaceTable ST;
	uint64_t proc = 0;
	
	while ( fqin_1.getNextPatternUnlocked(element_1) )
	{
		bool const ok_2 = fqin_2.getNextPatternUnlocked(element_2);
		
		if ( ! ok_2 )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "fastq input files contain a different number of reads" << std::endl;
			ex.finish();
			throw ex;
		}

		std::string const & name_1 = element_1.sid;
		std::string const & name_2 = element_2.sid;
		NameInfo const NI_1(name_1,ST);
		NameInfo const NI_2(name_2,ST);
		
		if ( (!NI_1.ispair) || (!NI_1.isfirst) )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "name " << name_1 << " does not look like a first mate read name" << std::endl;
			ex.finish();
			throw ex;		
		}
		if ( NI_2.isfirst )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "name " << name_2 << " does not look like a second mate read name" << std::endl;
			ex.finish();
			throw ex;		
		}

		bamwr.encodeAlignment(
			name_1.substr(NI_1.gl,NI_1.gr-NI_1.gl-2),
			-1,
			-1,
			0,
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
			(NI_1.isfirst ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
			std::string(),
			-1,
			-1,
			0,
			element_1.spattern,
			element_1.quality
		);
		bamwr.commit();

		bamwr.encodeAlignment(
			name_2.substr(NI_2.gl,NI_2.gr-NI_2.gl-2),
			-1,
			-1,
			0,
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
			(NI_2.isfirst ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
			std::string(),
			-1,
			-1,
			0,
			element_2.spattern,
			element_2.quality
		);
		bamwr.commit();
		
		proc += 2;
		if ( verbose && ((proc & (1024*1024-1)) == 0) )
			std::cerr << "[V] " << proc << std::endl;
	}

	if ( verbose )
		std::cerr << "[V] " << proc << std::endl;
}


void fastqtobam(libmaus::util::ArgInfo const & arginfo)
{
	std::vector<std::string> filenames = arginfo.getPairValues("I");
	for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
		filenames.push_back(arginfo.restargs[i]);

	int const level = getLevel(arginfo);
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const gz = arginfo.getValue<int>("gz",getDefaultGz());

	::libmaus::bambam::BamHeader bamheader;
	std::ostringstream headerostr;
	headerostr << "@HD\tVN:1.4\tSO:unknown\n";
	headerostr 
		<< "@PG"<< "\t" 
		<< "ID:" << "fastqtobam" << "\t" 
		<< "PN:" << "fastqtobam" << "\t"
		<< "CL:" << arginfo.commandline << "\t"
		<< "VN:" << std::string(PACKAGE_VERSION)
		<< std::endl;

	bamheader.text = headerostr.str();		

	/*
	 * start md5 callbacks
	 */
	std::string md5filename;

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
			md5filename = arginfo.getUnparsedValue("md5filename","");
		else
			std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

		if ( md5filename.size() )
		{
			::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
			Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
			cbs.push_back(Pmd5cb.get());
		}
	}
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5 callbacks
	 */

	::libmaus::bambam::BamWriter::unique_ptr_type bamwr(
		new ::libmaus::bambam::BamWriter(
			std::cout,bamheader,level,Pcbs
		)
	);

	if ( filenames.size() == 0 )
	{
		if ( gz )
		{
			libmaus::lz::BufferedGzipStream BGS(std::cin);
			fastqtobamSingle(BGS,*bamwr);
		}
		else
		{
			fastqtobamSingle(std::cin,*bamwr);		
		}
	}
	else if ( filenames.size() == 1 )
	{
		libmaus::aio::CheckedInputStream CIS(filenames[0]);		
		
		if ( gz )
		{
			libmaus::lz::BufferedGzipStream BGS(CIS);
			fastqtobamSingle(BGS,*bamwr);			
		}
		else
		{
			fastqtobamSingle(CIS,*bamwr);		
		}
	}
	else
	{
		if ( filenames.size() > 2 )
			std::cerr << "[D] warning, ignoring additional input files past first two" << std::endl;
	
		libmaus::aio::CheckedInputStream CIS_1(filenames[0]);
		libmaus::aio::CheckedInputStream CIS_2(filenames[1]);

		if ( gz )
		{
			libmaus::lz::BufferedGzipStream BGS_1(CIS_1);
			libmaus::lz::BufferedGzipStream BGS_2(CIS_2);
			fastqtobamPair(BGS_1,BGS_2,*bamwr);			
		}
		else
		{
			fastqtobamPair(CIS_1,CIS_2,*bamwr);		
		}
	}

	bamwr.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::timing::RealTimeClock rtc; rtc.start();
		
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

				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", "compression setting for output BAM file (default: -1, zlib default settings)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "gz=<[0]>", "input is gzip compressed FastQ (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<[0]>", "print progress report (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[input file name]>", "input file names (standard input if unset)" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				
				std::cerr << "The I key can be given twice for a pair of synced FastQ files." << std::endl;
				std::cerr << "Any none key=value arguments will be considered as input file names." << std::endl;

				return EXIT_SUCCESS;
			}
		
			
		fastqtobam(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
