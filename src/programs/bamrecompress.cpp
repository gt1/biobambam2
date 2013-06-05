/*
    libmaus
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
*/

#include <libmaus/lz/BgzfInflateDeflateParallel.hpp>
#include <libmaus/lz/BgzfInflateDeflateParallelThread.hpp>
#include <libmaus/util/ArgInfo.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; };
static int getDefaultNumThreads() { return 1; };

uint64_t bamrecompress(libmaus::util::ArgInfo const & arginfo)
{
	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const numthreads = std::max(1,arginfo.getValue<int>("numthreads",getDefaultNumThreads()));
	
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

	libmaus::lz::BgzfInflateDeflateParallel BIDP(std::cin,std::cout,level,numthreads,4*numthreads);
	libmaus::autoarray::AutoArray<char> B(64*1024,false);
	int r;
	uint64_t t = 0;
	uint64_t last = std::numeric_limits<uint64_t>::max();
	uint64_t lcnt = 0;
	uint64_t const mod = 64*1024*1024;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	libmaus::timing::RealTimeClock lrtc; lrtc.start();

	while ( (r = BIDP.read(B.begin(),B.size())) )
	{
		BIDP.write(B.begin(),r);
		
		lcnt += r;
		t += r;
		
		if ( t/mod != last/mod )
		{
			if ( verbose )
			{
				if ( isatty(STDERR_FILENO) )
					std::cerr 
						<< "\r" << std::string(60,' ') << "\r";

				std::cerr
						<< rtc.formatTime(rtc.getElapsedSeconds()) << " " << t/(1024*1024) << "MB, " << (lcnt/lrtc.getElapsedSeconds())/(1024.0*1024.0) << "MB/s";
			
				if ( isatty(STDERR_FILENO) )
					std::cerr << std::flush;
				else
					std::cerr << std::endl;
			}
			
			lrtc.start();
			last = t;
			lcnt = 0;
		}
	}

	if ( verbose )
	{
		if ( isatty(STDERR_FILENO) )
			std::cerr 
				<< "\r" << std::string(60,' ') << "\r";
	
		std::cerr
				<< rtc.formatTime(rtc.getElapsedSeconds()) << " " << t/(1024*1024) << "MB, " << (t/rtc.getElapsedSeconds())/(1024.0*1024.0) << "MB/s";
			
		std::cerr << std::endl;
	}
	
	return 0;
}

int main(int argc, char *argv[])
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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "numthreads=<["+::biobambam::Licensing::formatNumber(getDefaultNumThreads())+"]>", "number of recoding threads" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamrecompress(arginfo);

	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
