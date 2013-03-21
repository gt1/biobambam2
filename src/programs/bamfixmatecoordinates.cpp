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

#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>
#include <libmaus/bambam/EncoderBase.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/lz/BufferedGzipStream.hpp>

#include <libmaus/fastx/CompactFastEncoder.hpp>
#include <libmaus/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus/util/GetObject.hpp>

#include <libmaus/lz/Deflate.hpp>
#include <libmaus/lz/SnappyCompress.hpp>
#include <libmaus/util/CountPutObject.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/timing/RealTimeClock.hpp>

// #include <biobambam/BamFile.hpp>

static std::vector < std::string > tmpfilenames;

static void removeTempFiles()
{
	for ( uint64_t i = 0; i < tmpfilenames.size(); ++i )
		remove ( tmpfilenames[i].c_str() );
}

#include <csignal>

#if defined(__APPLE__) || defined(__FreeBSD__)
typedef sig_t sighandler_t;
#endif

static sighandler_t siginthandler = 0;
static sighandler_t sigtermhandler = 0;
static sighandler_t sighuphandler = 0;
static sighandler_t sigpipehandler = 0;

static void sigIntHandler(int arg)
{
	removeTempFiles();
	if ( siginthandler )
	{
		siginthandler(arg);
	}
	else
	{
		signal(SIGINT,SIG_DFL);
		raise(SIGINT);
	}
}

static void sigHupHandler(int arg)
{
	removeTempFiles();
	if ( sighuphandler )
	{
		sighuphandler(arg);
	}
	else
	{
		signal(SIGHUP,SIG_DFL);
		raise(SIGHUP);
	}
}

static void sigTermHandler(int arg)
{
	removeTempFiles();
	if ( sigtermhandler )
	{
		sigtermhandler(arg);
	}
	else
	{
		signal(SIGTERM,SIG_DFL);
		raise(SIGTERM);
	}
}

static void sigPipeHandler(int arg)
{
	removeTempFiles();
	if ( sigpipehandler )
	{
		sigpipehandler(arg);
	}
	else
	{
		signal(SIGPIPE,SIG_DFL);
		raise(SIGPIPE);
	}
}

static void setupTempFileRemovalRoutines()
{
	siginthandler = signal(SIGINT,sigIntHandler);
	sigtermhandler = signal(SIGTERM,sigTermHandler);
	sigpipehandler = signal(SIGPIPE,sigPipeHandler);
	sighuphandler = signal(SIGHUP,sigHupHandler);
	atexit(removeTempFiles);
}

static uint64_t getMapCnt(::libmaus::bambam::CollatingBamDecoder::alignment_ptr_type const & p)
{
	if ( p && p->isMapped() )
		return 1;
	else
		return 0;
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		::libmaus::timing::RealTimeClock rtc; rtc.start();
		bool const verbose = arginfo.getValue<unsigned int>("verbose",1);
	
		setupTempFileRemovalRoutines();	
		
		unsigned int const colhashbits = arginfo.getValue<unsigned int>("colhashbits",20);
		unsigned int const collistsize = arginfo.getValue<unsigned int>("collistsize",512*1024);

		int const level = arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION);
		
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


		std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
		std::string const tmpfilename = tmpfilenamebase + "_bamcollate";
		tmpfilenames.push_back(tmpfilename);
		
		::libmaus::bambam::CollatingBamDecoder CBD(std::cin,tmpfilename,false /* put rank */,colhashbits/*hash bits*/,collistsize/*size of output list*/);
		::libmaus::bambam::BamFormatAuxiliary auxdata;
		::libmaus::bambam::BamHeader const & bamheader = CBD.bamdecoder.bamheader;
		
		// "reconstruct" command line
		std::string cl;
		for ( int i = 0; i < argc; ++i )
		{
			cl += argv[i];
			if ( i+1 < argc )
				cl += " ";
		}

		// add PG line to header
		std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
			bamheader.text,
			"bamfixmatecoordinates", // ID
			"bamfixmatecoordinates", // PN
			cl, // CL
			::libmaus::bambam::ProgramHeaderLineSet(bamheader.text).getLastIdInChain(), // PP
			std::string(PACKAGE_VERSION) // VN			
		);
		// construct new header
		::libmaus::bambam::BamHeader uphead(upheadtext);
		
		if ( uphead.getSortOrder() != "queryname" )
			uphead.changeSortOrder("unknown");
		
		// setup bam writer
		::libmaus::bambam::BamWriter writer(std::cout,uphead,level);
		
		#if 0
		::libmaus::bambam::ProgramHeaderLineSet PHLS(bamheader.text);
		std::cerr << "Last id in PG chain: " << PHLS.getLastIdInChain() << std::endl;
		#endif

		// std::cout << bamheader.text;

		typedef ::libmaus::bambam::CollatingBamDecoder::alignment_ptr_type alignment_ptr_type;
		std::pair<alignment_ptr_type,alignment_ptr_type> P;
		uint64_t const mod = 1024*1024;
		uint64_t proc = 0;
		uint64_t lastproc = 0;
		uint64_t paircnt = 0;
		
		while ( CBD.tryPair(P) )
		{
			uint64_t const mapcnt = getMapCnt(P.first) + getMapCnt(P.second);
			
			if ( mapcnt == 1 )
			{
				int32_t refid = -1;
				int32_t pos = -1;
				
				if ( P.first )
				{
					refid = P.first->getRefID();
					pos = P.first->getPos();
				}
				else
				{
					assert ( P.second );

					refid = P.second->getRefID();
					pos = P.second->getPos();
				}
				
				P.first->putRefId(refid);
				P.first->putPos(pos);
				P.first->putNextRefId(refid);
				P.first->putNextPos(pos);
				P.second->putRefId(refid);
				P.second->putPos(pos);
				P.second->putNextRefId(refid);
				P.second->putNextPos(pos);
			}
			
			if ( P.first )
			{
				P.first->serialise(writer.bgzfos);
				++proc;
			}
			if ( P.second )
			{
				P.second->serialise(writer.bgzfos);
				++proc;
			}
			if ( P.first && P.second )
			{
				paircnt++;
			}
			
			if ( verbose && (proc/mod != lastproc/mod) )
			{
				std::cerr 
					<< "Processed " << proc << " fragments, " << paircnt << " pairs, " 
					<< proc/rtc.getElapsedSeconds() << " al/s"
					<< std::endl;
				lastproc = proc;
			}
		}		

		if ( verbose )
			std::cerr 	
				<< "Processed " << proc << " fragments, " << paircnt << " pairs, " 
				<< proc/rtc.getElapsedSeconds() << " al/s"
				<< std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
