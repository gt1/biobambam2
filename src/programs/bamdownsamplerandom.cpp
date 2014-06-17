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
#include <biobambam/AttachRank.hpp>
#include <biobambam/ResetAlignment.hpp>

#include <iomanip>

#include <config.h>

#include <libmaus/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamToFastqOutputFileSet.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/random/Random.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static std::string getDefaultInputFormat() { return "bam"; }
static double getDefaultProb() { return 1.0; }

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

static int getLevel(libmaus::util::ArgInfo const & arginfo)
{
	return libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
}

struct BamDownsampleRandomInputFileStream
{
	std::string const fn;
	libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	std::istream & in;
	
	static libmaus::aio::CheckedInputStream::unique_ptr_type openFile(std::string const & fn)
	{
		libmaus::aio::CheckedInputStream::unique_ptr_type ptr(new libmaus::aio::CheckedInputStream(fn));
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	BamDownsampleRandomInputFileStream(libmaus::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  CIS(
		(fn != "-") ? openFile(fn) : (libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}

	BamDownsampleRandomInputFileStream(std::string const & rfn)
	: fn(rfn), CIS(
		(fn != "-") ? openFile(fn) : (libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}
};

template<typename decoder_type>
std::string getModifiedHeaderText(decoder_type const & bamdec, libmaus::util::ArgInfo const & arginfo, bool reset = false)
{
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	std::string const headertext(header.text);

	// add PG line to header
	std::string upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamdownsamplerandom", // ID
		"bamdownsamplerandom", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);

	if ( reset )
	{
		std::vector<libmaus::bambam::HeaderLine> allheaderlines = libmaus::bambam::HeaderLine::extractLines(upheadtext);

		std::ostringstream upheadstr;
		for ( uint64_t i = 0; i < allheaderlines.size(); ++i )
			if ( allheaderlines[i].type != "SQ" )
				upheadstr << allheaderlines[i].line << std::endl;
		upheadtext = upheadstr.str();
	}

	return upheadtext;
}

void bamdownsamplerandom(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	if ( arginfo.hasArg("seed") )
	{
		uint64_t const seed = arginfo.getValue<uint64_t>("seed",0);
		libmaus::random::Random::setup(seed);
	}
	else
	{
		libmaus::random::Random::setup();	
	}
	
	double const p = arginfo.getValue<double>("p",getDefaultProb());
	
	if ( p < 0.0 || p > 1.0 )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "Value of p must be in [0,1] but is " << p << std::endl;
		se.finish();
		throw se;
	}
	
	uint32_t const up =
		( p == 1 ) ? 
		std::numeric_limits<uint32_t>::max() :
		static_cast<uint32_t>(std::max(0.0,
			std::min(
				std::floor(p * static_cast<double>(std::numeric_limits<uint32_t>::max()) + 0.5),
				static_cast<double>(std::numeric_limits<uint32_t>::max())
			)
		))
		;

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments processed
	uint64_t cnt = 0;
	// number of bytes processed
	uint64_t bcnt = 0;
	// number of alignments written
	uint64_t ocnt = 0;
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();

	// construct new header
	::libmaus::bambam::BamHeader uphead(getModifiedHeaderText(CHCBD,arginfo));
	uphead.changeSortOrder("unknown");

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string md5filename;
	std::string indexfilename;

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
	libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	if ( arginfo.getValue<unsigned int>("index",getDefaultIndex()) )
	{
		if ( arginfo.hasArg("indexfilename") &&  arginfo.getUnparsedValue("indexfilename","") != "" )
			indexfilename = arginfo.getUnparsedValue("indexfilename","");
		else
			std::cerr << "[V] no filename for index given, not creating index" << std::endl;

		if ( indexfilename.size() )
		{
			libmaus::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Tindex(new libmaus::bambam::BgzfDeflateOutputCallbackBamIndex(tmpfileindex));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
			cbs.push_back(Pindex.get());
		}
	}
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	// construct writer
	::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(std::cout,uphead,getLevel(arginfo),Pcbs));
	typedef libmaus::bambam::BamWriter::stream_type out_stream_type;
	out_stream_type & bgzfos = writer->getStream();
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		uint32_t const rv = ::libmaus::random::Random::rand32();
		
		if ( ob->fpair )
		{
			if ( rv <= up )
			{
				::libmaus::bambam::EncoderBase::putLE<out_stream_type,uint32_t>(bgzfos,ob->blocksizea);
				bgzfos.write(reinterpret_cast<char const *>(ob->Da),ob->blocksizea);
				::libmaus::bambam::EncoderBase::putLE<out_stream_type,uint32_t>(bgzfos,ob->blocksizeb);
				bgzfos.write(reinterpret_cast<char const *>(ob->Db),ob->blocksizeb);
				ocnt += 2;
			}

			cnt += 2;
			bcnt += (ob->blocksizea+ob->blocksizeb);
		}
		else if ( ob->fsingle )
		{
			if ( rv <= up )
			{
				::libmaus::bambam::EncoderBase::putLE<out_stream_type,uint32_t>(bgzfos,ob->blocksizea);
				bgzfos.write(reinterpret_cast<char const *>(ob->Da),ob->blocksizea);
				ocnt += 1;
			}

			cnt += 1;
			bcnt += (ob->blocksizea);
		}
		else if ( ob->forphan1 )
		{
			if ( rv <= up )
			{
				::libmaus::bambam::EncoderBase::putLE<out_stream_type,uint32_t>(bgzfos,ob->blocksizea);
				bgzfos.write(reinterpret_cast<char const *>(ob->Da),ob->blocksizea);
				ocnt += 1;
			}

			cnt += 1;
			bcnt += (ob->blocksizea);
		}
		else if ( ob->forphan2 )
		{
			if ( rv <= up )
			{
				::libmaus::bambam::EncoderBase::putLE<out_stream_type,uint32_t>(bgzfos,ob->blocksizea);
				bgzfos.write(reinterpret_cast<char const *>(ob->Da),ob->blocksizea);
				ocnt += 1;
			}

			cnt += 1;
			bcnt += (ob->blocksizea);
		}
		
		if ( precnt >> verbshift != cnt >> verbshift )
		{
			std::cerr 
				<< "[V] "
				<< (cnt >> 20) 
				<< "\t"
				<< (static_cast<double>(bcnt)/(1024.0*1024.0))/rtc.getElapsedSeconds() << "MB/s"
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() 
				<< " kept " << ocnt << " (" << static_cast<double>(ocnt)/cnt << ")"
				<< std::endl;
		}
	}
	
	std::cerr << "[V] " << cnt << std::endl;

	writer.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
	}
}

void bamdownsamplerandom(libmaus::util::ArgInfo const & arginfo)
{
	if ( arginfo.hasArg("ranges") && arginfo.getValue("inputformat", getDefaultInputFormat()) != "bam" )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "ranges are only supported for inputformat=bam" << std::endl;
		se.finish();
		throw se;
	}

	if ( arginfo.hasArg("ranges") && ((!arginfo.hasArg("filename")) || arginfo.getValue<std::string>("filename","-") == "-") )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "ranges are not supported for reading via standard input" << std::endl;
		se.finish();
		throw se;
	}

	if ( arginfo.hasArg("ranges") && arginfo.getValue<uint64_t>("collate",1) > 1 )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "ranges are not supported for collate > 1" << std::endl;
		se.finish();
		throw se;
	}
	
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	uint64_t const numthreads = arginfo.getValue<uint64_t>("threads",0);

	unsigned int const hlog = arginfo.getValue<unsigned int>("colhlog",18);
	uint64_t const sbs = arginfo.getValueUnsignedNumeric<uint64_t>("colsbs",128ull*1024ull*1024ull);

	if ( inputformat == "bam" )
	{
		if ( arginfo.hasArg("ranges") )
		{
			libmaus::bambam::BamRangeCircularHashCollatingBamDecoder CHCBD(inputfilename,arginfo.getUnparsedValue("ranges",""),tmpfilename,excludeflags,false,hlog,sbs);
			bamdownsamplerandom(arginfo,CHCBD);
		}
		else
		{
			BamDownsampleRandomInputFileStream bamin(inputfilename);

			if ( numthreads > 0 )
			{
				libmaus::bambam::BamParallelCircularHashCollatingBamDecoder CHCBD(
					bamin.in,
					numthreads,
					tmpfilename,excludeflags,
					false, /* put rank */
					hlog,
					sbs
					);
				bamdownsamplerandom(arginfo,CHCBD);
			}
			else
			{
				libmaus::bambam::BamCircularHashCollatingBamDecoder CHCBD(
					bamin.in,
					tmpfilename,excludeflags,
					false, /* put rank */
					hlog,
					sbs
					);
				bamdownsamplerandom(arginfo,CHCBD);
			}
		}
	}
	#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"r","",
			tmpfilename,excludeflags,
			false, /* put rank */
			hlog,sbs
		);
		bamdownsamplerandom(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			false, /* put rank */
			hlog,sbs
		);
		bamdownsamplerandom(arginfo,CHCBD);
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
				
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("p=<[")+libmaus::util::NumberSerialisation::formatNumber(getDefaultProb(),0)+"]>", "probability for keeping read" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("seed=<[]>"), "random seed" ) );
				V.push_back ( std::pair<std::string,std::string> ( "filename=<[stdin]>", "input filename (default: read file from standard input)" ) );
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "ranges=<[]>", "input ranges (bam input only, default: read complete file)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "exclude=<[SECONDARY,SUPPLEMENTARY]>", "exclude alignments matching any of the given flags" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<[0]>", "disable validation of input data" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhlog=<[18]>", "base 2 logarithm of hash table size used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("colsbs=<[")+libmaus::util::NumberSerialisation::formatNumber(128ull*1024*1024,0)+"]>", "size of hash table overflow list in bytes" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("T=<[") + arginfo.getDefaultTmpFileName() + "]>" , "temporary file name" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				std::cerr << "Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" << std::endl;

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
		
			
		bamdownsamplerandom(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
