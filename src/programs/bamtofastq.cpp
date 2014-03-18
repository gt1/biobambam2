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
#include <libmaus/lz/GzipOutputStream.hpp>
#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/aio/PosixFdOutputStream.hpp>

static std::string getDefaultInputFormat()
{
	return "bam";
}

bool getDefaultFastA()
{
	return 0;
}

bool getDefaultOutputPerReadgroup()
{
	return 0;
}

std::string getDefaultReadGroupSuffixF()
{
	return "_1.fq";
}

std::string getDefaultReadGroupSuffixF2()
{
	return "_2.fq";
}

std::string getDefaultReadGroupSuffixO()
{
	return "_o1.fq";
}

std::string getDefaultReadGroupSuffixO2()
{
	return "_o2.fq";
}

std::string getDefaultReadGroupSuffixS()
{
	return "_s.fq";
}

struct BamToFastQInputFileStream
{
	std::string const fn;
	uint64_t const inputbuffersize;
	libmaus::aio::PosixFdInputStream::unique_ptr_type CIS;
	std::istream & in;
	
	static uint64_t getDefaultBufferSize()
	{
		return 2*1024*1024;
	}
	
	static libmaus::aio::PosixFdInputStream::unique_ptr_type openFile(std::string const & fn, uint64_t const bufsize = getDefaultBufferSize())
	{
		libmaus::aio::PosixFdInputStream::unique_ptr_type ptr(new libmaus::aio::PosixFdInputStream(fn,bufsize,0));
		return UNIQUE_PTR_MOVE(ptr);
	}

	static libmaus::aio::PosixFdInputStream::unique_ptr_type openFile(int const fd, uint64_t const bufsize = getDefaultBufferSize())
	{
		libmaus::aio::PosixFdInputStream::unique_ptr_type ptr(new libmaus::aio::PosixFdInputStream(fd,bufsize,0));
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	BamToFastQInputFileStream(libmaus::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  inputbuffersize(arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",getDefaultBufferSize())),
	  CIS(
		(fn != "-") ? openFile(fn) : openFile(STDIN_FILENO)
	  ), in(*CIS) {}

	BamToFastQInputFileStream(std::string const & rfn, uint64_t const rinputbuffersize = getDefaultBufferSize())
	: fn(rfn),
	  inputbuffersize(rinputbuffersize), 
	  CIS(
		(fn != "-") ? openFile(fn) : openFile(STDIN_FILENO)
	), in(*CIS) {}
};

void bamtofastqNonCollating(libmaus::util::ArgInfo const & arginfo, libmaus::bambam::BamAlignmentDecoder & bamdec)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		bamdec.disableValidation();

	libmaus::timing::RealTimeClock rtc; rtc.start();
	bool const gz = arginfo.getValue<int>("gz",0);
	int const level = std::min(9,std::max(-1,arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION)));
	bool const fasta = arginfo.getValue<int>("fasta",getDefaultFastA());
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus::bambam::BamAlignment const & algn = bamdec.getAlignment();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	uint64_t cnt = 0;
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus::lz::GzipOutputStream::unique_ptr_type Pgzos;
	
	if ( gz )
	{
		libmaus::lz::GzipOutputStream::unique_ptr_type tPgzos(new libmaus::lz::GzipOutputStream(std::cout,libmaus::bambam::BamToFastqOutputFileSet::getGzipBufferSize(),level));
		Pgzos = UNIQUE_PTR_MOVE(tPgzos);
	}
	
	std::ostream & outputstream = gz ? *Pgzos : std::cout;
		
	while ( bamdec.readAlignment() )
	{
		uint64_t const precnt = cnt++;
		
		if ( ! (algn.getFlags() & excludeflags) )
		{
			uint64_t la = 
				fasta ?
				libmaus::bambam::BamAlignmentDecoderBase::putFastA(algn.D.begin(),T)
				:
				libmaus::bambam::BamAlignmentDecoderBase::putFastQ(algn.D.begin(),T)				
				;
			outputstream.write(reinterpret_cast<char const *>(T.begin()),la);
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

	outputstream.flush();
	if ( Pgzos )
		Pgzos.reset();
	std::cout.flush();
}

void bamtofastqNonCollating(libmaus::util::ArgInfo const & arginfo)
{
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	uint64_t const numthreads = arginfo.getValue<uint64_t>("threads",0);
	uint64_t const inputbuffersize = arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",BamToFastQInputFileStream::getDefaultBufferSize());
	
	if ( inputformat == "bam" )
	{
		if ( arginfo.hasArg("ranges") )
		{
			libmaus::bambam::BamRangeDecoder bamdec(inputfilename,arginfo.getUnparsedValue("ranges",""));
			bamtofastqNonCollating(arginfo,bamdec);
		}
		else
		{
			BamToFastQInputFileStream bamin(inputfilename,inputbuffersize);
			
			if ( numthreads > 0 )
			{
				libmaus::bambam::BamParallelDecoderWrapper bamdecwrap(bamin.in,numthreads);
				bamtofastqNonCollating(arginfo,bamdecwrap.getDecoder());	
			}
			else
			{
				libmaus::bambam::BamDecoder bamdec(bamin.in);
				bamtofastqNonCollating(arginfo,bamdec);
			}
		}
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

struct CollateCombs
{
	uint64_t pairs;
	uint64_t orphans1;
	uint64_t orphans2;
	uint64_t single;
	uint64_t alignments;

	CollateCombs()
	: pairs(0), orphans1(0), orphans2(0), single(0), alignments(0)
	{
	
	}
};

std::ostream & operator<<(std::ostream & out, CollateCombs const & combs)
{
	out << "[C]\tAlignments:\t" << combs.alignments << std::endl;
	out << "[C]\tComplete pairs:\t" << combs.pairs << std::endl;
	out << "[C]\tSingle:\t" << combs.single << std::endl;
	out << "[C]\tOrphans:\t" << (combs.orphans1 + combs.orphans2) 
		<< "\t" << combs.orphans1 << "\t" << combs.orphans2 << std::endl;
	return out;
};

void bamtofastqCollating(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	bool const fasta = arginfo.getValue<int>("fasta",getDefaultFastA());


	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	CollateCombs combs;

	bool const outputperreadgroup = arginfo.getValue<unsigned int>("outputperreadgroup",getDefaultOutputPerReadgroup());

	if ( outputperreadgroup )
	{
		std::string const Fsuffix = arginfo.getUnparsedValue("outputperreadgroupsuffixF",getDefaultReadGroupSuffixF());
		std::string const F2suffix = arginfo.getUnparsedValue("outputperreadgroupsuffixF2",getDefaultReadGroupSuffixF2());
		std::string const Osuffix = arginfo.getUnparsedValue("outputperreadgroupsuffixO",getDefaultReadGroupSuffixO());
		std::string const O2suffix = arginfo.getUnparsedValue("outputperreadgroupsuffixO2",getDefaultReadGroupSuffixO2());
		std::string const Ssuffix = arginfo.getUnparsedValue("outputperreadgroupsuffixS",getDefaultReadGroupSuffixS());
		
		// collect set of all suffixes
		std::set<std::string> suffixset;
		suffixset.insert(Fsuffix);
		suffixset.insert(F2suffix);
		suffixset.insert(Osuffix);
		suffixset.insert(O2suffix);
		suffixset.insert(Ssuffix);

		// assign id to each suffix
		std::map<std::string,uint64_t> suffixmap;
		std::vector<std::string> suffixremap;
		for ( std::set<std::string>::const_iterator ita = suffixset.begin(); ita != suffixset.end(); ++ita )
		{
			uint64_t const id = suffixmap.size();
			suffixmap[*ita] = id;
			suffixremap.push_back(*ita);
		}
		uint64_t const Fmap = suffixmap.find(Fsuffix)->second;
		uint64_t const F2map = suffixmap.find(F2suffix)->second;
		uint64_t const Omap = suffixmap.find(Osuffix)->second;
		uint64_t const O2map = suffixmap.find(O2suffix)->second;
		uint64_t const Smap = suffixmap.find(Ssuffix)->second;
		
		// get bam header
		libmaus::bambam::BamHeader const & header = CHCBD.getHeader();
		
		// get read group vector
		std::vector<libmaus::bambam::ReadGroup> const & readgroups = header.getReadGroups();
		
		// compute set of read group ids
		std::set<std::string> readgroupsidset;
		for ( uint64_t i = 0; i < readgroups.size(); ++i )
			readgroupsidset.insert(readgroups[i].ID);
			
		// check that read group ids are unique
		if ( readgroupsidset.size() != readgroups.size() )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "read group ids are not unique." << std::endl;
			ex.finish();
			throw ex;
		}
			
		// construct id for default read group
		std::string defaultidprefix = "default";
		std::string defaultid = defaultidprefix;
		uint64_t defaultidex = 0;
		while ( readgroupsidset.find(defaultid) != readgroupsidset.end() )
		{
			std::ostringstream ostr;
			ostr << defaultidprefix << "_" << (defaultidex++);
			defaultid = ostr.str();
		}
		assert ( readgroupsidset.find(defaultid) == readgroupsidset.end() );
		
		// number of output files
		uint64_t const filesperrg = suffixset.size();
		uint64_t const numoutputfiles = (readgroups.size()+1) * filesperrg;
		
		// output directory
		std::string outputdir = arginfo.getUnparsedValue("outputdir","");
		if ( outputdir.size() )
			outputdir += "/";
		
		// construct output file names
		std::vector<std::string> outputfilenamevector;	
		for ( std::set<std::string>::const_iterator ita = suffixset.begin(); ita != suffixset.end(); ++ita )
			outputfilenamevector.push_back(outputdir + defaultid + *ita);
		for ( uint64_t i = 0; i < readgroups.size(); ++i )
			for ( std::set<std::string>::const_iterator ita = suffixset.begin(); ita != suffixset.end(); ++ita )
				outputfilenamevector.push_back(outputdir + readgroups[i].ID + *ita);
				
		assert ( outputfilenamevector.size() == numoutputfiles );
		uint64_t const posixoutbufsize = 256*1024;
		libmaus::autoarray::AutoArray< ::libmaus::aio::PosixFdOutputStream::unique_ptr_type > APFOS(numoutputfiles);
		libmaus::autoarray::AutoArray< libmaus::lz::GzipOutputStream::unique_ptr_type > AGZOS(numoutputfiles);
		libmaus::autoarray::AutoArray< std::ostream * > AOS(numoutputfiles);
		libmaus::autoarray::AutoArray< uint64_t > filefrags(numoutputfiles);
		bool const gz = arginfo.getValue<int>("gz",0);
		int const level = std::min(9,std::max(-1,arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION)));
		for ( uint64_t i = 0; i < numoutputfiles; ++i )
		{
			::libmaus::aio::PosixFdOutputStream::unique_ptr_type tptr(
				new ::libmaus::aio::PosixFdOutputStream(outputfilenamevector[i],posixoutbufsize)
			);
			APFOS[i] = UNIQUE_PTR_MOVE(tptr);
			
			if ( gz )
			{
				libmaus::lz::GzipOutputStream::unique_ptr_type tPgzos(new libmaus::lz::GzipOutputStream(*(APFOS[i]),libmaus::bambam::BamToFastqOutputFileSet::getGzipBufferSize(),level));
				AGZOS[i] = UNIQUE_PTR_MOVE(tPgzos);
				AOS[i] = AGZOS[i].get();
			}
			else
			{
				AOS[i] = APFOS[i].get();
			}
		}

		while ( (ob = CHCBD.process()) )
		{
			uint64_t const precnt = cnt;
			
			int64_t const rg = ob->Da ?
				header.getReadGroupId(libmaus::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Da,ob->blocksizea))
				:
				-1
				;
			int64_t const rgfbase = rg + 1;
			assert ( rgfbase < readgroups.size() + 1 );
			uint64_t const rgfshift = rgfbase * filesperrg;
			
			if ( ob->fpair )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
				;
				AOS[rgfshift + Fmap]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Fmap]++;
				uint64_t lb = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Db,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Db,T)
					;
				AOS[rgfshift + F2map]->write(reinterpret_cast<char const *>(T.begin()),lb);
				filefrags[rgfshift + F2map]++;

				combs.pairs += 1;
				cnt += 2;
				bcnt += (la+lb);
			}
			else if ( ob->fsingle )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
					;
				AOS[rgfshift + Smap]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Smap]++;

				combs.single += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan1 )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
					;
				AOS[rgfshift + Omap]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Omap]++;

				combs.orphans1 += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan2 )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
					;
				AOS[rgfshift + O2map]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Omap]++;

				combs.orphans2 += 1;
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

		
		// close files
		for ( uint64_t i = 0; i < numoutputfiles; ++i )
		{
			if ( AGZOS[i] )
				AGZOS[i].reset();
			APFOS[i].reset();
			
			if ( ! filefrags[i] )
				remove(outputfilenamevector[i].c_str());
		}
	}
	else
	{		
		libmaus::bambam::BamToFastqOutputFileSet OFS(arginfo);
		
		while ( (ob = CHCBD.process()) )
		{
			uint64_t const precnt = cnt;
			
			if ( ob->fpair )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
				;
				OFS.Fout.write(reinterpret_cast<char const *>(T.begin()),la);
				uint64_t lb = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Db,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Db,T)
					;
				OFS.F2out.write(reinterpret_cast<char const *>(T.begin()),lb);

				combs.pairs += 1;
				cnt += 2;
				bcnt += (la+lb);
			}
			else if ( ob->fsingle )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
					;
				OFS.Sout.write(reinterpret_cast<char const *>(T.begin()),la);

				combs.single += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan1 )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
					;
				OFS.Oout.write(reinterpret_cast<char const *>(T.begin()),la);

				combs.orphans1 += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan2 )
			{
				uint64_t la = 
					fasta
					?
					libmaus::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T)
					:
					libmaus::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T)
					;
				OFS.O2out.write(reinterpret_cast<char const *>(T.begin()),la);

				combs.orphans2 += 1;
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
	}

	std::cerr << "[V] " << cnt << std::endl;
	
	combs.alignments = CHCBD.getRank();

	if ( arginfo.getValue<unsigned int>("combs",0) )
		std::cerr << combs;
}

void bamtofastqCollating(libmaus::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	uint64_t const numthreads = arginfo.getValue<uint64_t>("threads",0);
	uint64_t const inputbuffersize = arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",BamToFastQInputFileStream::getDefaultBufferSize());

	unsigned int const hlog = arginfo.getValue<unsigned int>("colhlog",18);
	uint64_t const sbs = arginfo.getValueUnsignedNumeric<uint64_t>("colsbs",32ull*1024ull*1024ull);

	if ( inputformat == "bam" )
	{
		if ( arginfo.hasArg("ranges") )
		{
			libmaus::bambam::BamRangeCircularHashCollatingBamDecoder CHCBD(inputfilename,arginfo.getUnparsedValue("ranges",""),tmpfilename,excludeflags,false,hlog,sbs);
			bamtofastqCollating(arginfo,CHCBD);
		}
		else
		{
			BamToFastQInputFileStream bamin(inputfilename,inputbuffersize);

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
				bamtofastqCollating(arginfo,CHCBD);
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
				bamtofastqCollating(arginfo,CHCBD);
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
		bamtofastqCollating(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			false, /* put rank */
			hlog,sbs
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

void bamtofastqCollatingRanking(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	libmaus::bambam::BamToFastqOutputFileSet OFS(arginfo);

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	CollateCombs combs;
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			uint64_t const rankb = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Db,ob->blocksizeb);
			uint64_t const la = libmaus::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,rankb,T);		
			OFS.Fout.write(reinterpret_cast<char const *>(T.begin()),la);
			uint64_t lb = libmaus::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Db,ranka,rankb,T);
			OFS.F2out.write(reinterpret_cast<char const *>(T.begin()),lb);

			combs.pairs += 1;
			cnt += 2;
			bcnt += (la+lb);
		}
		else if ( ob->fsingle )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,ranka,T);
			OFS.Sout.write(reinterpret_cast<char const *>(T.begin()),la);

			combs.single += 1;
			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan1 )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,ranka,T);
			OFS.Oout.write(reinterpret_cast<char const *>(T.begin()),la);

			combs.orphans1 += 1;
			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan2 )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			uint64_t la = libmaus::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,ranka,T);
			OFS.O2out.write(reinterpret_cast<char const *>(T.begin()),la);

			combs.orphans2 += 1;
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

	combs.alignments = CHCBD.getRank();

	if ( arginfo.getValue<unsigned int>("combs",0) )
		std::cerr << combs;
}

void bamtofastqCollatingRanking(libmaus::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	uint64_t const numthreads = arginfo.getValue<uint64_t>("threads",0);
	uint64_t const inputbuffersize = arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",BamToFastQInputFileStream::getDefaultBufferSize());

	unsigned int const hlog = arginfo.getValue<unsigned int>("colhlog",18);
	uint64_t const sbs = arginfo.getValueUnsignedNumeric<uint64_t>("colsbs",32ull*1024ull*1024ull);

	if ( inputformat == "bam" )
	{
		BamToFastQInputFileStream bamin(inputfilename,inputbuffersize);

		if ( numthreads > 0 )
		{
			libmaus::bambam::BamParallelCircularHashCollatingBamDecoder CHCBD(
				bamin.in,
				numthreads,
				tmpfilename,excludeflags,
				true, /* put rank */
				hlog,
				sbs
				);
			bamtofastqCollatingRanking(arginfo,CHCBD);
		}
		else
		{
			libmaus::bambam::BamCircularHashCollatingBamDecoder CHCBD(
				bamin.in,
				tmpfilename,excludeflags,
				true, /* put rank */
				hlog,
				sbs
				);
			bamtofastqCollatingRanking(arginfo,CHCBD);
		}
	}
	#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"r","",
			tmpfilename,excludeflags,
			true, /* put rank */
			hlog,sbs
		);
		bamtofastqCollatingRanking(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			true, /* put rank */
			hlog,sbs
		);
		bamtofastqCollatingRanking(arginfo,CHCBD);
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

	switch ( arginfo.getValue<uint64_t>("collate",1) )
	{
		case 0:
			bamtofastqNonCollating(arginfo);
			break;
		case 1:
			bamtofastqCollating(arginfo);
			break;
		case 2:		
			bamtofastqCollatingRanking(arginfo);
			break;
		default:
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "unknown collate argument " << arginfo.getValue<uint64_t>("collate",1) << std::endl;
			se.finish();
			throw se;
		}
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
				
				V.push_back ( std::pair<std::string,std::string> ( "F=<[stdout]>", "matched pairs first mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "F2=<[stdout]>", "matched pairs second mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "S=<[stdout]>", "single end" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "unmatched pairs first mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O2=<[stdout]>", "unmatched pairs second mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collate=<[1]>", "collate pairs" ) );
				V.push_back ( std::pair<std::string,std::string> ( "combs=<[0]>", "print some counts after collation based processing" ) );
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
				V.push_back ( std::pair<std::string,std::string> ( std::string("colsbs=<[")+libmaus::util::NumberSerialisation::formatNumber(32ull*1024*1024,0)+"]>", "size of hash table overflow list in bytes" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("T=<[") + arginfo.getDefaultTmpFileName() + "]>" , "temporary file name" ) );
				V.push_back ( std::pair<std::string,std::string> ( "gz=<[0]>", "compress output streams in gzip format (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", "compression setting if gz=1 (default: -1, zlib default settings)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("fasta=<[")+libmaus::util::NumberSerialisation::formatNumber(getDefaultFastA(),0)+"]>", "output FastA instead of FastQ" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputbuffersize=<["+::biobambam::Licensing::formatNumber(BamToFastQInputFileStream::getDefaultBufferSize())+"]>", "size of input buffer" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroup=<["+::biobambam::Licensing::formatNumber(getDefaultOutputPerReadgroup())+"]>", "split output per read group (for collate=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputdir=<>", "directory for output (default: in current directory)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixF=<["+getDefaultReadGroupSuffixF()+"]>", "suffix for F category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixF2=<["+getDefaultReadGroupSuffixF2()+"]>", "suffix for F2 category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixO=<["+getDefaultReadGroupSuffixO()+"]>", "suffix for O category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixO2=<["+getDefaultReadGroupSuffixO2()+"]>", "suffix for O2 category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixS=<["+getDefaultReadGroupSuffixS()+"]>", "suffix for S category when outputperreadgroup=1" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				std::cerr << "Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" << std::endl;

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
		
			
		bamtofastq(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
