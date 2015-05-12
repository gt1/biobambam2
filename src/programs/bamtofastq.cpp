/**
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
**/

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

#include <iomanip>

#include <config.h>

#include <libmaus2/aio/PosixFdInputStream.hpp>
#include <libmaus2/aio/PosixFdOutputStream.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamToFastqOutputFileSet.hpp>
#include <libmaus2/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus2/lz/GzipOutputStream.hpp>
#include <libmaus2/util/MemUsage.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <libmaus2/aio/LineSplittingPosixFdOutputStream.hpp>
#include <libmaus2/lz/LineSplittingGzipOutputStream.hpp>

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

bool getDefaultTryOQ()
{
	return 0;
}

std::string getDefaultReadGroupSuffixF(bool const gz, uint64_t const split)
{
	return std::string("_1.fq") + ((gz && split) ? ".gz" : "");
}

std::string getDefaultReadGroupSuffixF2(bool const gz, uint64_t const split)
{
	return std::string("_2.fq") + ((gz && split) ? ".gz" : "");
}

std::string getDefaultReadGroupSuffixO(bool const gz, uint64_t const split)
{
	return std::string("_o1.fq") + ((gz && split) ? ".gz" : "");
}

std::string getDefaultReadGroupSuffixO2(bool const gz, uint64_t const split)
{
	return std::string("_o2.fq") + ((gz && split) ? ".gz" : "");
}

std::string getDefaultReadGroupSuffixS(bool const gz, uint64_t const split)
{
	return std::string("_s.fq") + ((gz && split) ? ".gz" : "");
}

uint64_t getDefaultSplit()
{
	return 0;
}

std::string getDefaultSplitPrefix()
{
	return "bamtofastq_split";
}

struct BamToFastQInputFileStream
{
	std::string const fn;
	int64_t const inputbuffersize;
	libmaus2::aio::PosixFdInputStream::unique_ptr_type CIS;
	std::istream & in;
	
	static int64_t getDefaultBufferSize()
	{
		return -1;
	}
	
	static libmaus2::aio::PosixFdInputStream::unique_ptr_type openFile(std::string const & fn, int64_t const bufsize = getDefaultBufferSize())
	{
		libmaus2::aio::PosixFdInputStream::unique_ptr_type ptr(new libmaus2::aio::PosixFdInputStream(fn,bufsize,0));
		return UNIQUE_PTR_MOVE(ptr);
	}

	static libmaus2::aio::PosixFdInputStream::unique_ptr_type openFile(int const fd, int64_t const bufsize = getDefaultBufferSize())
	{
		libmaus2::aio::PosixFdInputStream::unique_ptr_type ptr(new libmaus2::aio::PosixFdInputStream(fd,bufsize,0));
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	BamToFastQInputFileStream(libmaus2::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  inputbuffersize(arginfo.hasArg("inputbuffersize") ? arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",getDefaultBufferSize()) : -1),
	  CIS(
		(fn != "-") ? openFile(fn) : openFile(STDIN_FILENO)
	  ), in(*CIS) 
	{
        }

	BamToFastQInputFileStream(std::string const & rfn, int64_t const rinputbuffersize = getDefaultBufferSize())
	: fn(rfn),
	  inputbuffersize(rinputbuffersize),
	  CIS(
		(fn != "-") ? openFile(fn) : openFile(STDIN_FILENO)
	), in(*CIS) 
	{
	}
};

enum bamtofastq_conversion_type { bamtofastq_conversion_type_fastq, bamtofastq_conversion_type_fasta, bamtofastq_conversion_type_fastq_try_oq };

template<bamtofastq_conversion_type conversion_type>
uint64_t getSplitMultiplier()
{
	switch ( conversion_type )
	{
		case bamtofastq_conversion_type_fasta:
			return 2;
		case bamtofastq_conversion_type_fastq:
		case bamtofastq_conversion_type_fastq_try_oq:
		default:
			return 4;
	}
}

template<bamtofastq_conversion_type conversion_type>
void bamtofastqNonCollating(libmaus2::util::ArgInfo const & arginfo, libmaus2::bambam::BamAlignmentDecoder & bamdec)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		bamdec.disableValidation();

	libmaus2::timing::RealTimeClock rtc; rtc.start();
	bool const gz = arginfo.getValue<int>("gz",0);
	uint64_t const split = arginfo.getValueUnsignedNumeric("split",getDefaultSplit());
	std::string splitprefix = arginfo.getUnparsedValue("splitprefix",getDefaultSplitPrefix());
	int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION));
	uint32_t const excludeflags = libmaus2::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
	::libmaus2::autoarray::AutoArray<uint8_t> T;
	uint64_t cnt = 0;
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus2::lz::GzipOutputStream::unique_ptr_type Pgzos;
	std::ostream * poutputstream = &std::cout;
	libmaus2::aio::LineSplittingPosixFdOutputStream::unique_ptr_type Psplitos;
	libmaus2::lz::LineSplittingGzipOutputStream::unique_ptr_type Psplitgzos;
	
	if ( split )
	{
		uint64_t const mult = getSplitMultiplier<conversion_type>();
		
		if ( gz )
		{
			libmaus2::lz::LineSplittingGzipOutputStream::unique_ptr_type Tsplitgzos(
				new libmaus2::lz::LineSplittingGzipOutputStream(splitprefix,split*mult,64*1024,level)
			);
			Psplitgzos = UNIQUE_PTR_MOVE(Tsplitgzos);
			poutputstream = Psplitgzos.get();
		
		}
		else
		{
			libmaus2::aio::LineSplittingPosixFdOutputStream::unique_ptr_type Tsplitos(
				new libmaus2::aio::LineSplittingPosixFdOutputStream(splitprefix,split*mult)
			);
			Psplitos = UNIQUE_PTR_MOVE(Tsplitos);
			poutputstream = Psplitos.get();
		}
	}
	else if ( gz )
	{
		libmaus2::lz::GzipOutputStream::unique_ptr_type tPgzos(new libmaus2::lz::GzipOutputStream(std::cout,libmaus2::bambam::BamToFastqOutputFileSet::getGzipBufferSize(),level));
		Pgzos = UNIQUE_PTR_MOVE(tPgzos);
		poutputstream = Pgzos.get();
	}
	
	std::ostream & outputstream = *poutputstream;
		
	while ( bamdec.readAlignment() )
	{
		uint64_t const precnt = cnt++;
		
		if ( ! (algn.getFlags() & excludeflags) )
		{
			uint64_t la;
			
			switch ( conversion_type )
			{
				case bamtofastq_conversion_type_fastq:
					la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(algn.D.begin(),T);
					break;
				case bamtofastq_conversion_type_fasta:
					la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(algn.D.begin(),T);
					break;
				case bamtofastq_conversion_type_fastq_try_oq:
					la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(algn.D.begin(),algn.blocksize,T);
					break;
			}

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

	Psplitos.reset();
	Psplitgzos.reset();
	outputstream.flush();
	if ( Pgzos )
		Pgzos.reset();
	std::cout.flush();
}

void bamtofastqNonCollating(libmaus2::util::ArgInfo const & arginfo, libmaus2::bambam::BamAlignmentDecoder & bamdec, bamtofastq_conversion_type const conversion_type)
{
	switch ( conversion_type )
	{
		case bamtofastq_conversion_type_fasta:
			bamtofastqNonCollating<bamtofastq_conversion_type_fasta>(arginfo,bamdec);
			break;
		case bamtofastq_conversion_type_fastq:					
			bamtofastqNonCollating<bamtofastq_conversion_type_fastq>(arginfo,bamdec);
			break;
		case bamtofastq_conversion_type_fastq_try_oq:
			bamtofastqNonCollating<bamtofastq_conversion_type_fastq_try_oq>(arginfo,bamdec);
			break;
	}
}

void bamtofastqNonCollating(libmaus2::util::ArgInfo const & arginfo)
{
	int64_t const inputbuffersize =
		arginfo.hasArg("inputbuffersize") ? arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",BamToFastQInputFileStream::getDefaultBufferSize()) : -1;

	bool const fasta = arginfo.getValue<int>("fasta",getDefaultFastA());
	bool const tryoq = arginfo.getValue<int>("tryoq",getDefaultTryOQ());	
	bamtofastq_conversion_type const conversion_type = fasta ? bamtofastq_conversion_type_fasta : (tryoq ? bamtofastq_conversion_type_fastq_try_oq : bamtofastq_conversion_type_fastq);

	libmaus2::aio::PosixFdInputStream PFIS(STDIN_FILENO,inputbuffersize);
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false /* put rank */, 0 /* copy stream */, PFIS
		)
	);
	libmaus2::bambam::BamAlignmentDecoder & decoder = decwrapper->getDecoder();

	bamtofastqNonCollating(arginfo,decoder,conversion_type);

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

template<bamtofastq_conversion_type conversion_type>
void bamtofastqCollating(
	libmaus2::util::ArgInfo const & arginfo,
	libmaus2::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	libmaus2::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus2::timing::RealTimeClock rtc; rtc.start();
	::libmaus2::autoarray::AutoArray<uint8_t> T;
	CollateCombs combs;

	bool const outputperreadgroup = arginfo.getValue<unsigned int>("outputperreadgroup",getDefaultOutputPerReadgroup());

	if ( outputperreadgroup )
	{
		bool const gz = arginfo.getValue<int>("gz",0);
		uint64_t const split = arginfo.getValueUnsignedNumeric("split",getDefaultSplit());
		std::string const Fsuffix = arginfo.getUnparsedValue("outputperreadgroupsuffixF",getDefaultReadGroupSuffixF(gz,split));
		std::string const F2suffix = arginfo.getUnparsedValue("outputperreadgroupsuffixF2",getDefaultReadGroupSuffixF2(gz,split));
		std::string const Osuffix = arginfo.getUnparsedValue("outputperreadgroupsuffixO",getDefaultReadGroupSuffixO(gz,split));
		std::string const O2suffix = arginfo.getUnparsedValue("outputperreadgroupsuffixO2",getDefaultReadGroupSuffixO2(gz,split));
		std::string const Ssuffix = arginfo.getUnparsedValue("outputperreadgroupsuffixS",getDefaultReadGroupSuffixS(gz,split));
		
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
		libmaus2::bambam::BamHeader const & header = CHCBD.getHeader();
		
		// get read group vector
		std::vector<libmaus2::bambam::ReadGroup> const & readgroups = header.getReadGroups();
		
		// compute set of read group ids
		std::set<std::string> readgroupsidset;
		for ( uint64_t i = 0; i < readgroups.size(); ++i )
			readgroupsidset.insert(readgroups[i].ID);
			
		// check that read group ids are unique
		if ( readgroupsidset.size() != readgroups.size() )
		{
			libmaus2::exception::LibMausException ex;
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
		int64_t const posixoutbufsize = 256*1024;
		libmaus2::autoarray::AutoArray< ::libmaus2::aio::PosixFdOutputStream::unique_ptr_type > APFOS(numoutputfiles);
		libmaus2::autoarray::AutoArray< libmaus2::lz::GzipOutputStream::unique_ptr_type > AGZOS(numoutputfiles);
		
		libmaus2::autoarray::AutoArray< ::libmaus2::aio::LineSplittingPosixFdOutputStream::unique_ptr_type > ALSPFDOS(numoutputfiles);
		libmaus2::autoarray::AutoArray< ::libmaus2::lz::LineSplittingGzipOutputStream::unique_ptr_type > ALSGZOS(numoutputfiles);
		
		libmaus2::autoarray::AutoArray< std::ostream * > AOS(numoutputfiles);
		libmaus2::autoarray::AutoArray< uint64_t > filefrags(numoutputfiles);
		int const level = std::min(9,std::max(-1,arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION)));
		
		if ( split )
		{
			uint64_t const mult = getSplitMultiplier<conversion_type>();
			
			if ( gz )
			{
				for ( uint64_t i = 0; i < numoutputfiles; ++i )
				{
					::libmaus2::lz::LineSplittingGzipOutputStream::unique_ptr_type tptr(
						new ::libmaus2::lz::LineSplittingGzipOutputStream(outputfilenamevector[i],mult*split,64*1024,level)
					);
					ALSGZOS[i] = UNIQUE_PTR_MOVE(tptr);
					AOS[i] = ALSGZOS[i].get();
				}	
			}
			else
			{
				for ( uint64_t i = 0; i < numoutputfiles; ++i )
				{
					::libmaus2::aio::LineSplittingPosixFdOutputStream::unique_ptr_type tptr(
						new ::libmaus2::aio::LineSplittingPosixFdOutputStream(outputfilenamevector[i],mult*split)
					);
					ALSPFDOS[i] = UNIQUE_PTR_MOVE(tptr);
					AOS[i] = ALSPFDOS[i].get();
				}				
			}
		}
		else
		{
			for ( uint64_t i = 0; i < numoutputfiles; ++i )
			{
				::libmaus2::aio::PosixFdOutputStream::unique_ptr_type tptr(
					new ::libmaus2::aio::PosixFdOutputStream(outputfilenamevector[i],posixoutbufsize)
				);
				APFOS[i] = UNIQUE_PTR_MOVE(tptr);
				
				if ( gz )
				{
					libmaus2::lz::GzipOutputStream::unique_ptr_type tPgzos(new libmaus2::lz::GzipOutputStream(*(APFOS[i]),libmaus2::bambam::BamToFastqOutputFileSet::getGzipBufferSize(),level));
					AGZOS[i] = UNIQUE_PTR_MOVE(tPgzos);
					AOS[i] = AGZOS[i].get();
				}
				else
				{
					AOS[i] = APFOS[i].get();
				}
			}
		}

		while ( (ob = CHCBD.process()) )
		{
			uint64_t const precnt = cnt;
			
			int64_t const rg = ob->Da ?
				header.getReadGroupId(libmaus2::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Da,ob->blocksizea))
				:
				-1
				;
			int64_t const rgfbase = rg + 1;
			assert ( static_cast<int64_t>(rgfbase) < static_cast<int64_t>(readgroups.size() + 1) );
			uint64_t const rgfshift = rgfbase * filesperrg;
			
			if ( ob->fpair )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}

				AOS[rgfshift + Fmap]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Fmap]++;
				
				uint64_t lb;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Db,T);
						break;
					case bamtofastq_conversion_type_fastq:
						lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Db,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Db,ob->blocksizeb,T);
						break;
				}

				AOS[rgfshift + F2map]->write(reinterpret_cast<char const *>(T.begin()),lb);
				filefrags[rgfshift + F2map]++;

				combs.pairs += 1;
				cnt += 2;
				bcnt += (la+lb);
			}
			else if ( ob->fsingle )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}

				AOS[rgfshift + Smap]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Smap]++;

				combs.single += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan1 )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}

				AOS[rgfshift + Omap]->write(reinterpret_cast<char const *>(T.begin()),la);
				filefrags[rgfshift + Omap]++;

				combs.orphans1 += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan2 )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}

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

		if ( split )
		{
			for ( uint64_t i = 0; i < numoutputfiles; ++i )
			{
				ALSPFDOS[i].reset();
				ALSGZOS[i].reset();
			}
		}
		else
		{
			// close files
			for ( uint64_t i = 0; i < numoutputfiles; ++i )
			{
				if ( AGZOS[i] )
					AGZOS[i].reset();
				APFOS[i].reset();
			
				// remove empty files
				if ( ! filefrags[i] )
					remove(outputfilenamevector[i].c_str());
			}
		}
	}
	else
	{		
		libmaus2::bambam::BamToFastqOutputFileSet OFS(arginfo);
		
		while ( (ob = CHCBD.process()) )
		{
			uint64_t const precnt = cnt;
			
			if ( ob->fpair )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}
				OFS.Fout.write(reinterpret_cast<char const *>(T.begin()),la);

				uint64_t lb;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Db,T);
						break;
					case bamtofastq_conversion_type_fastq:
						lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Db,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Db,ob->blocksizeb,T);
						break;
				}
				OFS.F2out.write(reinterpret_cast<char const *>(T.begin()),lb);

				combs.pairs += 1;
				cnt += 2;
				bcnt += (la+lb);
			}
			else if ( ob->fsingle )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}
				OFS.Sout.write(reinterpret_cast<char const *>(T.begin()),la);

				combs.single += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan1 )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}
				OFS.Oout.write(reinterpret_cast<char const *>(T.begin()),la);

				combs.orphans1 += 1;
				cnt += 1;
				bcnt += (la);
			}
			else if ( ob->forphan2 )
			{
				uint64_t la;
				switch ( conversion_type )
				{
					case bamtofastq_conversion_type_fasta:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastA(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQ(ob->Da,T);
						break;
					case bamtofastq_conversion_type_fastq_try_oq:
						la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQTryOQ(ob->Da,ob->blocksizea,T);
						break;
				}
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

void bamtofastqCollating(
	libmaus2::util::ArgInfo const & arginfo,
	libmaus2::bambam::CircularHashCollatingBamDecoder & CHCBD,
	bamtofastq_conversion_type const conversion_type
)
{
	switch ( conversion_type )
	{
		case bamtofastq_conversion_type_fasta:
			bamtofastqCollating<bamtofastq_conversion_type_fasta>(arginfo,CHCBD);
			break;
		case bamtofastq_conversion_type_fastq:
			bamtofastqCollating<bamtofastq_conversion_type_fastq>(arginfo,CHCBD);
			break;
		case bamtofastq_conversion_type_fastq_try_oq:
			bamtofastqCollating<bamtofastq_conversion_type_fastq_try_oq>(arginfo,CHCBD);
			break;
	}
}

void bamtofastqCollating(libmaus2::util::ArgInfo const & arginfo)
{
	// set up for temp file
	libmaus2::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());	
	libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	
	// exclude flags for collation
	uint32_t const excludeflags = libmaus2::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	
	// input format
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());
	// input filename
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");

	// input buffer size (if input is not via io_lib)
	int64_t const inputbuffersize =
		arginfo.hasArg("inputbuffersize") ? arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",BamToFastQInputFileStream::getDefaultBufferSize()) : -1;
	
	bool const fasta = arginfo.getValue<int>("fasta",getDefaultFastA());
	bool const tryoq = arginfo.getValue<int>("tryoq",getDefaultTryOQ());	
	bamtofastq_conversion_type const conversion_type = fasta ? bamtofastq_conversion_type_fasta : (tryoq ? bamtofastq_conversion_type_fastq_try_oq : bamtofastq_conversion_type_fastq);

	// table size
	unsigned int const hlog = arginfo.getValue<unsigned int>("colhlog",18);
	// overflow list length
	uint64_t const sbs = arginfo.getValueUnsignedNumeric<uint64_t>("colsbs",32ull*1024ull*1024ull);

	// adapter for standard input
	libmaus2::aio::PosixFdInputStream PFIS(STDIN_FILENO,inputbuffersize);
	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false /* put rank */, 0 /* copy stream */, PFIS
		)
	);
	libmaus2::bambam::BamAlignmentDecoder & decoder = decwrapper->getDecoder();
	// collator
	libmaus2::bambam::CircularHashCollatingBamDecoder CHCBD(decoder,tmpfilename,excludeflags,hlog,sbs);
	bamtofastqCollating(arginfo,CHCBD,conversion_type);
	
	std::cout.flush();
}

void bamtofastqCollatingRanking(
	libmaus2::util::ArgInfo const & arginfo,
	libmaus2::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	libmaus2::bambam::BamToFastqOutputFileSet OFS(arginfo);

	libmaus2::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus2::timing::RealTimeClock rtc; rtc.start();
	::libmaus2::autoarray::AutoArray<uint8_t> T;
	CollateCombs combs;
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair )
		{
			uint64_t const ranka = libmaus2::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			uint64_t const rankb = libmaus2::bambam::BamAlignmentDecoderBase::getRank(ob->Db,ob->blocksizeb);
			uint64_t const la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,rankb,T);		
			OFS.Fout.write(reinterpret_cast<char const *>(T.begin()),la);
			uint64_t lb = libmaus2::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Db,ranka,rankb,T);
			OFS.F2out.write(reinterpret_cast<char const *>(T.begin()),lb);

			combs.pairs += 1;
			cnt += 2;
			bcnt += (la+lb);
		}
		else if ( ob->fsingle )
		{
			uint64_t const ranka = libmaus2::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			uint64_t la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,ranka,T);
			OFS.Sout.write(reinterpret_cast<char const *>(T.begin()),la);

			combs.single += 1;
			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan1 )
		{
			uint64_t const ranka = libmaus2::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			uint64_t la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,ranka,T);
			OFS.Oout.write(reinterpret_cast<char const *>(T.begin()),la);

			combs.orphans1 += 1;
			cnt += 1;
			bcnt += (la);
		}
		else if ( ob->forphan2 )
		{
			uint64_t const ranka = libmaus2::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			uint64_t la = libmaus2::bambam::BamAlignmentDecoderBase::putFastQRanks(ob->Da,ranka,ranka,T);
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

void bamtofastqCollatingRanking(libmaus2::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus2::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus2::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const inputformat = arginfo.getValue<std::string>("inputformat",getDefaultInputFormat());
	std::string const inputfilename = arginfo.getValue<std::string>("filename","-");
	uint64_t const numthreads = arginfo.getValue<uint64_t>("threads",0);
	int64_t const inputbuffersize =
		arginfo.hasArg("inputbuffersize") ? arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",BamToFastQInputFileStream::getDefaultBufferSize()) : -1;

	unsigned int const hlog = arginfo.getValue<unsigned int>("colhlog",18);
	uint64_t const sbs = arginfo.getValueUnsignedNumeric<uint64_t>("colsbs",32ull*1024ull*1024ull);

	if ( inputformat == "bam" )
	{
		BamToFastQInputFileStream bamin(inputfilename,inputbuffersize);

		if ( numthreads > 0 )
		{
			libmaus2::bambam::BamParallelCircularHashCollatingBamDecoder CHCBD(
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
			libmaus2::bambam::BamCircularHashCollatingBamDecoder CHCBD(
				bamin.in,
				tmpfilename,excludeflags,
				true, /* put rank */
				hlog,
				sbs
				);
			bamtofastqCollatingRanking(arginfo,CHCBD);
		}
	}
	#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
	else if ( inputformat == "sam" )
	{
		libmaus2::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"r","",
			tmpfilename,excludeflags,
			true, /* put rank */
			hlog,sbs
		);
		bamtofastqCollatingRanking(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus2::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			true, /* put rank */
			hlog,sbs
		);
		bamtofastqCollatingRanking(arginfo,CHCBD);
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

void bamtofastq(libmaus2::util::ArgInfo const & arginfo)
{
	// if range is requested than check whether they are supported
	if ( 
		arginfo.hasArg("ranges") && arginfo.getValue("inputformat", getDefaultInputFormat()) != "bam"
		&&		
		arginfo.hasArg("ranges") && arginfo.getValue("inputformat", getDefaultInputFormat()) != "cram"
	)
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "ranges are only supported for inputformat=bam" << std::endl;
		se.finish();
		throw se;
	}
	// ranges are not supported via stdin
	if ( arginfo.hasArg("ranges") && ((!arginfo.hasArg("filename")) || arginfo.getValue<std::string>("filename","-") == "-") )
	{
		libmaus2::exception::LibMausException se;
		se.getStream() << "ranges are not supported for reading via standard input" << std::endl;
		se.finish();
		throw se;
	}
	// ranges are not supported for collate>1 (ranks do not make sense)
	if ( arginfo.hasArg("ranges") && arginfo.getValue<uint64_t>("collate",1) > 1 )
	{
		libmaus2::exception::LibMausException se;
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
			libmaus2::exception::LibMausException se;
			se.getStream() << "unknown collate argument " << arginfo.getValue<uint64_t>("collate",1) << std::endl;
			se.finish();
			throw se;
		}
	}
}

#if defined(LIBMAUS2_HAVE_IRODS)
#include <libmaus2/irods/IRodsInputStreamFactory.hpp>
#endif

int main(int argc, char * argv[])
{
	try
	{
		#if defined(LIBMAUS2_HAVE_IRODS)
		libmaus2::irods::IRodsInputStreamFactory::registerHandler();
		#endif

		libmaus2::timing::RealTimeClock rtc; rtc.start();
		
		::libmaus2::util::ArgInfo arginfo(argc,argv);
		
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
				
				bool const gz = arginfo.getValue<unsigned int>("gz",0);
				uint64_t const split = arginfo.getValueUnsignedNumeric("split",0);
				
				V.push_back ( std::pair<std::string,std::string> ( "F=<[stdout]>", "matched pairs first mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "F2=<[stdout]>", "matched pairs second mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "S=<[stdout]>", "single end" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "unmatched pairs first mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O2=<[stdout]>", "unmatched pairs second mates" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collate=<[1]>", "collate pairs" ) );
				V.push_back ( std::pair<std::string,std::string> ( "combs=<[0]>", "print some counts after collation based processing" ) );
				V.push_back ( std::pair<std::string,std::string> ( "filename=<[stdin]>", "input filename (default: read file from standard input)" ) );
				#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif
				#if defined(BIOBAMBAM_LIBMAUS2_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( "ranges=<[]>", "input ranges (bam and cram input only, default: read complete file)" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "ranges=<[]>", "input ranges (bam input only, default: read complete file)" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "exclude=<[SECONDARY,SUPPLEMENTARY]>", "exclude alignments matching any of the given flags" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<[0]>", "disable validation of input data" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhlog=<[18]>", "base 2 logarithm of hash table size used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("colsbs=<[")+libmaus2::util::NumberSerialisation::formatNumber(32ull*1024*1024,0)+"]>", "size of hash table overflow list in bytes" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("T=<[") + arginfo.getDefaultTmpFileName() + "]>" , "temporary file name" ) );
				V.push_back ( std::pair<std::string,std::string> ( "gz=<[0]>", "compress output streams in gzip format (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", std::string("compression setting if gz=1 (") + libmaus2::bambam::BamBlockWriterBaseFactory::getLevelHelpText() + std::string(")")  ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("fasta=<[")+libmaus2::util::NumberSerialisation::formatNumber(getDefaultFastA(),0)+"]>", "output FastA instead of FastQ" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputbuffersize=<["+::biobambam2::Licensing::formatNumber(BamToFastQInputFileStream::getDefaultBufferSize())+"]>", "size of input buffer" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroup=<["+::biobambam2::Licensing::formatNumber(getDefaultOutputPerReadgroup())+"]>", "split output per read group (for collate=1 only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputdir=<>", "directory for output if outputperreadgroup=1 (default: current directory)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixF=<["+getDefaultReadGroupSuffixF(gz,split)+"]>", "suffix for F category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixF2=<["+getDefaultReadGroupSuffixF2(gz,split)+"]>", "suffix for F2 category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixO=<["+getDefaultReadGroupSuffixO(gz,split)+"]>", "suffix for O category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixO2=<["+getDefaultReadGroupSuffixO2(gz,split)+"]>", "suffix for O2 category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputperreadgroupsuffixS=<["+getDefaultReadGroupSuffixS(gz,split)+"]>", "suffix for S category when outputperreadgroup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("tryoq=<[") + ::biobambam2::Licensing::formatNumber(getDefaultTryOQ()) + "]>", "use OQ field instead of quality field if present (collate={0,1} only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("split=<[")+libmaus2::util::NumberSerialisation::formatNumber(getDefaultSplit(),0)+"]>", "split named output files into chunks of this amount of reads (0: do not split)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("splitprefix=<[")+getDefaultSplitPrefix()+"]>", "file name prefix if collate=0 and split>0" ) );
				
				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				std::cerr << "Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" << std::endl;

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
		
		if ( arginfo.hasArg("filename") )
		{
			std::string const fn = arginfo.getUnparsedValue("filename",std::string());
			arginfo.argmap["I"] = fn;
			arginfo.argmultimap.insert(std::pair<std::string,std::string>(std::string("I"),fn));
		}
		if ( arginfo.hasArg("I") && !arginfo.hasArg("filename") )
		{
			std::string const fn = arginfo.getUnparsedValue("I",std::string());		
			arginfo.argmap["filename"] = fn;
			arginfo.argmultimap.insert(std::pair<std::string,std::string>(std::string("filename"),fn));
		}
		if ( arginfo.hasArg("ranges") )
		{
			std::string const range = arginfo.getUnparsedValue("ranges",std::string());
			arginfo.argmap["range"] = range;
			arginfo.argmultimap.insert(std::pair<std::string,std::string>(std::string("range"),range));
		}
		if ( arginfo.hasArg("range") && !arginfo.hasArg("ranges") )
		{
			std::string const range = arginfo.getUnparsedValue("range",std::string());		
			arginfo.argmap["ranges"] = range;
			arginfo.argmultimap.insert(std::pair<std::string,std::string>(std::string("ranges"),range));
		}
			
		bamtofastq(arginfo);
		
		std::cerr << "[V] " << libmaus2::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
