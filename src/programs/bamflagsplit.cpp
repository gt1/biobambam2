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
#include "config.h"

#include <iostream>
#include <queue>

#include <libmaus/aio/CheckedOutputStream.hpp>

#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamAlignmentNameComparator.hpp>
#include <libmaus/bambam/BamAlignmentPosComparator.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamEntryContainer.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/aio/PosixFdOutputStream.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GetObject.hpp>
#include <libmaus/util/PutObject.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }

#include <biobambam/BamBamConfig.hpp>

#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
#include <libmaus/bambam/ScramDecoder.hpp>
#endif

#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static bool getDefaultDisableValidation() { return false; }
static std::string getDefaultInputFormat() { return "bam"; }

struct AlPairCmp
{
	bool operator()(libmaus::bambam::BamAlignment const * A, libmaus::bambam::BamAlignment const * B)
	{
		bool const apaired = A->isPaired();
		bool const bpaired = B->isPaired();
		
		// A smaller when in pair but B is not
		if ( apaired != bpaired )
			return apaired;
			
		bool const asec = A->isSecondary();
		bool const bsec = B->isSecondary();

		// A smaller when not secondary
		if ( asec != bsec )
			return !asec;

		bool const asup = A->isSupplementary();
		bool const bsup = B->isSupplementary();

		// A smaller when not supplementary
		if ( asup != bsup )
			return !asup;
			
		bool const afirst = A->isRead1();
		bool const bfirst = B->isRead1();

		// A smaller when first				
		if ( afirst != bfirst )
			return afirst;
		
		return false;
	}
};

void handleAlignmentVector(
	std::vector<libmaus::bambam::BamAlignment *> & Calgn,
	libmaus::bambam::BamBlockWriterBase * singlewr,
	libmaus::bambam::BamBlockWriterBase * orphanwr,
	libmaus::bambam::BamBlockWriterBase * supplementarywr,
	libmaus::bambam::BamBlockWriterBase * unmappedwr,
	libmaus::bambam::BamBlockWriterBase * splitwr,
	libmaus::bambam::BamBlockWriterBase * samestrandwr,
	libmaus::bambam::BamBlockWriterBase * improperwr,
	libmaus::bambam::BamBlockWriterBase * properwr
)	
{
	std::sort(Calgn.begin(),Calgn.end(),AlPairCmp());
	
	#if 0
	std::cerr << Aalgn.size() << std::string(80,'-') << std::endl;
	for ( uint64_t i = 0; i < Calgn.size(); ++i )
		std::cerr << Calgn[i]->formatAlignment(header) << std::endl;
	#endif
		
	for ( uint64_t i = 1; i < Calgn.size(); ++i )
		if ( Calgn[i]->isPaired() != Calgn[0]->isPaired() )
		{
			::libmaus::exception::LibMausException se;
			se.getStream() << "[E] file is broken, read " << Calgn[i]->getName() << " is in a pair and not in a pair" << std::endl;
			se.finish();
			throw se;	
		}
	
	// single end
	if ( ! Calgn[0]->isPaired() )
	{
		for ( uint64_t i = 0; i < Calgn.size(); ++i )
			singlewr->writeAlignment(*(Calgn[i]));
	}
	// paired
	else
	{
		// first and second vectors
		std::vector< libmaus::bambam::BamAlignment * > R1;
		std::vector< libmaus::bambam::BamAlignment * > R2;
		for ( uint64_t i = 0; i < Calgn.size(); ++i )
		{
			int const isr1 = Calgn[i]->isRead1();
			int const isr2 = Calgn[i]->isRead2();
			if ( isr1 + isr2 != 1 )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[E] cannot handle read " << Calgn[i]->getName() << " which is not either read 1 or read 2" << std::endl;
				se.finish();
				throw se;	
			}
			
			if ( isr1 )
				R1.push_back(Calgn[i]);
			else
				R2.push_back(Calgn[i]);
		}
		// are the reads orphans?
		if ( R1.size() == 0 || R2.size() == 0 )
		{
			for ( uint64_t i = 0; i < Calgn.size(); ++i )
				orphanwr->writeAlignment(*(Calgn[i]));				
		}
		else
		{
			// check that first entries in both lists are primary
			if ( 
				R1[0]->isSupplementary()
				||
				R1[0]->isSecondary()
			)
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[E] first read for name " << R1[0]->getName() << " is not primary" << std::endl;
				se.finish();
				throw se;								
			}
			if ( 
				R2[0]->isSupplementary()
				||
				R2[0]->isSecondary()
			)
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[E] second read for name " << R2[0]->getName() << " is not primary" << std::endl;
				se.finish();
				throw se;								
			}
		
			// extract and remove supplementary reads	
			uint64_t o = 0;
			for ( uint64_t i = 0; i < R1.size(); ++i )
				if ( R1[i]->isSupplementary() )
					supplementarywr->writeAlignment(*(R1[i]));
				else
					R1[o++] = R1[i];
			R1.resize(o);
			o = 0;
			for ( uint64_t i = 0; i < R2.size(); ++i )
				if ( R2[i]->isSupplementary() )
					supplementarywr->writeAlignment(*(R2[i]));
				else
					R2[o++] = R2[i];
			R2.resize(o);
			
			for ( uint64_t i = 1; i < R1.size(); ++i )
				if ( ! R1[i]->isSecondary() )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "[E] multiple primary mappings for read 1 of name " << R1[0]->getName() << std::endl;
					se.finish();
					throw se;											
				}
			
			for ( uint64_t i = 1; i < R2.size(); ++i )
				if ( ! R2[i]->isSecondary() )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "[E] multiple primary mappings for read 2 of name " << R1[0]->getName() << std::endl;
					se.finish();
					throw se;
				}

			for ( uint64_t i = 1; i < R1.size(); ++i )
				if ( R1[i]->isMapped() != R1[0]->isMapped() )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "[E] read 1 for name " << R1[0]->getName() << " is mapped and unmapped" << std::endl;
					se.finish();
					throw se;							
				}
			if ( R1.size() > 1 && (!R1[0]->isMapped()) )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[E] read 1 for name " << R1[0]->getName() << " has multiple unmapped versions" << std::endl;
				se.finish();
				throw se;													
			}
			for ( uint64_t i = 1; i < R1.size(); ++i )
				if ( R2[i]->isMapped() != R2[0]->isMapped() )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "[E] read 2 for name " << R2[0]->getName() << " is mapped and unmapped" << std::endl;
					se.finish();
					throw se;							
				}
			if ( R2.size() > 1 && (!R2[0]->isMapped()) )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[E] read 2 for name " << R1[0]->getName() << " has multiple unmapped versions" << std::endl;
				se.finish();
				throw se;													
			}
			
			#if 0
			bool properfirst = true;
			bool improperfirst = true;
			bool samestrandfirst = true;
			bool splitfirst = true;
			bool unmappedfirst = true;
			#endif
			
			if ( R1.size() == 1 && R2.size() == 1 )
			{
				libmaus::bambam::BamAlignment * r1 = R1[0];
				libmaus::bambam::BamAlignment * r2 = R2[0];
				
				// at least one unmapped
				if ( (!(r1->isMapped())) || (!(r1->isMapped())) )
				{
					unmappedwr->writeAlignment(*r1);
					unmappedwr->writeAlignment(*r2);
					// unmappedfirst = false;
				}
				// not on same reference sequence
				else if ( r1->getRefID() != r2->getRefID() )
				{
					splitwr->writeAlignment(*r1);
					splitwr->writeAlignment(*r2);
					// splitfirst = false;							
				}
				// both on same strand
				else if ( r1->isReverse() == r2->isReverse() )
				{							
					samestrandwr->writeAlignment(*r1);
					samestrandwr->writeAlignment(*r2);
					// samestrandfirst = false;							
				}
				// improper
				else if ( ! r1->isProper() )
				{
					if ( r2->isProper() )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "[E] mate information for name " << r1->getName() << " is not consistent" << std::endl;
						se.finish();
						throw se;
					}
					improperwr->writeAlignment(*r1);
					improperwr->writeAlignment(*r2);
					// improperfirst = false;							
				}
				else
				{
					if ( !r2->isProper() )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "[E] mate information for name " << r1->getName() << " is not consistent" << std::endl;
						se.finish();
						throw se;
					}
					properwr->writeAlignment(*r1);
					properwr->writeAlignment(*r2);
					// properfirst = false;
				}
			}
			else
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[E] secondary mapping are not yet supported (name " << R1[0]->getName() << ")" << std::endl;
				se.finish();
				throw se;
				
				#if 0
				for ( uint64_t i = 0; i < R1.size(); ++i )
					for ( uint64_t j = 0; j < R2.size(); ++j )
					{
						libmaus::bambam::BamAlignment * r1 = R1[i];
						libmaus::bambam::BamAlignment * r2 = R2[j];
					}
				#endif
			}
		}
	}

}

int bamflagsplit(::libmaus::util::ArgInfo const & arginfo)
{
	::libmaus::util::TempFileRemovalContainer::setup();
	
	bool const inputisstdin = (!arginfo.hasArg("I")) || (arginfo.getUnparsedValue("I","-") == "-");

	if ( isatty(STDIN_FILENO) && inputisstdin && (arginfo.getValue<std::string>("inputformat","bam") != "sam") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	std::string const inputformat = arginfo.getUnparsedValue("inputformat",getDefaultInputFormat());
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	#if 0
	int const insertmin = arginfo.getValue<int>("insertmin",0);
	int const insertmax = arginfo.getValue<int>("insertmax",500);
	#endif
	bool const disablevalidation = arginfo.getValue<int>("disablevalidation",getDefaultDisableValidation());

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	if ( disablevalidation )
		dec.disableValidation();
	::libmaus::bambam::BamHeader const & header = dec.getHeader();

	
	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilenameout = tmpfilenamebase + "_bamflagsplit";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenameout);
	
	::libmaus::bambam::BamHeader::unique_ptr_type genuphead(
		libmaus::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bamflagsplit",std::string(PACKAGE_VERSION))
	);
	::libmaus::bambam::BamHeader::unique_ptr_type splituphead(	
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\tsplit reads\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type singleuphead(	
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\tsingle reads\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type orphanuphead(	
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\torphan reads\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type unmappeduphead(	
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\tunmapped reads\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type supplementaryuphead(	
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\tsupplementary reads\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type improperuphead(
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\timproperly mapped reads\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type samestranduphead(
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\treads with both ends mapping to same strand\n")
	);
	::libmaus::bambam::BamHeader::unique_ptr_type properuphead(
		new ::libmaus::bambam::BamHeader(genuphead->text + "@CO\tproperly mapped reads\n")
	);

	if ( ! arginfo.hasArg("split") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for split alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}
	
	std::string const splitfn = arginfo.getUnparsedValue("split","notset");

	if ( ! arginfo.hasArg("single") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for single alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}
	
	std::string const singlefn = arginfo.getUnparsedValue("single","notset");

	if ( ! arginfo.hasArg("orphan") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for orphan alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}
	
	std::string const orphanfn = arginfo.getUnparsedValue("orphan","notset");

	if ( ! arginfo.hasArg("unmapped") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for unmapped alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}
	
	std::string const unmappedfn = arginfo.getUnparsedValue("unmapped","notset");

	if ( ! arginfo.hasArg("supplementary") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for supplementary alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}
	
	std::string const supplementaryfn = arginfo.getUnparsedValue("supplementary","notset");

	if ( ! arginfo.hasArg("improper") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for improper alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}
	
	std::string const improperfn = arginfo.getUnparsedValue("improper","notset");

	if ( ! arginfo.hasArg("samestrand") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for samestrand alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}

	std::string const samestrandfn = arginfo.getUnparsedValue("samestrand","notset");

	if ( ! arginfo.hasArg("proper") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "File name for proper alignments is missing" << std::endl;
		se.finish();
		throw se;		
	}

	std::string const properfn = arginfo.getUnparsedValue("proper","notset");
	
	remove(splitfn.c_str());
	remove(singlefn.c_str());
	remove(orphanfn.c_str());
	remove(unmappedfn.c_str());
	remove(supplementaryfn.c_str());
	remove(samestrandfn.c_str());
	remove(improperfn.c_str());
	remove(properfn.c_str());
	
	libmaus::aio::PosixFdOutputStream::unique_ptr_type splitfile(new libmaus::aio::PosixFdOutputStream(splitfn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type singlefile(new libmaus::aio::PosixFdOutputStream(singlefn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type orphanfile(new libmaus::aio::PosixFdOutputStream(orphanfn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type unmappedfile(new libmaus::aio::PosixFdOutputStream(unmappedfn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type supplementaryfile(new libmaus::aio::PosixFdOutputStream(supplementaryfn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type samestrandfile(new libmaus::aio::PosixFdOutputStream(samestrandfn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type improperfile(new libmaus::aio::PosixFdOutputStream(improperfn));
	libmaus::aio::PosixFdOutputStream::unique_ptr_type properfile(new libmaus::aio::PosixFdOutputStream(properfn));

	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Psplitmd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Psinglemd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Porphanmd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Punmappedmd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Psupplementarymd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pimpropermd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Psamestrandmd5;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Ppropermd5;

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > splitcbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > singlecbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > orphancbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > unmappedcbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > supplementarycbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > samestrandcbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > impropercbs;
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > propercbs;
	
	if ( arginfo.hasArg("splitmd5") && arginfo.hasArg("splitmd5filename") && arginfo.getValue<unsigned int>("splitmd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tsplitmd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Psplitmd5 = UNIQUE_PTR_MOVE(Tsplitmd5);
		splitcbs.push_back(Psplitmd5.get());
	}
	if ( arginfo.hasArg("singlemd5") && arginfo.hasArg("singlemd5filename") && arginfo.getValue<unsigned int>("singlemd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tsinglemd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Psinglemd5 = UNIQUE_PTR_MOVE(Tsinglemd5);
		singlecbs.push_back(Psinglemd5.get());
	}
	if ( arginfo.hasArg("orphanmd5") && arginfo.hasArg("orphanmd5filename") && arginfo.getValue<unsigned int>("orphanmd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Torphanmd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Porphanmd5 = UNIQUE_PTR_MOVE(Torphanmd5);
		orphancbs.push_back(Porphanmd5.get());
	}
	if ( arginfo.hasArg("unmappedmd5") && arginfo.hasArg("unmappedmd5filename") && arginfo.getValue<unsigned int>("unmappedmd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tunmappedmd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Punmappedmd5 = UNIQUE_PTR_MOVE(Tunmappedmd5);
		unmappedcbs.push_back(Punmappedmd5.get());
	}
	if ( arginfo.hasArg("supplementarymd5") && arginfo.hasArg("supplementarymd5filename") && arginfo.getValue<unsigned int>("supplementarymd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tsupplementarymd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Psupplementarymd5 = UNIQUE_PTR_MOVE(Tsupplementarymd5);
		supplementarycbs.push_back(Psupplementarymd5.get());
	}
	if ( arginfo.hasArg("impropermd5") && arginfo.hasArg("impropermd5filename") && arginfo.getValue<unsigned int>("impropermd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Timpropermd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Pimpropermd5 = UNIQUE_PTR_MOVE(Timpropermd5);
		impropercbs.push_back(Pimpropermd5.get());
	}
	if ( arginfo.hasArg("samestrandmd5") && arginfo.hasArg("samestrandmd5filename") && arginfo.getValue<unsigned int>("samestrandmd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tsamestrandmd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Psamestrandmd5 = UNIQUE_PTR_MOVE(Tsamestrandmd5);
		samestrandcbs.push_back(Psamestrandmd5.get());
	}
	if ( arginfo.hasArg("propermd5") && arginfo.hasArg("propermd5filename") && arginfo.getValue<unsigned int>("propermd5",0) )
	{
		::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tpropermd5(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
		Ppropermd5 = UNIQUE_PTR_MOVE(Tpropermd5);
		propercbs.push_back(Ppropermd5.get());
	}
	
	#if 0
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	#endif


	libmaus::bambam::BamBlockWriterBase::unique_ptr_type splitwr(
		new libmaus::bambam::BamWriter(*splitfile,*splituphead,level,splitcbs.size() ? (&splitcbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type singlewr(
		new libmaus::bambam::BamWriter(*singlefile,*singleuphead,level,singlecbs.size() ? (&singlecbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type orphanwr(
		new libmaus::bambam::BamWriter(*orphanfile,*orphanuphead,level,orphancbs.size() ? (&orphancbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type unmappedwr(
		new libmaus::bambam::BamWriter(*unmappedfile,*unmappeduphead,level,unmappedcbs.size() ? (&unmappedcbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type supplementarywr(
		new libmaus::bambam::BamWriter(*supplementaryfile,*supplementaryuphead,level,supplementarycbs.size() ? (&supplementarycbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type improperwr(
		new libmaus::bambam::BamWriter(*improperfile,*improperuphead,level,impropercbs.size() ? (&impropercbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type samestrandwr(
		new libmaus::bambam::BamWriter(*samestrandfile,*samestranduphead,level,samestrandcbs.size() ? (&samestrandcbs) : 0)
	);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type properwr(
		new libmaus::bambam::BamWriter(*properfile,*properuphead,level,propercbs.size() ? (&propercbs) : 0)
	);
	
	libmaus::bambam::BamAlignment & curalgn = dec.getAlignment();
	libmaus::bambam::BamAuxFilterVector MQfilter;
	MQfilter.set("MQ");
	uint64_t c = 0;
	
	std::vector < libmaus::bambam::BamAlignment::shared_ptr_type > Aalgn;
	std::vector < libmaus::bambam::BamAlignment * > Falgn;
	std::vector < libmaus::bambam::BamAlignment * > Calgn;

	while ( dec.readAlignment() )
	{
		// no more free alignment blocks?
		if ( ! Falgn.size() )
		{
			// allocate new block
			libmaus::bambam::BamAlignment::shared_ptr_type Tptr(new libmaus::bambam::BamAlignment);
			Aalgn.push_back(Tptr);
			Falgn.push_back(Aalgn.back().get());
		}

		libmaus::bambam::BamAlignment * const calgn = Falgn.back();
		Falgn.pop_back();
		calgn->swap(curalgn);
	
		// new name?
		if ( !Calgn.size() || (strcmp(Calgn.back()->getName(),calgn->getName()) != 0) )
		{
			if ( Calgn.size() )
			{
				handleAlignmentVector(Calgn,singlewr.get(),orphanwr.get(),supplementarywr.get(),
					unmappedwr.get(),splitwr.get(),samestrandwr.get(),
					improperwr.get(),properwr.get());

				// move alignments to free list				
				while ( Calgn.size() )
				{
					Falgn.push_back(Calgn.back());
					Calgn.pop_back();
				}
			}
		}

		Calgn.push_back(calgn);
		
		if ( verbose && ( ( ++c & ((1ull<<20)-1) ) == 0 ) )
			std::cerr << "[V] " << c << std::endl;
	}

	if ( Calgn.size() )
	{
		handleAlignmentVector(Calgn,singlewr.get(),orphanwr.get(),supplementarywr.get(),
			unmappedwr.get(),splitwr.get(),samestrandwr.get(),
			improperwr.get(),properwr.get());
	}
	
	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

	splitwr.reset();
	singlewr.reset();
	orphanwr.reset();
	unmappedwr.reset();
	supplementarywr.reset();
	improperwr.reset();
	samestrandwr.reset();
	properwr.reset();

	splitfile->flush(); splitfile.reset();
	singlefile->flush(); singlefile.reset();
	orphanfile->flush(); orphanfile.reset();
	unmappedfile->flush(); unmappedfile.reset();
	supplementaryfile->flush(); supplementaryfile.reset();
	improperfile->flush(); improperfile.reset();
	samestrandfile->flush(); samestrandfile.reset();
	properfile->flush(); properfile.reset();

	if ( Psplitmd5 )
		Psplitmd5->saveDigestAsFile(arginfo.getUnparsedValue("splitmd5filename","not set"));
	if ( Psinglemd5 )
		Psinglemd5->saveDigestAsFile(arginfo.getUnparsedValue("singlemd5filename","not set"));
	if ( Porphanmd5 )
		Porphanmd5->saveDigestAsFile(arginfo.getUnparsedValue("orphanmd5filename","not set"));
	if ( Punmappedmd5 )
		Punmappedmd5->saveDigestAsFile(arginfo.getUnparsedValue("unmappedmd5filename","not set"));
	if ( Psupplementarymd5 )
		Psupplementarymd5->saveDigestAsFile(arginfo.getUnparsedValue("supplementarymd5filename","not set"));
	if ( Pimpropermd5 )
		Pimpropermd5->saveDigestAsFile(arginfo.getUnparsedValue("impropermd5filename","not set"));
	if ( Psamestrandmd5 )
		Psamestrandmd5->saveDigestAsFile(arginfo.getUnparsedValue("samestrandmd5filename","not set"));
	if ( Ppropermd5 )
		Ppropermd5->saveDigestAsFile(arginfo.getUnparsedValue("propermd5filename","not set"));
	
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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );				
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				
				V.push_back ( std::pair<std::string,std::string> ( "single=<filename>", "output file name for single file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "singlemd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for single file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "singlemd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "unmapped=<filename>", "output file name for unmapped file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "unmappedmd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for unmapped file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "unmappedmd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "orphan=<filename>", "output file name for orphan file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "orphanmd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for orphan file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "orphanmd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "split=<filename>", "output file name for split file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "splitmd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for split file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "splitmd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "supplementary=<filename>", "output file name for supplementary file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "supplementarymd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for supplementary file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "supplementarymd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "samestrand=<filename>", "output file name for samestrand file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "samestrandmd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for samestrand file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "samestrandmd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "improper=<filename>", "output file name for improper file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "impropermd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for improper file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "impropermd5filename=<filename>", "file name for md5 check sum of sigle file" ) );

				V.push_back ( std::pair<std::string,std::string> ( "proper=<filename>", "output file name for proper file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "propermd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for proper file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "propermd5filename=<filename>", "file name for md5 check sum of sigle file" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamflagsplit(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
