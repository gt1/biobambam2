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
#include <libmaus/bambam/BamToFastqOutputFileSet.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/bambam/BamAuxFilterVector.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/aio/PosixFdInputStream.hpp>

static int getDefaultVerbose() { return 1; }
static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static std::string getDefaultInputFormat() { return "bam"; }
static uint64_t getDefaultInputBufferSize() { return 64*1024; }

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }
static int getDefaultMapQThreshold() { return -1; }
static std::string getDefaultClassFilter() { return "F,F2,O,O2,S"; }
static bool getDefaultResetAux() { return true; }

static uint32_t const classmask_F  = (1ull << 0);
static uint32_t const classmask_F2 = (1ull << 1);
static uint32_t const classmask_O  = (1ull << 2);
static uint32_t const classmask_O2 = (1ull << 3);
static uint32_t const classmask_S  = (1ull << 4);
static uint32_t const classmask_all = classmask_F | classmask_F2 | classmask_O | classmask_O2 | classmask_S;

static std::string despace(std::string const & s)
{
	std::deque<char> Q(s.begin(),s.end());
	while ( Q.size() && isspace(Q.front()) )
		Q.pop_front();
	while ( Q.size() && isspace(Q.back()) )
		Q.pop_back();
	return s;
}

static uint32_t parseClassList(std::string s)
{
	s = despace(s);
	
	if ( ! s.size() )
		return 0;
	
	std::deque<std::string> tokens = libmaus::util::stringFunctions::tokenize<std::string>(s,std::string(","));
	
	uint32_t mask = 0;
	
	for ( uint64_t i = 0; i < tokens.size(); ++i )
	{
		tokens[i] = despace(tokens[i]);
		
		if ( tokens[i] == "F" )
			mask |= classmask_F;
		else if ( tokens[i] == "F2" )
			mask |= classmask_F2;
		else if ( tokens[i] == "O" )
			mask |= classmask_O;
		else if ( tokens[i] == "O2" )
			mask |= classmask_O2;
		else if ( tokens[i] == "S" )
			mask |= classmask_S;
		else
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Unknown read class " << tokens[i] << std::endl;
			se.finish();
			throw se;
		}
	}
	
	return mask;
}

struct BamToFastQInputFileStream
{
	std::string const fn;
	libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	std::istream & in;
	
	static libmaus::aio::CheckedInputStream::unique_ptr_type openFile(std::string const & fn)
	{
		libmaus::aio::CheckedInputStream::unique_ptr_type ptr(new libmaus::aio::CheckedInputStream(fn));
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	BamToFastQInputFileStream(libmaus::util::ArgInfo const & arginfo)
	: fn(arginfo.getValue<std::string>("filename","-")),
	  CIS(
		(fn != "-") ? openFile(fn) : (libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}

	BamToFastQInputFileStream(std::string const & rfn)
	: fn(rfn), CIS(
		(fn != "-") ? openFile(fn) : (libmaus::aio::CheckedInputStream::unique_ptr_type())
	), in((fn != "-") ? (*CIS) : std::cin) {}
};

template<typename decoder_type>
std::string getModifiedHeaderText(decoder_type const & bamdec, libmaus::util::ArgInfo const & arginfo, bool reset = false)
{
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	std::string headertext = header.text;

	// reset header if requested
	if ( reset )
	{
		// no replacement header file given
		if ( ! arginfo.hasArg("resetheadertext") )
		{
			// remove SQ lines
			std::vector<libmaus::bambam::HeaderLine> allheaderlines = libmaus::bambam::HeaderLine::extractLines(headertext);

			std::ostringstream upheadstr;
			for ( uint64_t i = 0; i < allheaderlines.size(); ++i )
				if ( allheaderlines[i].type != "SQ" )
					upheadstr << allheaderlines[i].line << std::endl;

			headertext = upheadstr.str();
		}
		// replace header given in file
		else
		{
			std::string const headerfilename = arginfo.getUnparsedValue("resetheadertext","");
			uint64_t const headerlen = libmaus::util::GetFileSize::getFileSize(headerfilename);
			libmaus::aio::CheckedInputStream CIS(headerfilename);
			libmaus::autoarray::AutoArray<char> ctext(headerlen,false);
			CIS.read(ctext.begin(),headerlen);
			headertext = std::string(ctext.begin(),ctext.end());		
		}
	}

	// add PG line to header
	headertext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamcollate2", // ID
		"bamcollate2", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);

	return headertext;
}

::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type getRGTrie(libmaus::util::ArgInfo const & arginfo)
{
	std::vector < std::string > vreadgroups;
	std::string const readgroups = arginfo.getValue<std::string>("readgroups",std::string());
	::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type LHTsnofailure;
	
	if ( readgroups.size() )
	{
		std::deque<std::string> qreadgroups = ::libmaus::util::stringFunctions::tokenize(readgroups,std::string(","));
		vreadgroups = std::vector<std::string>(qreadgroups.begin(),qreadgroups.end());
		::libmaus::trie::Trie<char> trienofailure;
		trienofailure.insertContainer(vreadgroups);
		::libmaus::trie::LinearHashTrie<char,uint32_t>::unique_ptr_type LHTnofailure(trienofailure.toLinearHashTrie<uint32_t>());
		LHTsnofailure = ::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type(LHTnofailure.release());
	}
	
	return LHTsnofailure;
}

void bamcollate2NonCollating(libmaus::util::ArgInfo const & arginfo, libmaus::bambam::BamAlignmentDecoder & bamdec)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		bamdec.disableValidation();

	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	libmaus::timing::RealTimeClock rtc; rtc.start();
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	uint64_t cnt = 0;
	unsigned int const verbshift = 20;
	bool const reset = arginfo.getValue<unsigned int>("reset",false);
	bool const resetaux = arginfo.getValue<unsigned int>("resetaux",getDefaultResetAux());
	libmaus::bambam::BamAuxFilterVector::unique_ptr_type const prgfilter(libmaus::bambam::BamAuxFilterVector::parseAuxFilterList(arginfo));
	libmaus::bambam::BamAuxFilterVector const * rgfilter = prgfilter.get();

	// construct new header
	::libmaus::bambam::BamHeader uphead(getModifiedHeaderText(bamdec,arginfo,reset));

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
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead,arginfo,Pcbs)
	);

	// read group trie
	::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type const rgtrie = getRGTrie(arginfo);
		
	while ( bamdec.readAlignment() )
	{
		uint64_t const precnt = cnt++;
		
		if ( 
			(! (algn.getFlags() & excludeflags)) 
			&&
			( (!rgtrie) || (algn.getReadGroup() && rgtrie->searchCompleteNoFailureZ(algn.getReadGroup()) != -1) )
		)
		{
			if ( reset )
				resetAlignment(algn,
					resetaux /* reset aux */,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
				        rgfilter /* rg filter */
				);
		
			writer->writeAlignment(algn);
		}

		if ( (precnt >> verbshift != cnt >> verbshift) && verbose )
			std::cerr 
				<< (cnt >> 20) 
				<< "\t"
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() << std::endl;
	}

	if ( verbose )
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

void bamcollate2NonCollating(libmaus::util::ArgInfo const & arginfo)
{
	uint64_t const inputbuffersize = arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",getDefaultInputBufferSize());
	libmaus::aio::PosixFdInputStream PFIS(STDIN_FILENO,inputbuffersize);
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false /* put rank */, 0 /* copy stream */, PFIS
		)
	);
	bamcollate2NonCollating(arginfo,decwrapper->getDecoder());

	std::cout.flush();
}

void bamcollate2Collating(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry * ob = 0;
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	bool const reset = arginfo.getValue<unsigned int>("reset",false);
	bool const resetaux = arginfo.getValue<unsigned int>("resetaux",getDefaultResetAux());
	libmaus::bambam::BamAuxFilterVector::unique_ptr_type const prgfilter(libmaus::bambam::BamAuxFilterVector::parseAuxFilterList(arginfo));
	libmaus::bambam::BamAuxFilterVector const * rgfilter = prgfilter.get();
	libmaus::bambam::BamAlignment Ralgna, Ralgnb;
	std::string const sclassfilter = arginfo.getValue<std::string>("classes",getDefaultClassFilter());
	uint32_t const classmask = parseClassList(sclassfilter);

	// construct new header
	::libmaus::bambam::BamHeader uphead(getModifiedHeaderText(CHCBD,arginfo,reset));
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
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead,arginfo,Pcbs)
	);

	// read group trie
	::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type rgtrie = getRGTrie(arginfo);
	int const mapqthres = arginfo.getValue<int>("mapqthres",getDefaultMapQThreshold());
	
	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair && ((classmask & classmask_F) || (classmask & classmask_F2)) )
		{
			bool pass = true;
			
			if ( rgtrie )
			{
				char const * rga = libmaus::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Da,ob->blocksizea);
				char const * rgb = libmaus::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Db,ob->blocksizeb);
				
				pass = 
					rga && rgb &&
					(rgtrie->searchCompleteNoFailureZ(rga) != -1) &&
					(rgtrie->searchCompleteNoFailureZ(rgb) != -1);
			}
			
			if ( mapqthres >= 0 )
			{
				bool const amapped = !libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Da));
				bool const aok = amapped && static_cast<int>(libmaus::bambam::BamAlignmentDecoderBase::getMapQ(ob->Da)) >= mapqthres;

				bool const bmapped = !libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Db));
				bool const bok = bmapped && static_cast<int>(libmaus::bambam::BamAlignmentDecoderBase::getMapQ(ob->Db)) >= mapqthres;
				
				pass = pass && (aok || bok);
			}
			
			if ( pass )
			{
				if ( reset )
				{
					if ( classmask & classmask_F )
					{
						Ralgna.copyFrom(ob->Da,ob->blocksizea);
						resetAlignment(
							Ralgna,
							resetaux /* reset aux */,
							libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
						        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
						        rgfilter /* RG filter */
							);
						writer->writeAlignment(Ralgna);
						bcnt += (Ralgna.blocksize);
					}
					
					if ( classmask & classmask_F2 )
					{
						Ralgnb.copyFrom(ob->Db,ob->blocksizeb);
						resetAlignment(
							Ralgnb,
							resetaux /* reset aux */,
							libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
						        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
						        rgfilter /* RG filter */
						);
						writer->writeAlignment(Ralgnb);
						bcnt += (Ralgnb.blocksize);
					}
				}
				else
				{
					if ( classmask & classmask_F )
					{
						writer->writeBamBlock(ob->Da,ob->blocksizea);
						bcnt += (ob->blocksizea);
					}
					
					if ( classmask & classmask_F2 )
					{
						writer->writeBamBlock(ob->Db,ob->blocksizeb);
						bcnt += (ob->blocksizeb);
					}
				}

				cnt += 2;
			}
		}
		else if ( ob->fsingle && (classmask & classmask_S) )
		{
			bool pass = true;
			
			if ( rgtrie )
			{
				char const * rga = libmaus::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Da,ob->blocksizea);
				pass = rga && (rgtrie->searchCompleteNoFailureZ(rga) != -1);
			}

			if ( mapqthres >= 0 )
			{
				bool const amapped = !libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Da));
				bool const aok = amapped && static_cast<int>(libmaus::bambam::BamAlignmentDecoderBase::getMapQ(ob->Da)) >= mapqthres;
				
				pass = pass && aok;
			}
			
			if ( pass )
			{
				if ( reset )
				{				
					Ralgna.copyFrom(ob->Da,ob->blocksizea);
					resetAlignment(
						Ralgna,
						resetaux /* reset aux */,
						libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
					        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
					        rgfilter /* RG filter */
					);
					writer->writeAlignment(Ralgna);

					bcnt += (Ralgna.blocksize);
				}
				else
				{
					writer->writeBamBlock(ob->Da,ob->blocksizea);
					bcnt += (ob->blocksizea);
				}

				cnt += 1;
			}
		}
		else if ( ob->forphan1 && (classmask & classmask_O) )
		{
			bool pass = true;
			
			if ( rgtrie )
			{
				char const * rga = libmaus::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Da,ob->blocksizea);
				pass = rga && (rgtrie->searchCompleteNoFailureZ(rga) != -1);
			}

			if ( mapqthres >= 0 )
			{
				bool const amapped = !libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Da));
				bool const aok = amapped && static_cast<int>(libmaus::bambam::BamAlignmentDecoderBase::getMapQ(ob->Da)) >= mapqthres;
				
				pass = pass && aok;
			}
			
			if ( pass )
			{
				if ( reset )
				{				
					Ralgna.copyFrom(ob->Da,ob->blocksizea);
					resetAlignment(
						Ralgna,
						resetaux /* reset aux */,
						libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
					        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
					        rgfilter /* RG filter */					);
					writer->writeAlignment(Ralgna);
					bcnt += (Ralgna.blocksize);
				}
				else
				{
					writer->writeBamBlock(ob->Da,ob->blocksizea);
					bcnt += (ob->blocksizea);
				}

				cnt += 1;
			}
		}
		else if ( ob->forphan2 && (classmask & classmask_O2) )
		{
			bool pass = true;
			
			if ( rgtrie )
			{
				char const * rgb = libmaus::bambam::BamAlignmentDecoderBase::getReadGroup(ob->Db,ob->blocksizeb);
				pass = rgb && (rgtrie->searchCompleteNoFailureZ(rgb) != -1);
			}

			if ( mapqthres >= 0 )
			{
				bool const bmapped = !libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(ob->Db));
				bool const bok = bmapped && static_cast<int>(libmaus::bambam::BamAlignmentDecoderBase::getMapQ(ob->Db)) >= mapqthres;
				
				pass = pass && bok;
			}
			
			if ( pass )
			{
				if ( reset )
				{				
					Ralgna.copyFrom(ob->Da,ob->blocksizea);
					resetAlignment(
						Ralgna,
						resetaux /* reset aux */,
						libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
					        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
					        rgfilter /* RG filter */
					);
					writer->writeAlignment(Ralgna);
					bcnt += (Ralgna.blocksize);
				}
				else
				{
					writer->writeBamBlock(ob->Da,ob->blocksizea);
					bcnt += (ob->blocksizea);
				}

				cnt += 1;
			}
		}
		
		if ( (precnt >> verbshift != cnt >> verbshift) && verbose )
		{
			std::cerr 
				<< "[V] "
				<< (cnt >> 20) 
				<< "\t"
				<< (static_cast<double>(bcnt)/(1024.0*1024.0))/rtc.getElapsedSeconds() << "MB/s"
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() << std::endl;
		}
	}
	
	if ( verbose )
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

void bamcollate2Collating(libmaus::util::ArgInfo const & arginfo)
{
	uint32_t const excludeflags = libmaus::bambam::BamFlagBase::stringToFlags(arginfo.getValue<std::string>("exclude","SECONDARY,SUPPLEMENTARY"));
	libmaus::util::TempFileRemovalContainer::setup();
	std::string const tmpfilename = arginfo.getValue<std::string>("T",arginfo.getDefaultTmpFileName());
	libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	
	unsigned int const hlog = arginfo.getValue<unsigned int>("colhlog",18);
	uint64_t const sbs = arginfo.getValueUnsignedNumeric<uint64_t>("colsbs",128ull*1024ull*1024ull);

	uint64_t const inputbuffersize = arginfo.getValueUnsignedNumeric<uint64_t>("inputbuffersize",getDefaultInputBufferSize());
	libmaus::aio::PosixFdInputStream PFIS(STDIN_FILENO,inputbuffersize);
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false /* put rank */, 0 /* copy stream */, PFIS
		)
	);
	libmaus::bambam::CircularHashCollatingBamDecoder CHCBD(decwrapper->getDecoder(),tmpfilename,excludeflags,hlog,sbs);
	bamcollate2Collating(arginfo,CHCBD);
	
	std::cout.flush();
}

void bamcollate2CollatingRanking(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	// number of bytes written to files
	uint64_t bcnt = 0;
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	unsigned int const verbshift = 20;
	libmaus::timing::RealTimeClock rtc; rtc.start();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	libmaus::bambam::BamAlignment algn;

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
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead,arginfo,Pcbs)
	);

	libmaus::bambam::BamAuxFilterVector zrtag;
	zrtag.set('Z','R');

	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			uint64_t const rankb = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Db,ob->blocksizeb);
							
			std::ostringstream nameostr;
			nameostr << ranka << "_" << rankb << "_" << libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			std::string const name = nameostr.str();
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(name.c_str(),name.size());
			algn.filterOutAux(zrtag);
			writer->writeAlignment(algn);
			
			if ( algn.D.size() < ob->blocksizeb )
				algn.D.resize(ob->blocksizeb);
			std::copy ( ob->Db, ob->Db + ob->blocksizeb, algn.D.begin() );
			algn.blocksize = ob->blocksizeb;
			algn.replaceName(name.c_str(),name.size());
			algn.filterOutAux(zrtag);
			writer->writeAlignment(algn);
			
			cnt += 2;
		}
		else if ( ob->fsingle )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);

				
			std::ostringstream nameostr;
			nameostr << ranka << "_" << libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			std::string const name = nameostr.str();
				
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(name.c_str(),name.size());
			algn.filterOutAux(zrtag);
			writer->writeAlignment(algn);

			cnt += 1;
		}
		else if ( ob->forphan1 )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);

			std::ostringstream nameostr;
			nameostr << ranka << "_" << ranka << "_" << libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			std::string const name = nameostr.str();
				
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(name.c_str(),name.size());
			algn.filterOutAux(zrtag);
			writer->writeAlignment(algn);

			cnt += 1;
		}
		else if ( ob->forphan2 )
		{
			uint64_t const ranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);

			std::ostringstream nameostr;
			nameostr << ranka << "_" << ranka << "_" << libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);
			std::string const name = nameostr.str();
				
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(name.c_str(),name.size());
			algn.filterOutAux(zrtag);
			writer->writeAlignment(algn);

			cnt += 1;
		}
		
		if ( (precnt >> verbshift != cnt >> verbshift) && verbose )
		{
			std::cerr 
				<< "[V] "
				<< (cnt >> 20) 
				<< "\t"
				<< (static_cast<double>(bcnt)/(1024.0*1024.0))/rtc.getElapsedSeconds() << "MB/s"
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() << std::endl;
		}
	}
	
	if ( verbose )
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

void bamcollate2CollatingRanking(libmaus::util::ArgInfo const & arginfo)
{
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
		BamToFastQInputFileStream bamin(inputfilename);

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
			bamcollate2CollatingRanking(arginfo,CHCBD);
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
			bamcollate2CollatingRanking(arginfo,CHCBD);
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
		bamcollate2CollatingRanking(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			true, /* put rank */
			hlog,sbs
		);
		bamcollate2CollatingRanking(arginfo,CHCBD);
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

void bamcollate2CollatingPostRanking(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bambam::CircularHashCollatingBamDecoder & CHCBD
)
{
	int const reset = arginfo.getValue<int>("reset",1);
	bool const resetaux = arginfo.getValue<int>("resetaux",0);
	libmaus::bambam::BamAuxFilterVector::unique_ptr_type const prgfilter(libmaus::bambam::BamAuxFilterVector::parseAuxFilterList(arginfo));
	libmaus::bambam::BamAuxFilterVector const * rgfilter = prgfilter.get();

	if ( arginfo.getValue<unsigned int>("disablevalidation",0) )
		CHCBD.disableValidation();

	libmaus::bambam::CircularHashCollatingBamDecoder::OutputBufferEntry const * ob = 0;
	
	// number of alignments written to files
	uint64_t cnt = 0;
	unsigned int const verbshift = 20;
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	libmaus::timing::RealTimeClock rtc; rtc.start();
	::libmaus::autoarray::AutoArray<uint8_t> T;
	libmaus::bambam::BamAlignment algn;

	// construct new header
	::libmaus::bambam::BamHeader uphead(getModifiedHeaderText(CHCBD,arginfo,reset));
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
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead,arginfo,Pcbs)
	);
	uint64_t r = 0;

	libmaus::bambam::BamAuxFilterVector zrtag;
	zrtag.set('Z','R');

	libmaus::bambam::BamAuxFilterVector zzbafv;
	zzbafv.set('z','z');

	libmaus::autoarray::AutoArray<char> namebuffer;

	while ( (ob = CHCBD.process()) )
	{
		uint64_t const precnt = cnt;
		
		if ( ob->fpair )
		{
			uint64_t const ranka = r++;
			uint64_t const rankb = r++;
			uint64_t const zranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			uint64_t const zrankb = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Db,ob->blocksizeb);

			uint64_t const orignamelen = (::libmaus::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1);
			char const * origname = libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);			
			uint64_t const namelen = 
				libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(ranka) + 1 +
				libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(rankb) + 1 +
				orignamelen
				;

			if ( namelen > namebuffer.size() )
				namebuffer = libmaus::autoarray::AutoArray<char>(namelen);
				
			char * np = namebuffer.begin();
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,ranka);
			*(np++) = '_';
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,rankb);
			*(np++) = '_';
			std::copy(origname,origname + orignamelen,np);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(namebuffer.begin(),namelen);
			algn.filterOutAux(zrtag);
			if ( reset )
				resetAlignment(
					algn,
					resetaux,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
				        rgfilter /* RG filter */
				);
			attachRank(algn,zranka,zzbafv);
			writer->writeAlignment(algn);
			
			if ( algn.D.size() < ob->blocksizeb )
				algn.D.resize(ob->blocksizeb);
			std::copy ( ob->Db, ob->Db + ob->blocksizeb, algn.D.begin() );
			algn.blocksize = ob->blocksizeb;
			algn.replaceName(namebuffer.begin(),namelen);
			algn.filterOutAux(zrtag);
			if ( reset )
				resetAlignment(
					algn,
					resetaux,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
				        rgfilter /* RG filter */
				);
			attachRank(algn,zrankb,zzbafv);
			writer->writeAlignment(algn);
			
			cnt += 2;
		}
		else if ( ob->fsingle )
		{
			uint64_t const ranka = r++;
			uint64_t const zranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);
				
			uint64_t const orignamelen = (::libmaus::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1);
			char const * origname = libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);	
			uint64_t const namelen = libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(ranka) + 1 + orignamelen;

			if ( namelen > namebuffer.size() )
				namebuffer = libmaus::autoarray::AutoArray<char>(namelen);
				
			char * np = namebuffer.begin();
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,ranka);
			*(np++) = '_';
			std::copy(origname,origname + orignamelen,np);
	
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(namebuffer.begin(),namelen);
			algn.filterOutAux(zrtag);
			if ( reset )
				resetAlignment(
					algn,
					resetaux,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
				        rgfilter /* RG filter */
				);
			attachRank(algn,zranka,zzbafv);
			writer->writeAlignment(algn);

			cnt += 1;
		}
		else if ( ob->forphan1 )
		{
			uint64_t const ranka = r++;
			uint64_t const zranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);

			uint64_t const orignamelen = (::libmaus::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1);
			char const * origname = libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);			
			uint64_t const namelen = 
				libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(ranka) + 1 +
				libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(ranka) + 1 +
				orignamelen
				;

			if ( namelen > namebuffer.size() )
				namebuffer = libmaus::autoarray::AutoArray<char>(namelen);
				
			char * np = namebuffer.begin();
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,ranka);
			*(np++) = '_';
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,ranka);
			*(np++) = '_';
			std::copy(origname,origname + orignamelen,np);

			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(namebuffer.begin(),namelen);
			algn.filterOutAux(zrtag);
			if ( reset )
				resetAlignment(
					algn,
					resetaux,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
				        rgfilter /* RG filter */				
				);
			attachRank(algn,zranka,zzbafv);
			writer->writeAlignment(algn);

			cnt += 1;
		}
		else if ( ob->forphan2 )
		{
			uint64_t const ranka = r++;
			uint64_t const zranka = libmaus::bambam::BamAlignmentDecoderBase::getRank(ob->Da,ob->blocksizea);
			
			if ( algn.D.size() < ob->blocksizea )
				algn.D.resize(ob->blocksizea);

			uint64_t const orignamelen = (::libmaus::bambam::BamAlignmentDecoderBase::getLReadName(ob->Da)-1);
			char const * origname = libmaus::bambam::BamAlignmentDecoderBase::getReadName(ob->Da);			
			uint64_t const namelen = 
				libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(ranka) + 1 +
				libmaus::bambam::BamAlignmentDecoderBase::getDecimalNumberLength(ranka) + 1 +
				orignamelen
				;

			if ( namelen > namebuffer.size() )
				namebuffer = libmaus::autoarray::AutoArray<char>(namelen);
				
			char * np = namebuffer.begin();
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,ranka);
			*(np++) = '_';
			np = libmaus::bambam::BamAlignmentDecoderBase::putNumberDecimal(np,ranka);
			*(np++) = '_';
			std::copy(origname,origname + orignamelen,np);
		
			std::copy ( ob->Da, ob->Da + ob->blocksizea, algn.D.begin() );
			algn.blocksize = ob->blocksizea;
			algn.replaceName(namebuffer.begin(),namelen);
			algn.filterOutAux(zrtag);
			if ( reset )
				resetAlignment(
					algn,
					resetaux,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
				        libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY,
				        rgfilter /* RG filter */
				);
			attachRank(algn,zranka,zzbafv);
			writer->writeAlignment(algn);

			cnt += 1;
		}
		
		if ( (precnt >> verbshift != cnt >> verbshift) && verbose )
		{
			std::cerr 
				<< "[V] "
				<< (cnt >> 20) 
				<< "\t" << static_cast<double>(cnt)/rtc.getElapsedSeconds() << std::endl;
		}
	}
	
	if ( verbose )
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

void bamcollate2CollatingPostRanking(libmaus::util::ArgInfo const & arginfo)
{
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
		BamToFastQInputFileStream bamin(inputfilename);

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
			bamcollate2CollatingPostRanking(arginfo,CHCBD);
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
			bamcollate2CollatingPostRanking(arginfo,CHCBD);
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
		bamcollate2CollatingPostRanking(arginfo,CHCBD);
	}
	else if ( inputformat == "cram" )
	{
		std::string const reference = arginfo.getValue<std::string>("reference","");
		libmaus::bambam::ScramCircularHashCollatingBamDecoder CHCBD(inputfilename,"rc",reference,
			tmpfilename,excludeflags,
			true, /* put rank */
			hlog,sbs
		);
		bamcollate2CollatingPostRanking(arginfo,CHCBD);
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

void bamcollate2(libmaus::util::ArgInfo const & arginfo)
{
	if ( 
		arginfo.hasArg("ranges") && arginfo.getValue("inputformat", getDefaultInputFormat()) != "bam" 
		&&
		arginfo.hasArg("ranges") && arginfo.getValue("inputformat", getDefaultInputFormat()) != "cram" 
	)
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
			bamcollate2NonCollating(arginfo);
			break;
		case 1:
			bamcollate2Collating(arginfo);
			break;
		case 2:		
			bamcollate2CollatingRanking(arginfo);
			break;
		case 3:	
			bamcollate2CollatingPostRanking(arginfo);
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
		
		::libmaus::util::ArgInfo arginfo(argc,argv);
	
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
				
				V.push_back ( std::pair<std::string,std::string> ( "collate=<[1]>", "collate pairs" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reset=<>", "reset alignments and header like bamreset (for collate=0,1 or 3 only, default enabled for 3)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "filename=<[stdin]>", "input filename (default: read file from standard input)" ) );
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( "ranges=<[]>", "input ranges (bam/cram input only, collate<2 only, default: read complete file)" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "ranges=<[]>", "input ranges (bam input only, collate<2 only, default: read complete file)" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "exclude=<[SECONDARY,SUPPLEMENTARY]>", "exclude alignments matching any of the given flags" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<[0]>", "disable validation of input data" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhlog=<[18]>", "base 2 logarithm of hash table size used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("colsbs=<[")+libmaus::util::NumberSerialisation::formatNumber(128ull*1024*1024,0)+"]>", "size of hash table overflow list in bytes" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("T=<[") + arginfo.getDefaultTmpFileName() + "]>" , "temporary file name" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "readgroups=[<>]", "read group filter (default: keep all)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("mapqthres=<[")+::biobambam::Licensing::formatNumber(getDefaultMapQThreshold())+"]>", "mapping quality threshold (collate=1 only, default: keep all)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("classes=[") + getDefaultClassFilter() + std::string("]"), "class filter (collate=1 only, default: keep all)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "resetheadertext=[<>]", "replacement SAM header text file for reset=1 (default: filter header in source BAM file)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("resetaux=<[")+::biobambam::Licensing::formatNumber(getDefaultResetAux())+"]>", "reset auxiliary fields (collate=0,1 only with reset=1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "auxfilter=[<>]", "comma separated list of aux tags to keep if reset=1 and resetaux=0 (default: keep all)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );				
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputbuffersize=<[")+::biobambam::Licensing::formatNumber(getDefaultInputBufferSize())+"]>", "input buffer size" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

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
		if ( arginfo.hasArg("threads") )
		{
			std::string const threads = arginfo.getUnparsedValue("threads",std::string());		
			arginfo.argmap["inputthreads"] = threads;
			arginfo.argmultimap.insert(std::pair<std::string,std::string>(std::string("inputthreads"),threads));
		}
					
		bamcollate2(arginfo);
		
		if ( arginfo.getValue<unsigned int>("verbose",getDefaultVerbose()) )
			std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
