/**
    bambam
    Copyright (C) 2009-2016 German Tischler
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
#include <config.h>

#include <libmaus2/bambam/OptNameReader.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/aio/OutputStreamInstance.hpp>
#include <libmaus2/bambam/BamAlignmentInputCallbackBam.hpp>
#include <libmaus2/bambam/BamAlignmentInputCallbackSnappy.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamHeaderUpdate.hpp>
#include <libmaus2/bambam/BamParallelRewrite.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/BamAlignmentSnappyInput.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus2/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus2/bambam/CollatingBamDecoder.hpp>
#include <libmaus2/bambam/DuplicationMetrics.hpp>
#include <libmaus2/bambam/DupMarkBase.hpp>
#include <libmaus2/bambam/DupSetCallbackVector.hpp>
#include <libmaus2/bambam/OpticalComparator.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/bambam/ReadEndsContainer.hpp>
#include <libmaus2/bambam/SortedFragDecoder.hpp>
#include <libmaus2/bitio/BitVector.hpp>
#include <libmaus2/fastx/FastATwoBitTable.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/lz/BgzfInflateDeflateParallel.hpp>
#include <libmaus2/lz/BgzfParallelRecodeDeflateBase.hpp>
#include <libmaus2/lz/BgzfRecode.hpp>
#include <libmaus2/lz/BgzfRecodeParallel.hpp>
#include <libmaus2/math/iabs.hpp>
#include <libmaus2/math/numbits.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/trie/SimpleTrie.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ContainerGetObject.hpp>
#include <libmaus2/util/MemUsage.hpp>
#include <biobambam2/Licensing.hpp>
#include <libmaus2/lz/GzipOutputStream.hpp>
#include <libmaus2/lz/BufferedGzipStream.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>

// static std::string formatNumber(int64_t const n) { std::ostringstream ostr; ostr << n; return ostr.str(); }
static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1048576;  }
static bool getDefaultRewriteBam() { return 0; }
static int getDefaultRewriteBamLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultColHashBits() { return 20; }
static uint64_t getDefaultColListSize() { return 32*1024*1024; }
static uint64_t getDefaultFragBufSize() { return 48*1024*1024; }
static bool getDefaultRmDup() { return 0; }
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }
static int getDefaultOptMinPixelDif() { return 100; }
static std::string getDefaultODTag() { return "od"; }
static std::string getProgId() { return "bammarkduplicatesopt"; }
static std::string getDefaultInputFormat() { return "bam"; }
static int getDefaultAddMateCigar() { return 0; }

struct OptEntry
{
	typedef OptEntry this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t rank;
	uint64_t refrank;

	OptEntry()
	{
	}

	OptEntry(uint64_t const rrank, uint64_t const rrefrank) : rank(rrank), refrank(rrefrank) {}
	OptEntry(std::istream & in)
	: rank(libmaus2::util::NumberSerialisation::deserialiseNumber(in)), refrank(libmaus2::util::NumberSerialisation::deserialiseNumber(in))
	{

	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,rank);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,refrank);
		return out;
	}

	std::istream & deserialise(std::istream & in)
	{
		*this = OptEntry(in);
		return in;
	}
};

struct OptEntryRefRankCmp
{
	bool operator()(OptEntry const & A, OptEntry const & B) const
	{
		if ( A.refrank != B.refrank )
			return A.refrank < B.refrank;
		else
			return A.rank < B.rank;
	}
};

struct OptMark : public libmaus2::bambam::DupMarkBase::MarkOptical
{
	typedef OptMark this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::string const fn;
	typedef libmaus2::sorting::SerialisingSortingBufferedOutputFile<OptEntry,OptEntryRefRankCmp> sort_file_type;
	sort_file_type::unique_ptr_type Psortfile;
	typedef sort_file_type::merger_ptr_type merger_ptr_type;

	OptMark(std::string const & rfn)
	: fn(rfn), Psortfile(new sort_file_type(fn))
	{

	}

	virtual void operator()(uint64_t const readid, uint64_t const refreadid)
	{
		// std::cerr << "[V] adding opt entry " << readid << "," << refreadid << std::endl;
		Psortfile->put(OptEntry(readid,refreadid));
	}

	merger_ptr_type getMerger()
	{
		merger_ptr_type merger(Psortfile->getMerger(64*1024));
		return UNIQUE_PTR_MOVE(merger);
	}
};

struct NameInput
{
	typedef NameInput this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::string const fn;
	libmaus2::aio::InputStreamInstance ISI;
	libmaus2::lz::BufferedGzipStream BGS;

	NameInput(std::string const & rfn) : fn(rfn), ISI(fn), BGS(ISI)
	{

	}

	bool get(std::string & s)
	{
		if ( BGS.peek() != std::istream::traits_type::eof() )
		{
			s = libmaus2::util::StringSerialisation::deserialiseString(BGS);
			return true;
		}
		else
		{
			return false;
		}
	}
};

struct NameArray
{
	typedef NameArray this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	NameInput NI;
	int64_t rank;
	bool rankvalid;
	std::string name;

	NameArray(std::string const & fn)
	: NI(fn), rank(-1), rankvalid(false)
	{

	}

	std::string operator[](uint64_t const r)
	{
		while ( rank < static_cast<int64_t>(r) )
		{
			bool const ok = NI.get(name);
			assert ( ok );
			rank += 1;
		}
		assert ( rank == static_cast<int64_t>(r) );
		return name;
	}
};

struct BamInputCallback : public libmaus2::bambam::CollatingBamDecoderAlignmentInputCallback
{
	typedef BamInputCallback this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::string const namefn;
	libmaus2::aio::OutputStreamInstance::unique_ptr_type POSI;
	libmaus2::lz::GzipOutputStream::unique_ptr_type PGZ;
	std::ostream & OSI;
	uint64_t rank;

	BamInputCallback(std::string const & rnamefn)
	: namefn(rnamefn), POSI(new libmaus2::aio::OutputStreamInstance(namefn)), PGZ(new libmaus2::lz::GzipOutputStream(*POSI)), OSI(*PGZ), rank(0)
	{

	}

	~BamInputCallback()
	{
		OSI.flush();
		PGZ.reset();
		POSI.reset();
	}

	virtual void operator()(::libmaus2::bambam::BamAlignment const & A)
	{
		// uint64_t const lrank = rank++;
		libmaus2::util::StringSerialisation::serialiseString(OSI,A.getName());
		// std::cerr << lrank << " " << A.getName() << std::endl;
	}
};

static int markDuplicatesOpt(::libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::timing::RealTimeClock globrtc; globrtc.start();

	::libmaus2::util::TempFileRemovalContainer::setup();

	if ( (!(arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != ""))) && isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "refusing to read compressed data from terminal. please use I=<filename> or redirect standard input to a file" << std::endl;
		se.finish();
		throw se;
	}

	if (
		(!(arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != ""))) && isatty(STDOUT_FILENO)
		&&
		(!arginfo.hasArg("outputformat") || arginfo.getUnparsedValue("outputformat",std::string()) != "sam" )
	)
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "refusing to write compressed data to terminal. please use O=<filename> or redirect standard output to a file" << std::endl;
		se.finish();
		throw se;
	}

	// logarithm of collation hash table size
	unsigned int const colhashbits = arginfo.getValue<unsigned int>("colhashbits",getDefaultColHashBits());
	// length of collation output list
	uint64_t const collistsize = arginfo.getValueUnsignedNumeric<uint64_t>("collistsize",getDefaultColListSize());
	// buffer size for fragment and pair data
	uint64_t const fragbufsize = arginfo.getValueUnsignedNumeric<uint64_t>("fragbufsize",getDefaultFragBufSize());
	// print verbosity messages
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	// rewritten file should be in bam format, if input is given via stdin
	unsigned int const rewritebam = arginfo.getValue<unsigned int>("rewritebam",getDefaultRewriteBam());
	int const rewritebamlevel = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("rewritebamlevel",getDefaultRewriteBamLevel()));
	int const optminpixeldif = arginfo.getValue<unsigned int>("optminpixeldif",getDefaultOptMinPixelDif());
	std::string const odtag = arginfo.getUnparsedValue("odtag",getDefaultODTag());
	int const addmatecigar = arginfo.getValue<unsigned int>("addmatecigar",getDefaultAddMateCigar());

	if ( odtag.size() != 2 )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "[E] invalid odtag value (length is not two)" << std::endl;
		se.finish();
		throw se;
	}
	bool const odfirstvalid = (odtag[0] >= 'a' && odtag[0] <= 'z') || (odtag[0] >= 'A' && odtag[0] <= 'Z');
	bool const odsecondvalid = (odtag[1] >= 'a' && odtag[1] <= 'z') || (odtag[1] >= 'A' && odtag[1] <= 'Z') || (odtag[1] >= '0' && odtag[1] <= '9');
	bool const odvalid = odfirstvalid && odsecondvalid;
	if ( ! odvalid )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "[E] invalid odtag value (see SAM spec)" << std::endl;
		se.finish();
		throw se;
	}

	std::cerr << "[V] optminpixeldif=" << optminpixeldif << std::endl;

	enum tag_type_enum
	{
		tag_type_none,
		tag_type_string,
		tag_type_nucleotide
	};

	// tag field
	bool const havetag = arginfo.hasArg("tag");
	std::string const tag = arginfo.getUnparsedValue("tag","no tag");
	libmaus2::trie::SimpleTrie::unique_ptr_type Ptagtrie;

	if ( havetag && (tag.size() != 2 || (!isalpha(tag[0])) || (!isalnum(tag[1])) ) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "tag " << tag << " is invalid" << std::endl;
		se.finish();
		throw se;
	}

	if ( havetag )
	{
		libmaus2::trie::SimpleTrie::unique_ptr_type Ttagtrie(new libmaus2::trie::SimpleTrie);
		Ptagtrie = UNIQUE_PTR_MOVE(Ttagtrie);

		// allocate tag id 0 for empty tag
		uint8_t const * p = 0;
		Ptagtrie->insert(p,p);
	}

	// nucl tag field
	bool const havenucltag = arginfo.hasArg("nucltag");
	std::string const nucltag = arginfo.getUnparsedValue("nucltag","no tag");

	if ( havenucltag && (nucltag.size() != 2 || (!isalpha(nucltag[0])) || (!isalnum(nucltag[1])) ) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "nucltag " << tag << " is invalid" << std::endl;
		se.finish();
		throw se;
	}

	if ( havetag && havenucltag )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "tag and nucltag are mutually exclusive" << std::endl;
		se.finish();
		throw se;
	}

	tag_type_enum tag_type;

	if ( havetag )
		tag_type = tag_type_string;
	else if ( havenucltag )
		tag_type = tag_type_nucleotide;
	else
		tag_type = tag_type_none;

	char const * ctag = havetag ? tag.c_str() : 0;
	char const * cnucltag = havenucltag ? nucltag.c_str() : 0;
	char const * tag1 = 0;
	char const * tag2 = 0;
	libmaus2::autoarray::AutoArray<char> tagbuffer;
	uint64_t taglen = 0;
	uint64_t tagid = 0;
	libmaus2::fastx::FastATwoBitTable const FATBT;

	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilename = tmpfilenamebase + "_bamcollate";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const tmpfilenamereadfrags = tmpfilenamebase + "_readfrags";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilenamereadfrags);
	std::string const tmpfilenamereadpairs = tmpfilenamebase + "_readpairs";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilenamereadpairs);
	std::string const tmpfilesnappyreads = tmpfilenamebase + "_alignments";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfilesnappyreads);
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string const nametmpfile = tmpfilenamebase + "_bamnames";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(nametmpfile);
	std::string const matecigartmpfile = tmpfilenamebase + "_matecigar";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(matecigartmpfile);
	std::string const matecigartmptmpfile = tmpfilenamebase + "_matecigartmp";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(matecigartmptmpfile);

	std::string const opttmpfile = tmpfilenamebase + "_opttmpfile";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(opttmpfile);
	std::string const optnametmpfile = tmpfilenamebase + "_optnametmpfile";
	std::string const optnametmptmpfile = tmpfilenamebase + "_optnametmptmpfile";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(optnametmpfile);
	::libmaus2::util::TempFileRemovalContainer::addTempFile(optnametmptmpfile);

	int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));

	if ( verbose )
		std::cerr << "[V] output compression level " << level << std::endl;

	::libmaus2::timing::RealTimeClock fragrtc; fragrtc.start();

	libmaus2::bambam::BamAlignmentInputCallbackSnappy<libmaus2::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type SRC;
	libmaus2::bambam::BamAlignmentInputCallbackBam<libmaus2::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type BWR;
	::libmaus2::aio::InputStreamInstance::unique_ptr_type CIS;
	libmaus2::aio::OutputStreamInstance::unique_ptr_type copybamstr;

	typedef ::libmaus2::bambam::CircularHashCollatingBamDecoder col_base_type;
	typedef col_base_type::unique_ptr_type col_base_ptr_type;
	col_base_ptr_type CBD;

	libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type inwrap;

	if (
		arginfo.getPairCount("I") >= 1
		||
		(arginfo.getPairCount("I") == 0 && rewritebam < 2)
	)
	{
		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type tinwrap(
			libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo,true /* put rank */,0 /* copystr */,std::cin /* istdin */,false,false /* stream */)
		);
		inwrap = UNIQUE_PTR_MOVE(tinwrap);
	}
	else
	{
		assert ( arginfo.getPairCount("I") == 0 && rewritebam >= 2 );

		libmaus2::aio::OutputStreamInstance::unique_ptr_type tcopybamstr(new libmaus2::aio::OutputStreamInstance(tmpfilesnappyreads));
		copybamstr = UNIQUE_PTR_MOVE(tcopybamstr);

		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type tinwrap(
			libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo,true /* put rank */,copybamstr.get() /* copystr */,std::cin /* istdin */,false,false /* stream */)
		);
		inwrap = UNIQUE_PTR_MOVE(tinwrap);
	}

	uint64_t const colexcludeflags =
		libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY |
		libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSUPPLEMENTARY |
		libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FQCFAIL;

	{
		col_base_ptr_type tCBD(new col_base_type(
			inwrap->getDecoder(),
			tmpfilename,
			colexcludeflags,
			colhashbits,
			collistsize
			)
		)
		;
		CBD = UNIQUE_PTR_MOVE(tCBD);
	}

	if ( arginfo.getPairCount("I") == 0 )
	{
		if ( rewritebam == 1 )
		{
			// rewrite file and mark duplicates
			::libmaus2::bambam::BamAlignmentInputCallbackBam<libmaus2::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type tBWR(new ::libmaus2::bambam::BamAlignmentInputCallbackBam<libmaus2::bambam::BamAlignmentInputPositionCallbackNull>(tmpfilesnappyreads,CBD->getHeader(),rewritebamlevel));
			BWR = UNIQUE_PTR_MOVE(tBWR);
			CBD->setInputCallback(BWR.get());

			if ( verbose )
				std::cerr << "[V] Writing bam compressed alignments to file " << tmpfilesnappyreads << std::endl;
		}
		else if ( rewritebam == 0 )
		{
			libmaus2::bambam::BamAlignmentInputCallbackSnappy<libmaus2::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type tSRC(
				new libmaus2::bambam::BamAlignmentInputCallbackSnappy<libmaus2::bambam::BamAlignmentInputPositionCallbackNull>(tmpfilesnappyreads,CBD->getHeader())
			);
			SRC = UNIQUE_PTR_MOVE(tSRC);
			CBD->setInputCallback(SRC.get());
			if ( verbose )
				std::cerr << "[V] Writing snappy compressed alignments to file " << tmpfilesnappyreads << std::endl;
		}
	}

	BamInputCallback::unique_ptr_type BIC(new BamInputCallback(nametmpfile));
	CBD->addInputCallback(BIC.get());

	OptMark::unique_ptr_type Poptmark(new OptMark(opttmpfile));

	::libmaus2::bambam::BamHeader const bamheader = CBD->getHeader();

	typedef col_base_type::alignment_ptr_type alignment_ptr_type;
	std::pair<alignment_ptr_type,alignment_ptr_type> P;
	uint64_t const mod = arginfo.getValue<unsigned int>("mod",getDefaultMod()); // modulus for verbosity
	uint64_t fragcnt = 0; // mapped fragments
	uint64_t paircnt = 0; // mapped pairs
	uint64_t lastproc = 0; // printed at last fragment count

	// bool const copyAlignments = false;
	#if defined(MARKDUPLICATECOPYALIGNMENTS)
	bool const copyAlignments = true;
	#else
	bool const copyAlignments = false;
	#endif
	::libmaus2::bambam::ReadEndsContainer::unique_ptr_type fragREC(new ::libmaus2::bambam::ReadEndsContainer(fragbufsize,tmpfilenamereadfrags,copyAlignments)); // fragment container
	::libmaus2::bambam::ReadEndsContainer::unique_ptr_type pairREC(new ::libmaus2::bambam::ReadEndsContainer(fragbufsize,tmpfilenamereadpairs,copyAlignments)); // pair container

	int64_t maxrank = -1; // maximal appearing rank
	uint64_t als = 0; // number of processed alignments (= mapped+unmapped fragments)
	std::map<uint64_t,::libmaus2::bambam::DuplicationMetrics> metrics;

	::libmaus2::timing::RealTimeClock rtc; rtc.start(); // clock

	// #define UNPAIREDDEBUG

	#if defined(UNPAIREDDEBUG)
	::libmaus2::bambam::BamFormatAuxiliary bamauxiliary;
	#endif

	libmaus2::timing::RealTimeClock readinrtc; readinrtc.start();

	libmaus2::aio::OutputStreamInstance::unique_ptr_type Pmatecigartmp;
	if ( addmatecigar )
	{
		libmaus2::aio::OutputStreamInstance::unique_ptr_type tptr(new libmaus2::aio::OutputStreamInstance(matecigartmpfile));
		Pmatecigartmp = UNIQUE_PTR_MOVE(tptr);
	}

	while ( CBD->tryPair(P) )
	{
		assert ( P.first || P.second );
		uint64_t const lib =
			P.first
			?
			P.first->getLibraryId(bamheader)
			:
			P.second->getLibraryId(bamheader)
			;
		::libmaus2::bambam::DuplicationMetrics & met = metrics[lib];

		if ( addmatecigar )
		{
			if ( P.first && P.second )
			{
				if ( P.first->isMapped() )
				{
					uint64_t const rank = P.second->getRank();
					std::string const scigar = P.first->getCigarString();
					libmaus2::bambam::OptName(rank,scigar).serialise(*Pmatecigartmp);
				}
				if ( P.second->isMapped() )
				{
					uint64_t const rank = P.first->getRank();
					std::string const scigar = P.second->getCigarString();
					libmaus2::bambam::OptName(rank,scigar).serialise(*Pmatecigartmp);
				}
			}
		}

		if ( P.first )
		{
			maxrank = std::max(maxrank,P.first->getRank());
			als++;

			if ( P.first->isUnmap() )
			{
				++met.unmapped;
			}
			else if ( (!P.first->isPaired()) || P.first->isMateUnmap() )
			{
				#if defined(UNPAIREDDEBUG)
				std::cerr << "[D]\t1\t" << P.first->formatAlignment(bamheader,bamauxiliary) << std::endl;
				#endif
				met.unpaired++;
			}
		}
		if ( P.second )
		{
			maxrank = std::max(maxrank,P.second->getRank());
			als++;

			if ( P.second->isUnmap() )
			{
				++met.unmapped;
			}
			else if ( (!P.second->isPaired()) || P.second->isMateUnmap() )
			{
				#if defined(UNPAIREDDEBUG)
				std::cerr << "[D]\t2\t" << P.first->formatAlignment(bamheader,bamauxiliary) << std::endl;
				#endif
				met.unpaired++;
			}
		}

		switch ( tag_type )
		{
			case tag_type_string:
			{
				// length of tags for read1 and read2
				uint64_t l1 = 0, l2 = 0;

				// aux lookup for read1
				if ( P.first )
				{
					tag1 = P.first->getAuxString(ctag);
					l1 = tag1 ? strlen(tag1) : 0;
				}
				// aux lookup for read2
				if ( P.second )
				{
					tag2 = P.second->getAuxString(ctag);
					l2 = tag2 ? strlen(tag2) : 0;
				}

				// length of concatenated tag
				taglen = l1 + l2 + 2;
				// expand buffer if necessary
				if ( taglen > tagbuffer.size() )
					tagbuffer = libmaus2::autoarray::AutoArray<char>(taglen,false);

				// concatenate tags
				char * outptr = tagbuffer.begin();

				memcpy(outptr,tag1,l1);
				outptr += l1;
				*(outptr++) = 0;

				memcpy(outptr,tag2,l2);
				outptr += l2;
				*(outptr++) = 0;

				assert ( outptr - tagbuffer.begin() == static_cast<ptrdiff_t>(taglen) );

				// look up tag id
				tagid = Ptagtrie->insert(
					tagbuffer.begin(),
					outptr
				);

				break;
			}
			case tag_type_nucleotide:
			{
				// aux lookup for read1
				tag1 = P.first ? P.first->getAuxString(cnucltag) : 0;
				// aux lookup for read2
				tag2 = P.second ? P.second->getAuxString(cnucltag) : 0;

				tagid = (FATBT(tag1) << 32) | FATBT(tag2);

				break;
			}
			default:
			{
				tagid = 0;
				break;
			}
		}


		// we are not interested in unmapped reads below, ignore them
		if ( P.first && P.first->isUnmap() )
		{
			P.first = 0;
		}
		if ( P.second && P.second->isUnmap() )
		{
			P.second = 0;
		}

		if ( P.first && P.second )
		{
			#if defined(DEBUG)
			std::cerr << "[V] Got pair for name " << P.first->getName()
				<< "," << P.first->getCoordinate()
				<< "," << P.first->getFlagsS()
				<< "," << P.second->getCoordinate()
				<< "," << P.second->getFlagsS()
				<< std::endl;
			#endif

			met.readpairsexamined++;

			assert ( ! P.first->isUnmap() );
			assert ( ! P.second->isUnmap() );

			// if first appears after second one, then swap the reads, otherwise leave
			if ( !libmaus2::bambam::ReadEndsBase::orderOK(*(P.first),*(P.second)) )
				std::swap(P.first,P.second);

			pairREC->putPair(*(P.first),*(P.second),bamheader,tagid);
			paircnt++;
		}

		if ( P.first )
		{
			fragREC->putFrag(*(P.first),bamheader,tagid);
			fragcnt++;
		}
		if ( P.second )
		{
			fragREC->putFrag(*(P.second),bamheader,tagid);
			fragcnt++;
		}

		if ( verbose && fragcnt/mod != lastproc/mod )
		{
			std::cerr
				<< "[V] "
				<< als << " als, "
				<< fragcnt << " mapped frags, "
				<< paircnt << " mapped pairs, "
				<< fragcnt/rtc.getElapsedSeconds() << " frags/s "
				<< ::libmaus2::util::MemUsage()
				<< " time "
				<< readinrtc.getElapsedSeconds()
				<< " total "
				<< fragrtc.formatTime(fragrtc.getElapsedSeconds())
				<< std::endl;
			readinrtc.start();
			lastproc = fragcnt;
		}
	}

	if ( copybamstr )
	{
		copybamstr->flush();
		copybamstr.reset();
	}

	if ( addmatecigar )
	{
		Pmatecigartmp->flush();
		Pmatecigartmp.reset();

		libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::bambam::OptName>::reduce(std::vector<std::string>(1,matecigartmpfile),matecigartmptmpfile);
		libmaus2::aio::OutputStreamFactoryContainer::rename(matecigartmptmpfile,matecigartmpfile);
	}

	#if 0
	if ( addmatecigar )
	{
		libmaus2::bambam::OptNameReader ONR(matecigartmpfile);
		libmaus2::bambam::OptName ON;
		while ( ONR.peekNext(ON) )
		{
			std::cerr << ON.rank << " " << ON.refreadname << std::endl;
			ONR.getNext(ON);
		}
	}
	#endif

	// number of lines in input file (due to dropping secondary etc. alignments this can be larger than
	// maxrank+1 and larger than als)
	uint64_t const numranks = CBD->getRank(); // maxrank+1;

	CBD.reset();
	CIS.reset();
	SRC.reset();
	BWR.reset();
	BIC.reset();

	fragREC->flush();
	pairREC->flush();
	fragREC->releaseArray();
	pairREC->releaseArray();

	if ( verbose )
		std::cerr << "[V] fragment and pair data computed in time " << fragrtc.getElapsedSeconds() << " (" << fragrtc.formatTime(fragrtc.getElapsedSeconds()) << ")" << std::endl;

	#if 0
	if ( numranks != als )
		std::cerr << "[D] numranks=" << numranks << " != als=" << als << std::endl;

	assert ( numranks == als );
	#endif

	if ( verbose )
		std::cerr
			<< "[V] "
			<< numranks << " lines, "
			<< als << " als, "
			<< fragcnt << " mapped frags, "
			<< paircnt << " mapped pairs, "
			<< fragcnt/rtc.getElapsedSeconds() << " frags/s "
			<< ::libmaus2::util::MemUsage()
			<< std::endl;

	::libmaus2::bambam::DupSetCallbackVector DSCV(numranks,metrics);

	/*
	 * process fragment and pair data to determine which reads are to be marked as duplicates
	 */
	::libmaus2::bambam::ReadEnds nextfrag;
	std::vector< ::libmaus2::bambam::ReadEnds > lfrags;
	uint64_t dupcnt = 0;

	if ( verbose )
		std::cerr << "[V] Checking pairs...";
	rtc.start();
	::libmaus2::bambam::SortedFragDecoder::unique_ptr_type pairDec(pairREC->getDecoder());
	pairREC.reset();
	pairDec->getNext(lfrags);

	while ( pairDec->getNext(nextfrag) )
	{
		if ( ! libmaus2::bambam::DupMarkBase::isDupPair(nextfrag,lfrags.front()) )
		{
			dupcnt += libmaus2::bambam::DupMarkBase::markDuplicatePairsVector(lfrags,DSCV,Poptmark.get(),optminpixeldif);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}

	dupcnt += libmaus2::bambam::DupMarkBase::markDuplicatePairsVector(lfrags,DSCV,Poptmark.get(),optminpixeldif);
	lfrags.resize(0);
	pairDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << paircnt/rtc.getElapsedSeconds() << std::endl;

	if ( verbose )
		std::cerr << "[V] Checking single fragments...";
	rtc.start();
	::libmaus2::bambam::SortedFragDecoder::unique_ptr_type fragDec(fragREC->getDecoder());
	fragREC.reset();
	fragDec->getNext(lfrags);
	while ( fragDec->getNext(nextfrag) )
	{
		if ( !libmaus2::bambam::DupMarkBase::isDupFrag(nextfrag,lfrags.front()) )
		{
			dupcnt += libmaus2::bambam::DupMarkBase::markDuplicateFrags(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	dupcnt += libmaus2::bambam::DupMarkBase::markDuplicateFrags(lfrags,DSCV);
	lfrags.resize(0);
	fragDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << fragcnt/rtc.getElapsedSeconds() << std::endl;

	if ( verbose )
		std::cerr << "[V] number of alignments marked as duplicates: " << DSCV.getNumDups() << " time " << fragrtc.getElapsedSeconds() << " (" << fragrtc.formatTime(fragrtc.getElapsedSeconds()) << ")" << std::endl;
	/*
	 * end of fragment processing
	 */

	libmaus2::aio::OutputStreamInstance::unique_ptr_type PoptnameOSI(new libmaus2::aio::OutputStreamInstance(optnametmpfile));
	NameArray::unique_ptr_type NA(new NameArray(nametmpfile));
	OptMark::merger_ptr_type optmerger(Poptmark->getMerger());
	OptEntry OE;
	while ( optmerger->getNext(OE) )
		libmaus2::bambam::OptName(OE.rank,(*NA)[OE.refrank]).serialise(*PoptnameOSI);
	PoptnameOSI->flush();
	PoptnameOSI.reset();
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::bambam::OptName>::reduce(std::vector<std::string>(1,optnametmpfile),optnametmptmpfile);
	libmaus2::aio::OutputStreamFactoryContainer::rename(optnametmptmpfile,optnametmpfile);

	#if 0
	{
		libmaus2::bambam::OptNameReader ONR(optnametmpfile);
		libmaus2::bambam::OptName ON;
		while ( ONR.peekNext(ON) )
		{
			std::cerr << ON.rank << " " << ON.refreadname << std::endl;
			ONR.getNext(ON);
		}
	}
	#endif

	/**
	 * write metrics
	 **/
	::libmaus2::aio::OutputStreamInstance::unique_ptr_type pM;
	std::ostream * pmetricstr = 0;

	if ( arginfo.hasArg("M") && (arginfo.getValue<std::string>("M","") != "") )
	{
		::libmaus2::aio::OutputStreamInstance::unique_ptr_type tpM(
                                new ::libmaus2::aio::OutputStreamInstance(arginfo.getValue<std::string>("M",std::string("M")))
                        );
		pM = UNIQUE_PTR_MOVE(tpM);
		pmetricstr = pM.get();
	}
	else
	{
		pmetricstr = & std::cerr;
	}

	std::ostream & metricsstr = *pmetricstr;

	::libmaus2::bambam::DuplicationMetrics::printFormatHeader(arginfo.commandline,metricsstr);
	for ( std::map<uint64_t,::libmaus2::bambam::DuplicationMetrics>::const_iterator ita = metrics.begin(); ita != metrics.end();
		++ita )
		ita->second.format(metricsstr, bamheader.getLibraryName(ita->first));

	if ( metrics.size() == 1 )
	{
		metricsstr << std::endl;
		metricsstr << "## HISTOGRAM\nBIN\tVALUE" << std::endl;
		metrics.begin()->second.printHistogram(metricsstr);
	}

	metricsstr.flush();
	pM.reset();
	/*
	 * end of metrics file writing
	 */


	/*
	 * mark the duplicates
	 */
	libmaus2::bambam::DupMarkBase::markOptDuplicatesInFile(
		arginfo,verbose,bamheader,maxrank,mod,DSCV,tmpfilesnappyreads,rewritebam,tmpfileindex,
		getProgId(),
		std::string(PACKAGE_VERSION),
		getDefaultRmDup(),
		getDefaultMD5(),
		getDefaultIndex(),
		optnametmpfile,
		odtag,
		addmatecigar ? matecigartmpfile : std::string()
	);

	if ( verbose )
		std::cerr << "[V] " << ::libmaus2::util::MemUsage() << " "
			<< globrtc.getElapsedSeconds()
			<< " ("
			<< globrtc.formatTime(globrtc.getElapsedSeconds())
			<< ")"
			<< std::endl;

	return EXIT_SUCCESS;
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
				std::cerr << ::biobambam2::Licensing::license();
				std::cerr << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;

				std::vector< std::pair<std::string,std::string> > V;

				V.push_back ( std::pair<std::string,std::string> ( "I=<filename>", "input file, stdin if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<filename>", "output file, stdout if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "M=<filename>", "metrics file, stderr if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "D=<filename>", "duplicates output file if rmdup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report (default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam2::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for each mod'th record/alignment" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebam=<["+::biobambam2::Licensing::formatNumber(getDefaultRewriteBam())+"]>", "compression of temporary alignment file when input is via stdin (0=snappy,1=gzip/bam,2=copy)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebamlevel=<["+::biobambam2::Licensing::formatNumber(getDefaultRewriteBamLevel())+"]>", std::string("compression setting for rewritten input file if rewritebam=1 (") + libmaus2::bambam::BamBlockWriterBaseFactory::getLevelHelpText() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rmdup=<["+::biobambam2::Licensing::formatNumber(getDefaultRmDup())+"]>", "remove duplicates (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tag=<[a-zA-Z][a-zA-Z0-9]>", "aux field id for tag string extraction" ) );
				V.push_back ( std::pair<std::string,std::string> ( "nucltag=<[a-zA-Z][a-zA-Z0-9]>", "aux field id for nucleotide tag extraction" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhashbits=<["+::biobambam2::Licensing::formatNumber(getDefaultColHashBits())+"]>", "log_2 of size of hash table used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collistsize=<["+::biobambam2::Licensing::formatNumber(getDefaultColListSize())+"]>", "output list size for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "fragbufsize=<["+::biobambam2::Licensing::formatNumber(getDefaultFragBufSize())+"]>", "size of each fragment/pair file buffer in bytes" ) );

				V.push_back ( std::pair<std::string,std::string> ( "D=<filename>", "duplicates output file if rmdup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupmd5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for duplicates output file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupmd5filename=<filename>", "file name for md5 check sum of dup file (default: extend duplicates output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupindex=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index for duplicates file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupindexfilename=<filename>", "file name for BAM index file for duplicates file (default: extend duplicates output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "optminpixeldif=<["+::biobambam2::Licensing::formatNumber(getDefaultOptMinPixelDif())+"]>", "pixel difference threshold for optical duplicates (default: 100)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "odtag=<["+getDefaultODTag()+"]>", "tag added for optical duplicates (default: od)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus2::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus2::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus2::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "addmatecigar=<["+::biobambam2::Licensing::formatNumber(getDefaultAddMateCigar())+"]>", "add mate cigar string field MC (default: 0)" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return markDuplicatesOpt(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
