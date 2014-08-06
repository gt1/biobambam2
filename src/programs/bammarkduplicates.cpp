/**
    bambam
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
#include <config.h>

#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/aio/CheckedInputStream.hpp>
#include <libmaus/aio/CheckedOutputStream.hpp>
#include <libmaus/bambam/BamAlignmentInputCallbackBam.hpp>
#include <libmaus/bambam/BamAlignmentInputCallbackSnappy.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamParallelRewrite.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamAlignmentSnappyInput.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/bambam/CircularHashCollatingBamDecoder.hpp>
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/DuplicationMetrics.hpp>
#include <libmaus/bambam/DupMarkBase.hpp>
#include <libmaus/bambam/DupSetCallbackVector.hpp>
#include <libmaus/bambam/OpticalComparator.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/bambam/ReadEndsContainer.hpp>
#include <libmaus/bambam/SortedFragDecoder.hpp>
#include <libmaus/bitio/BitVector.hpp>
#include <libmaus/fastx/FastATwoBitTable.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/lz/BgzfInflateDeflateParallel.hpp>
#include <libmaus/lz/BgzfParallelRecodeDeflateBase.hpp>
#include <libmaus/lz/BgzfRecode.hpp>
#include <libmaus/lz/BgzfRecodeParallel.hpp>
#include <libmaus/math/iabs.hpp>
#include <libmaus/math/numbits.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/trie/SimpleTrie.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/ContainerGetObject.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <biobambam/Licensing.hpp>

// static std::string formatNumber(int64_t const n) { std::ostringstream ostr; ostr << n; return ostr.str(); }
static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1048576;  }
static bool getDefaultRewriteBam() { return 0; }
static int getDefaultRewriteBamLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultColHashBits() { return 20; }
static uint64_t getDefaultColListSize() { return 32*1024*1024; }
static uint64_t getDefaultFragBufSize() { return 48*1024*1024; }
static uint64_t getDefaultMarkThreads() { return 1; }
static bool getDefaultRmDup() { return 0; }
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }
static std::string getProgId() { return "bammarkduplicates"; }

struct MarkDuplicatesRewriteRequest
{
	libmaus::util::ArgInfo const arginfo;
	bool const verbose;
	libmaus::bambam::BamHeader const bamheader;
	int64_t const maxrank;
	uint64_t const mod;
	int const level;
	::libmaus::bambam::DupSetCallbackVector const DSCV;
	std::string const tmpfilesnappyreads;
	unsigned int const rewritebam;
	std::string const tmpfileindex;
	std::string const progid;
	std::string const packageversion;
	bool const rmdup;
	bool const md5;
	bool const index;
	uint64_t const markthreads;
	
	MarkDuplicatesRewriteRequest(std::istream & in)
	:
		arginfo(in),
		verbose(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		bamheader(libmaus::util::StringSerialisation::deserialiseString(in)),
		maxrank(libmaus::util::NumberSerialisation::deserialiseSignedNumber(in)),
		mod(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		level(libmaus::util::NumberSerialisation::deserialiseSignedNumber(in)),
		DSCV(in),
		tmpfilesnappyreads(libmaus::util::StringSerialisation::deserialiseString(in)),
		rewritebam(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		tmpfileindex(libmaus::util::StringSerialisation::deserialiseString(in)),
		progid(libmaus::util::StringSerialisation::deserialiseString(in)),
		packageversion(libmaus::util::StringSerialisation::deserialiseString(in)),
		rmdup(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		md5(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		index(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		markthreads(libmaus::util::NumberSerialisation::deserialiseNumber(in))
	{	
	}

	static void serialise(
		std::ostream & out,
		libmaus::util::ArgInfo const & arginfo,
		bool const verbose,
		libmaus::bambam::BamHeader const & bamheader,
		int64_t const maxrank,
		uint64_t const mod,
		int const level,
		::libmaus::bambam::DupSetCallbackVector const & DSCV,
		std::string const & tmpfilesnappyreads,
		unsigned int const rewritebam,
		std::string const & tmpfileindex,
		std::string const & progid,
		std::string const & packageversion,
		bool const rmdup,
		bool const md5,
		bool const index,
		uint64_t const markthreads
	)
	{
		arginfo.serialise(out);
		libmaus::util::NumberSerialisation::serialiseNumber(out,verbose);
		libmaus::util::StringSerialisation::serialiseString(out,bamheader.text);
		libmaus::util::NumberSerialisation::serialiseSignedNumber(out,maxrank);
		libmaus::util::NumberSerialisation::serialiseNumber(out,mod);
		libmaus::util::NumberSerialisation::serialiseSignedNumber(out,level);
		DSCV.serialise(out);
		libmaus::util::StringSerialisation::serialiseString(out,tmpfilesnappyreads);
		libmaus::util::NumberSerialisation::serialiseNumber(out,rewritebam);
		libmaus::util::StringSerialisation::serialiseString(out,tmpfileindex);
		libmaus::util::StringSerialisation::serialiseString(out,progid);
		libmaus::util::StringSerialisation::serialiseString(out,packageversion);
		libmaus::util::NumberSerialisation::serialiseNumber(out,rmdup);
		libmaus::util::NumberSerialisation::serialiseNumber(out,md5);
		libmaus::util::NumberSerialisation::serialiseNumber(out,index);
		libmaus::util::NumberSerialisation::serialiseNumber(out,markthreads);
	}
	
	void dispatch() const
	{
		libmaus::bambam::DupMarkBase::markDuplicatesInFile(
			arginfo,verbose,bamheader,maxrank,mod,level,DSCV,tmpfilesnappyreads,rewritebam,tmpfileindex,
			progid,packageversion,rmdup,md5,index,markthreads);
	}
};

static int markDuplicates(::libmaus::util::ArgInfo const & arginfo)
{
	libmaus::timing::RealTimeClock globrtc; globrtc.start();

	::libmaus::util::TempFileRemovalContainer::setup();

	if ( (!(arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != ""))) && isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "refusing to read compressed data from terminal. please use I=<filename> or redirect standard input to a file" << std::endl;
		se.finish();
		throw se;
	}
	
	if ( (!(arginfo.hasArg("O") && (arginfo.getValue<std::string>("O","") != ""))) && isatty(STDOUT_FILENO) )
	{
		::libmaus::exception::LibMausException se;
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
	int const rewritebamlevel = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("rewritebamlevel",getDefaultRewriteBamLevel()));

	enum tag_type_enum
	{
		tag_type_none,
		tag_type_string,
		tag_type_nucleotide
	};

	// tag field
	bool const havetag = arginfo.hasArg("tag");
	std::string const tag = arginfo.getUnparsedValue("tag","no tag");
	libmaus::trie::SimpleTrie::unique_ptr_type Ptagtrie;

	if ( havetag && (tag.size() != 2 || (!isalpha(tag[0])) || (!isalnum(tag[1])) ) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "tag " << tag << " is invalid" << std::endl;
		se.finish();
		throw se;			
	}
	
	if ( havetag )
	{
		libmaus::trie::SimpleTrie::unique_ptr_type Ttagtrie(new libmaus::trie::SimpleTrie);
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
		::libmaus::exception::LibMausException se;
		se.getStream() << "nucltag " << tag << " is invalid" << std::endl;
		se.finish();
		throw se;			
	}
	
	if ( havetag && havenucltag )
	{
		::libmaus::exception::LibMausException se;
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
	libmaus::autoarray::AutoArray<char> tagbuffer;
	uint64_t taglen = 0;
	uint64_t tagid = 0;
	libmaus::fastx::FastATwoBitTable const FATBT;

	// prefix for tmp files
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfilename = tmpfilenamebase + "_bamcollate";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilename);
	std::string const tmpfilenamereadfrags = tmpfilenamebase + "_readfrags";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenamereadfrags);
	std::string const tmpfilenamereadpairs = tmpfilenamebase + "_readpairs";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilenamereadpairs);
	std::string const tmpfilesnappyreads = tmpfilenamebase + "_alignments";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfilesnappyreads);
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	
	if ( verbose )
		std::cerr << "[V] output compression level " << level << std::endl;

	::libmaus::timing::RealTimeClock fragrtc; fragrtc.start();

	libmaus::bambam::BamAlignmentInputCallbackSnappy<libmaus::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type SRC;
	libmaus::bambam::BamAlignmentInputCallbackBam<libmaus::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type BWR;
	::libmaus::aio::CheckedInputStream::unique_ptr_type CIS;
	libmaus::aio::CheckedOutputStream::unique_ptr_type copybamstr;

	typedef ::libmaus::bambam::BamCircularHashCollatingBamDecoder col_type;
	typedef ::libmaus::bambam::BamParallelCircularHashCollatingBamDecoder par_col_type;
	typedef ::libmaus::bambam::CircularHashCollatingBamDecoder col_base_type;
	typedef ::libmaus::bambam::BamMergeCoordinateCircularHashCollatingBamDecoder merge_col_type;
	typedef col_base_type::unique_ptr_type col_base_ptr_type;	
	col_base_ptr_type CBD;

	uint64_t const markthreads = 
		std::max(static_cast<uint64_t>(1),arginfo.getValue<uint64_t>("markthreads",getDefaultMarkThreads()));

	uint64_t const colexcludeflags =
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY |
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSUPPLEMENTARY |
		libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL;

	if ( arginfo.getPairCount("I") > 1 )
	{
		std::vector<std::string> const inputfilenames = arginfo.getPairValues("I");
		col_base_ptr_type tCBD(new merge_col_type(inputfilenames,tmpfilename,colexcludeflags,true,colhashbits,collistsize));
		CBD = UNIQUE_PTR_MOVE(tCBD);
	}
	// if we are reading the input from a file
	else if ( arginfo.hasArg("I") && (arginfo.getValue<std::string>("I","") != "") )
	{
		std::string const inputfilename = arginfo.getValue<std::string>("I","I");
		::libmaus::aio::CheckedInputStream::unique_ptr_type tCIS(new ::libmaus::aio::CheckedInputStream(inputfilename));
		CIS = UNIQUE_PTR_MOVE(tCIS);
		
		if ( markthreads > 1 )
		{
			col_base_ptr_type tCBD(new par_col_type(
                                *CIS,
                                markthreads,
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));
			CBD = UNIQUE_PTR_MOVE(tCBD);
		}
		else
		{
			col_base_ptr_type tCBD(new col_type(
                                *CIS,
                                // numthreads
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));	
			CBD = UNIQUE_PTR_MOVE(tCBD);
		}
	}
	// not a file, we are reading from standard input
	else
	{
		// rewrite to bam
		if ( rewritebam )
		{
			if ( rewritebam > 1 )
			{
				libmaus::aio::CheckedOutputStream::unique_ptr_type tcopybamstr(new libmaus::aio::CheckedOutputStream(tmpfilesnappyreads));
				copybamstr = UNIQUE_PTR_MOVE(tcopybamstr);

				if ( markthreads > 1 )
				{
					col_base_ptr_type tCBD(new par_col_type(std::cin,
                                                *copybamstr,
                                                markthreads,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}
				else
				{
					col_base_ptr_type tCBD(new col_type(std::cin,
                                                *copybamstr,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}

				if ( verbose )
					std::cerr << "[V] Copying bam compressed alignments to file " << tmpfilesnappyreads << std::endl;
			}
			else
			{
				if ( markthreads > 1 )
				{
					col_base_ptr_type tCBD(new par_col_type(std::cin,markthreads,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}
				else
				{
					col_base_ptr_type tCBD(new col_type(std::cin,
                                                tmpfilename,
                                                colexcludeflags,
                                                true /* put rank */,
                                                colhashbits,collistsize));
					CBD = UNIQUE_PTR_MOVE(tCBD);
				}

				// rewrite file and mark duplicates
				::libmaus::bambam::BamAlignmentInputCallbackBam<libmaus::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type tBWR(new ::libmaus::bambam::BamAlignmentInputCallbackBam<libmaus::bambam::BamAlignmentInputPositionCallbackNull>(tmpfilesnappyreads,CBD->getHeader(),rewritebamlevel));
				BWR = UNIQUE_PTR_MOVE(tBWR);
				CBD->setInputCallback(BWR.get());

				if ( verbose )
					std::cerr << "[V] Writing bam compressed alignments to file " << tmpfilesnappyreads << std::endl;
			}
		}
		else
		{
			col_base_ptr_type tCBD(new col_type(std::cin,
                                tmpfilename,
                                colexcludeflags,
                                true /* put rank */,
                                colhashbits,collistsize));
			CBD = UNIQUE_PTR_MOVE(tCBD);

			libmaus::bambam::BamAlignmentInputCallbackSnappy<libmaus::bambam::BamAlignmentInputPositionCallbackNull>::unique_ptr_type tSRC(
				new libmaus::bambam::BamAlignmentInputCallbackSnappy<libmaus::bambam::BamAlignmentInputPositionCallbackNull>(tmpfilesnappyreads,CBD->getHeader())
			);
			SRC = UNIQUE_PTR_MOVE(tSRC);
			CBD->setInputCallback(SRC.get());
			if ( verbose )
				std::cerr << "[V] Writing snappy compressed alignments to file " << tmpfilesnappyreads << std::endl;
		}
	}
	
	::libmaus::bambam::BamHeader const bamheader = CBD->getHeader();

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
	::libmaus::bambam::ReadEndsContainer::unique_ptr_type fragREC(new ::libmaus::bambam::ReadEndsContainer(fragbufsize,tmpfilenamereadfrags,copyAlignments)); // fragment container
	::libmaus::bambam::ReadEndsContainer::unique_ptr_type pairREC(new ::libmaus::bambam::ReadEndsContainer(fragbufsize,tmpfilenamereadpairs,copyAlignments)); // pair container
	
	int64_t maxrank = -1; // maximal appearing rank
	uint64_t als = 0; // number of processed alignments (= mapped+unmapped fragments)
	std::map<uint64_t,::libmaus::bambam::DuplicationMetrics> metrics;

	::libmaus::timing::RealTimeClock rtc; rtc.start(); // clock

	// #define UNPAIREDDEBUG

	#if defined(UNPAIREDDEBUG)
	::libmaus::bambam::BamFormatAuxiliary bamauxiliary;
	#endif
	
	libmaus::timing::RealTimeClock readinrtc; readinrtc.start();
	
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
		::libmaus::bambam::DuplicationMetrics & met = metrics[lib];
		
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
					tagbuffer = libmaus::autoarray::AutoArray<char>(taglen,false);

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
				if ( P.first )
					tag1 = P.first->getAuxString(cnucltag);
				// aux lookup for read2
				if ( P.second )
					tag2 = P.second->getAuxString(cnucltag);

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
			if ( 
				P.second->getRefID() > P.first->getRefID()	
				||
				(
					P.second->getRefID() == P.first->getRefID() 
					&&
					P.second->getCoordinate() >= P.first->getCoordinate()
				)
			)
			{
			
			}
			else
			{
				std::swap(P.first,P.second);
			}
		
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
				<< ::libmaus::util::MemUsage()
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

	// number of lines in input file (due to dropping secondary etc. alignments this can be larger than
	// maxrank+1 and larger than als)
	uint64_t const numranks = CBD->getRank(); // maxrank+1;
	
	CBD.reset();
	CIS.reset();
	SRC.reset();
	BWR.reset();
			
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
			<< ::libmaus::util::MemUsage()
			<< std::endl;

	::libmaus::bambam::DupSetCallbackVector DSCV(numranks,metrics);

	/*
	 * process fragment and pair data to determine which reads are to be marked as duplicates
	 */		
	::libmaus::bambam::ReadEnds nextfrag;
	std::vector< ::libmaus::bambam::ReadEnds > lfrags;
	uint64_t dupcnt = 0;

	if ( verbose )
		std::cerr << "[V] Checking pairs...";
	rtc.start();
	::libmaus::bambam::SortedFragDecoder::unique_ptr_type pairDec(pairREC->getDecoder());
	pairREC.reset();
	pairDec->getNext(lfrags);

	while ( pairDec->getNext(nextfrag) )
	{
		if ( ! libmaus::bambam::DupMarkBase::isDupPair(nextfrag,lfrags.front()) )
		{
			dupcnt += libmaus::bambam::DupMarkBase::markDuplicatePairsVector(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	dupcnt += libmaus::bambam::DupMarkBase::markDuplicatePairsVector(lfrags,DSCV);
	lfrags.resize(0);
	pairDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << paircnt/rtc.getElapsedSeconds() << std::endl;

	if ( verbose )
		std::cerr << "[V] Checking single fragments...";
	rtc.start();
	::libmaus::bambam::SortedFragDecoder::unique_ptr_type fragDec(fragREC->getDecoder());
	fragREC.reset();
	fragDec->getNext(lfrags);
	while ( fragDec->getNext(nextfrag) )
	{
		if ( !libmaus::bambam::DupMarkBase::isDupFrag(nextfrag,lfrags.front()) )
		{
			dupcnt += libmaus::bambam::DupMarkBase::markDuplicateFrags(lfrags,DSCV);
			lfrags.resize(0);
		}

		lfrags.push_back(nextfrag);
	}
	dupcnt += libmaus::bambam::DupMarkBase::markDuplicateFrags(lfrags,DSCV);
	lfrags.resize(0);
	fragDec.reset();
	if ( verbose )
		std::cerr << "done, rate " << fragcnt/rtc.getElapsedSeconds() << std::endl;		

	if ( verbose )
		std::cerr << "[V] number of alignments marked as duplicates: " << DSCV.getNumDups() << " time " << fragrtc.getElapsedSeconds() << " (" << fragrtc.formatTime(fragrtc.getElapsedSeconds()) << ")" << std::endl;
	/*
	 * end of fragment processing
	 */

	/**
	 * write metrics
	 **/
	::libmaus::aio::CheckedOutputStream::unique_ptr_type pM;
	std::ostream * pmetricstr = 0;
	
	if ( arginfo.hasArg("M") && (arginfo.getValue<std::string>("M","") != "") )
	{
		::libmaus::aio::CheckedOutputStream::unique_ptr_type tpM(
                                new ::libmaus::aio::CheckedOutputStream(arginfo.getValue<std::string>("M",std::string("M")))
                        );
		pM = UNIQUE_PTR_MOVE(tpM);
		pmetricstr = pM.get();
	}
	else
	{
		pmetricstr = & std::cerr;
	}

	std::ostream & metricsstr = *pmetricstr;

	::libmaus::bambam::DuplicationMetrics::printFormatHeader(arginfo.commandline,metricsstr);
	for ( std::map<uint64_t,::libmaus::bambam::DuplicationMetrics>::const_iterator ita = metrics.begin(); ita != metrics.end();
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
	libmaus::bambam::DupMarkBase::markDuplicatesInFile(
		arginfo,verbose,bamheader,maxrank,mod,level,DSCV,tmpfilesnappyreads,rewritebam,tmpfileindex,
		getProgId(),
		std::string(PACKAGE_VERSION),
		getDefaultRmDup(),
		getDefaultMD5(),
		getDefaultIndex(),
		getDefaultMarkThreads()
	);
		
	if ( verbose )
		std::cerr << "[V] " << ::libmaus::util::MemUsage() << " " 
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
				
				V.push_back ( std::pair<std::string,std::string> ( "I=<filename>", "input file, stdin if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<filename>", "output file, stdout if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "M=<filename>", "metrics file, stderr if unset" ) );
				V.push_back ( std::pair<std::string,std::string> ( "D=<filename>", "duplicates output file if rmdup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "markthreads=<["+::biobambam::Licensing::formatNumber(getDefaultMarkThreads())+"]>", "number of helper threads" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report (default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for each mod'th record/alignment" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebam=<["+::biobambam::Licensing::formatNumber(getDefaultRewriteBam())+"]>", "compression of temporary alignment file when input is via stdin (0=snappy,1=gzip/bam,2=copy)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rewritebamlevel=<["+::biobambam::Licensing::formatNumber(getDefaultRewriteBamLevel())+"]>", std::string("compression setting for rewritten input file if rewritebam=1 (") + libmaus::bambam::BamBlockWriterBaseFactory::getLevelHelpText() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rmdup=<["+::biobambam::Licensing::formatNumber(getDefaultRmDup())+"]>", "remove duplicates (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tag=<[a-zA-Z][a-zA-Z0-9]>", "aux field id for tag string extraction" ) );
				V.push_back ( std::pair<std::string,std::string> ( "nucltag=<[a-zA-Z][a-zA-Z0-9]>", "aux field id for nucleotide tag extraction" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhashbits=<["+::biobambam::Licensing::formatNumber(getDefaultColHashBits())+"]>", "log_2 of size of hash table used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collistsize=<["+::biobambam::Licensing::formatNumber(getDefaultColListSize())+"]>", "output list size for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "fragbufsize=<["+::biobambam::Licensing::formatNumber(getDefaultFragBufSize())+"]>", "size of each fragment/pair file buffer in bytes" ) );

				V.push_back ( std::pair<std::string,std::string> ( "D=<filename>", "duplicates output file if rmdup=1" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupmd5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum for duplicates output file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupmd5filename=<filename>", "file name for md5 check sum of dup file (default: extend duplicates output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupindex=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index for duplicates file (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "dupindexfilename=<filename>", "file name for BAM index file for duplicates file (default: extend duplicates output file name)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return markDuplicates(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

