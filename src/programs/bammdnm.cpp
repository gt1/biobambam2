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
#include <cstdlib>
#include <iostream>
#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/aio/PosixFdOutputStream.hpp>
#include <libmaus/fastx/FastAReader.hpp>
#include <libmaus/fastx/StreamFastAReader.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/GetFileSize.hpp>
#include <libmaus/util/OutputFileNameTools.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

static int getDefaultIndex() { return 0; }
static int getDefaultMD5() { return 0; }
static int getDefaultDisableValidation() { return 0; }
static int getDefaultLevel() { return -1; }
static int getDefaultVerbose() { return 1; }
static int getDefaultRecompIndetOnly() { return 0; }
static int getDefaultWarnChange() { return 0; }

static uint64_t const ioblocksize = 2*1024*1024;

struct MdNmRecalculation
{
	typedef MdNmRecalculation this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	::libmaus::aio::PosixFdInputStream fain;
	libmaus::fastx::StreamFastAReaderWrapper fareader;
	libmaus::fastx::StreamFastAReaderWrapper::pattern_type currefpat;

	bool validate;
	
	int64_t loadrefid;
	int64_t prevcheckrefid;
	int64_t prevcheckpos;
	
	uint64_t numrecalc;
	uint64_t numkept;

	::libmaus::bambam::MdStringComputationContext context;

	libmaus::autoarray::AutoArray<uint8_t> vmap;
	libmaus::autoarray::AutoArray<uint64_t> nar;
	libmaus::rank::ERank222B::unique_ptr_type Prank;
	
	bool recompindetonly;
	bool warnchange;

	MdNmRecalculation(std::string const & reference, bool const rvalidate, bool const rrecompindetonly, bool const rwarnchange)
	: fain(reference,ioblocksize), 
	  fareader(fain), validate(rvalidate), loadrefid(-1), prevcheckrefid(-1), prevcheckpos(-1),
	  numrecalc(0), numkept(0), vmap(256,false),
	  recompindetonly(rrecompindetonly),
	  warnchange(rwarnchange)
	{
		std::fill(vmap.begin(),vmap.end(),1);
		vmap['A'] = vmap['C'] = vmap['G'] = vmap['T'] =
		vmap['a'] = vmap['c'] = vmap['g'] = vmap['t'] = 0;
	}
	
	void checkValidity(uint8_t const * D, uint64_t const blocksize)
	{
		if ( validate )
		{
			libmaus::bambam::libmaus_bambam_alignment_validity const validity =
				libmaus::bambam::BamAlignmentDecoderBase::valid(D,blocksize);

			if ( validity != libmaus::bambam::libmaus_bambam_alignment_validity_ok )
			{
				libmaus::exception::LibMausException se;
				se.getStream() << "Invalid alignment " << validity << std::endl;
				se.finish();
				throw se;
			}	
		}
	}

	void checkOrder(uint8_t const * D)
	{
		// information for this new alignment
		int64_t const thisrefid = libmaus::bambam::BamAlignmentDecoderBase::getRefID(D);
		int64_t const thispos = libmaus::bambam::BamAlignmentDecoderBase::getPos(D);

		// map negative to maximum positive for checking order
		int64_t const thischeckrefid = (thisrefid >= 0) ? thisrefid : std::numeric_limits<int64_t>::max();
		int64_t const thischeckpos   = (thispos >= 0) ? thispos : std::numeric_limits<int64_t>::max();
		
		// true iff order is ok
		bool const orderok =
			(thischeckrefid > prevcheckrefid)
			||
			(thischeckrefid == prevcheckrefid && thischeckpos >= prevcheckpos);

		// throw exception if alignment stream is not sorted by coordinate
		if ( ! orderok )
		{
			::libmaus::exception::LibMausException se;
			se.getStream() << "File is not sorted by coordinate." << std::endl;
			se.finish();
			throw se;
		}
	
		prevcheckrefid = thischeckrefid;
		prevcheckpos = thischeckpos;
	}
	
	std::string const & loadReference(uint64_t const id)
	{
		while ( loadrefid < static_cast<int64_t>(id) )
		{
			bool const ok = fareader.getNextPatternUnlocked(currefpat);
			
			if ( ! ok )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "[D] Failed to load reference sequence id " << id << std::endl;
				se.finish();
				throw se;		
			}
			
			if ( recompindetonly )
			{
				std::string const & seq = currefpat.spattern;
				uint64_t const seqlen = seq.size();
				
				if ( nar.size() < (seqlen+1+63)/64 )
					nar = libmaus::autoarray::AutoArray<uint64_t>((seqlen+1+63)/64,false);
					
				uint8_t const * p = reinterpret_cast<uint8_t const *>(seq.c_str());
				for ( uint64_t i = 0; i < (seqlen/64); ++i )
				{
					uint64_t v = 0;
					for ( uint64_t j = 0; j < 64; ++j )
					{
						v <<= 1;
						v |= vmap[*(p++)];
					}

					nar[i] = v;
				}
				uint64_t const rest = seqlen-((seqlen/64)*64);
				if ( rest )
				{
					uint64_t v = 0;
					for ( uint64_t j = 0; j < rest; ++j )
					{
						v <<= 1;
						v |= vmap[*(p++)];				
					}
					
					nar[seqlen/64] = v << (63-rest);
				}
				
				#if 0
				p = reinterpret_cast<uint8_t const *>(seq.c_str());
				uint64_t nn = 0;
				for ( uint64_t i = 0; i < seqlen; ++i )
				{
					if ( *p == 'N' )
					{
						assert ( libmaus::bitio::getBit(nar.begin(),i) );
						++nn;
					}
					assert ( libmaus::bitio::getBit(nar.begin(),i) == vmap[*(p++)] );
				}
				#endif
				
				libmaus::rank::ERank222B::unique_ptr_type Trank(new libmaus::rank::ERank222B(nar.begin(),64*((seqlen+1+63)/64)));
				Prank = UNIQUE_PTR_MOVE(Trank);
			}
				
			loadrefid += 1;
		}
		
		return currefpat.spattern;
	}

	bool calmdnm(uint8_t const * D, uint64_t const blocksize)
	{
		checkValidity(D,blocksize);
		checkOrder(D);
		
		bool recalc = false;
		
		if ( ! libmaus::bambam::BamAlignmentDecoderBase::isUnmap(libmaus::bambam::BamAlignmentDecoderBase::getFlags(D)) )
		{
			int64_t const refid = libmaus::bambam::BamAlignmentDecoderBase::getRefID(D);
			int64_t const pos = libmaus::bambam::BamAlignmentDecoderBase::getPos(D);
			int64_t const refend = pos + libmaus::bambam::BamAlignmentDecoderBase::getReferenceLength(D);
			std::string const & ref = loadReference(refid);
			
			if ( 
				(
					recompindetonly 
					&& 
					(
						(Prank->rankm1(refend)-Prank->rankm1(pos)) 
						|| 
						libmaus::bambam::BamAlignmentDecoderBase::hasNonACGT(D)
					)
				)
				||
				(!recompindetonly)
			)
			{
				numrecalc += 1;
				libmaus::bambam::BamAlignmentDecoderBase::calculateMd(D,blocksize,context,ref.begin() + pos,warnchange);
				if ( context.diff )
					recalc = true;
			}
			else
			{
				numkept += 1;
			}
		}

		return recalc;
	}
	
	bool calmdnm(libmaus::bambam::BamAlignment const & algn)
	{
		return calmdnm(algn.D.begin(),algn.blocksize);
	}
};

static int bammdnm(libmaus::util::ArgInfo const & arginfo)
{
	libmaus::timing::RealTimeClock rtc;
	rtc.start();
	
	::libmaus::util::TempFileRemovalContainer::setup();

	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	std::string md5filename;
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
	std::string indexfilename;
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

	if ( !arginfo.hasArg("reference") )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Reference key not set, aborting." << std::endl;
		se.finish();
		throw se;
	}

	std::string const reference = arginfo.getUnparsedValue("reference","");

	::libmaus::aio::PosixFdInputStream PFIS(STDIN_FILENO,ioblocksize);
	::libmaus::lz::BgzfInflate< ::libmaus::aio::PosixFdInputStream > infl(PFIS);
	libmaus::lz::BgzfInflateInfo bgzfinfo;
	libmaus::autoarray::AutoArray<char> B(libmaus::lz::BgzfConstants::getBgzfMaxBlockSize(),false);
	bool haveheader = false;
	::libmaus::bambam::BamHeader header;
	::libmaus::bambam::BamHeader::BamHeaderParserState bamheaderparsestate;

	/* parser state types */
	enum parsestate { state_reading_blocklen,  state_post_skip };
	parsestate state = state_reading_blocklen;
	unsigned int blocklenred = 0;
	uint32_t blocklen = 0;
	uint64_t alcnt = 0;

	::libmaus::bambam::BamAlignment algn;
	uint8_t * copyptr = 0;
	
	bool const validate = !arginfo.getValue<unsigned int>("disablevalidation",getDefaultDisableValidation());
	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue("level",Z_DEFAULT_COMPRESSION));
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	bool const recompindetonly = arginfo.getValue<unsigned int>("recompindetonly",getDefaultRecompIndetOnly());
	bool const warnchange = arginfo.getValue<unsigned int>("warnchange",getDefaultWarnChange());

	libmaus::aio::PosixFdOutputStream::unique_ptr_type Ppfos;
	libmaus::bambam::BamWriter::unique_ptr_type Pout;
	libmaus::bambam::BamWriter::stream_type * bamout = 0;
	
	MdNmRecalculation::unique_ptr_type Precalc;

	while ( !(bgzfinfo = infl.readAndInfo(B.begin(),B.size())).streameof )
	{
		uint8_t const * pa = reinterpret_cast<uint8_t *>(B.begin()); // buffer current pointer
		uint8_t const * pc = pa + bgzfinfo.uncompressed; // buffer end pointer
		
		if ( (! haveheader) && (pa != pc) )
		{			
			::libmaus::util::GetObject<uint8_t const *> G(pa);
			std::pair<bool,uint64_t> const P = ::libmaus::bambam::BamHeader::parseHeader(G,bamheaderparsestate,bgzfinfo.uncompressed);

			// header complete?
			if ( P.first )
			{
				header.init(bamheaderparsestate);
				haveheader = true;
				pa = reinterpret_cast<uint8_t *>(B.begin()) + P.second;
				
				::libmaus::bambam::BamHeader::unique_ptr_type uphead(libmaus::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bamlocalrealign",std::string(PACKAGE_VERSION)));

				libmaus::aio::PosixFdOutputStream::unique_ptr_type Tpfos(new libmaus::aio::PosixFdOutputStream(STDOUT_FILENO,ioblocksize));
				Ppfos = UNIQUE_PTR_MOVE(Tpfos);

				libmaus::bambam::BamWriter::unique_ptr_type Tout(new libmaus::bambam::BamWriter(*Ppfos,*uphead,level,Pcbs));
				Pout = UNIQUE_PTR_MOVE(Tout);
				bamout = &(Pout->getStream());
				
				MdNmRecalculation::unique_ptr_type Trecalc(new MdNmRecalculation(reference,validate,recompindetonly,warnchange));
				Precalc = UNIQUE_PTR_MOVE(Trecalc);
			}
		}
		
		if ( (haveheader) && (pa != pc) )
		{
			while ( pa != pc )
			{
				switch ( state )
				{
					/* read length of next alignment block */
					case state_reading_blocklen:
						/* if this is a little endian machine allowing unaligned access */
						#if defined(LIBMAUS_HAVE_i386)
						if ( 
							(!blocklenred) && 
							((pc-pa) >= static_cast<ptrdiff_t>(sizeof(uint32_t))) 
						)
						{
							blocklen = *(reinterpret_cast<uint32_t const *>(pa));
							blocklenred = sizeof(uint32_t);
							pa += sizeof(uint32_t);
							
							
							if ( pc - pa >= blocklen )
							{
								bool const needupdate = Precalc->calmdnm(pa,blocklen);
								if ( needupdate )
								{
									if ( algn.D.size() < blocklen )
										algn.D = ::libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
									algn.blocksize = blocklen;
									
									std::copy(pa,pa+blocklen,algn.D.begin());
									algn.fillMd(Precalc->context);			
									algn.serialise(*bamout);
								}
								else
								{
									bamout->put((blocklen >> 0) & 0xFF);
									bamout->put((blocklen >> 8) & 0xFF);
									bamout->put((blocklen >> 16) & 0xFF);
									bamout->put((blocklen >> 24) & 0xFF);
									bamout->write(reinterpret_cast<char const *>(pa),blocklen);
								}

								alcnt++;
							
								if ( verbose && (alcnt % (1024*1024) == 0) )
									std::cerr << "[V] " << alcnt/(1024*1024) << " " << (alcnt / rtc.getElapsedSeconds()) << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " recalculated=" << Precalc->numrecalc << std::endl;
							
								pa += blocklen;						
								blocklen = 0;
								blocklenred = 0;
							}
							else
							{
								state = state_post_skip;
								if ( algn.D.size() < blocklen )
									algn.D = ::libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
								algn.blocksize = blocklen;
								copyptr = algn.D.begin();
							}
						}
						else
						#endif
						{
							while ( pa != pc && blocklenred < sizeof(uint32_t) )
								blocklen |= static_cast<uint32_t>(*(pa++)) << ((blocklenred++)*8);

							if ( blocklenred == sizeof(uint32_t) )
							{
								if ( pc - pa >= blocklen )
								{
									bool const needupdate = Precalc->calmdnm(pa,blocklen);
									if ( needupdate )
									{
										if ( algn.D.size() < blocklen )
											algn.D = ::libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
										algn.blocksize = blocklen;
										
										std::copy(pa,pa+blocklen,algn.D.begin());
										algn.fillMd(Precalc->context);			
										algn.serialise(*bamout);
									}
									else
									{
										bamout->put((blocklen >> 0) & 0xFF);
										bamout->put((blocklen >> 8) & 0xFF);
										bamout->put((blocklen >> 16) & 0xFF);
										bamout->put((blocklen >> 24) & 0xFF);
										bamout->write(reinterpret_cast<char const *>(pa),blocklen);
									}

									alcnt++;
								
									if ( verbose && (alcnt % (1024*1024) == 0) )
										std::cerr << "[V] " << alcnt/(1024*1024) << " " << (alcnt / rtc.getElapsedSeconds()) << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " recalculated=" << Precalc->numrecalc << std::endl;

									pa += blocklen;			
									blocklen = 0;
									blocklenred = 0;
								}
								else
								{
									state = state_post_skip;
									if ( algn.D.size() < blocklen )
										algn.D = ::libmaus::bambam::BamAlignment::D_array_type(blocklen,false);
									algn.blocksize = blocklen;
									copyptr = algn.D.begin();
								}
							}
						}
						
						break;
					case state_post_skip:
					{
						uint32_t const skip = std::min(
							pc-pa,static_cast<ptrdiff_t>(blocklen)
						);
						std::copy(pa,pa+skip,copyptr);
						copyptr += skip;
						pa += skip;
						blocklen -= skip;
						
						if ( ! blocklen )
						{
							bool const needupdate = Precalc->calmdnm(algn);
							if ( needupdate )
								algn.fillMd(Precalc->context);			

							algn.serialise(*bamout);
							
							// finished an alignment, set up for next one
							state = state_reading_blocklen;
							
							blocklenred = 0;
							blocklen = 0;

							alcnt++;

							if ( verbose && (alcnt % (1024*1024) == 0) )
								std::cerr << "[V] " << alcnt/(1024*1024) << " " << (alcnt / rtc.getElapsedSeconds()) << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " recalculated=" << Precalc->numrecalc << std::endl;
						}
						break;
					}
				}
			}
		}
	}

	Pout.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(std::string(md5filename));
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
	}

	std::cerr << "[V] " << alcnt/(1024*1024) << " " << (alcnt / rtc.getElapsedSeconds()) << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " recalculated=" << Precalc->numrecalc << std::endl;

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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA" ) );
				V.push_back ( std::pair<std::string,std::string> ( "recompindetonly=<["+::biobambam::Licensing::formatNumber(getDefaultRecompIndetOnly())+"]>", "only compute MD/NM fields in the presence of indeterminate bases (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "warnchange=<["+::biobambam::Licensing::formatNumber(getDefaultWarnChange())+"]>", "print a warning message when MD/NM field is present but different from the recomputed value (default: 0)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}


		return bammdnm(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
