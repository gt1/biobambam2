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
#include <libmaus2/aio/PosixFdInputStream.hpp>
#include <libmaus2/aio/PosixFdOutputStream.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/fastx/FastABgzfIndex.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/BamHeaderUpdate.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
#include <libmaus2/bambam/MdNmRecalculation.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/lz/RAZFIndex.hpp>
#include <libmaus2/lz/RAZFDecoder.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

static int getDefaultIndex() { return 0; }
static int getDefaultMD5() { return 0; }
static int getDefaultDisableValidation() { return 0; }
static int getDefaultLevel() { return -1; }
static int getDefaultVerbose() { return 1; }
static int getDefaultRecompIndetOnly() { return 0; }
static int getDefaultWarnChange() { return 0; }
static uint64_t getDefaultIOBlockSize() { return 128*1024; }


static int bammdnm(libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::timing::RealTimeClock rtc;
	rtc.start();
	
	::libmaus2::util::TempFileRemovalContainer::setup();

	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	std::string md5filename;
	if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
			md5filename = arginfo.getUnparsedValue("md5filename","");
		else
			std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

		if ( md5filename.size() )
		{
			::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus2::lz::BgzfDeflateOutputCallbackMD5);
			Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
			cbs.push_back(Pmd5cb.get());
		}
	}
	libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Pindex;
	std::string indexfilename;
	if ( arginfo.getValue<unsigned int>("index",getDefaultIndex()) )
	{
		if ( arginfo.hasArg("indexfilename") &&  arginfo.getUnparsedValue("indexfilename","") != "" )
			indexfilename = arginfo.getUnparsedValue("indexfilename","");
		else
			std::cerr << "[V] no filename for index given, not creating index" << std::endl;

		if ( indexfilename.size() )
		{
			libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex::unique_ptr_type Tindex(new libmaus2::bambam::BgzfDeflateOutputCallbackBamIndex(tmpfileindex));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
			cbs.push_back(Pindex.get());
		}
	}
	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;

	if ( !arginfo.hasArg("reference") )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Reference key not set, aborting." << std::endl;
		se.finish();
		throw se;
	}

	std::string const reference = arginfo.getUnparsedValue("reference","");
	uint64_t const ioblocksize = arginfo.getValueUnsignedNumeric<uint64_t>("ioblocksize",getDefaultIOBlockSize());

	::libmaus2::aio::PosixFdInputStream PFIS(STDIN_FILENO,ioblocksize);
	::libmaus2::lz::BgzfInflate< ::libmaus2::aio::PosixFdInputStream > infl(PFIS);
	libmaus2::lz::BgzfInflateInfo bgzfinfo;
	libmaus2::autoarray::AutoArray<char> B(libmaus2::lz::BgzfConstants::getBgzfMaxBlockSize(),false);
	bool haveheader = false;
	::libmaus2::bambam::BamHeader header;
	::libmaus2::bambam::BamHeaderParserState bamheaderparsestate;

	/* parser state types */
	enum parsestate { state_reading_blocklen,  state_post_skip };
	parsestate state = state_reading_blocklen;
	unsigned int blocklenred = 0;
	uint32_t blocklen = 0;
	uint64_t alcnt = 0;

	::libmaus2::bambam::BamAlignment algn;
	uint8_t * copyptr = 0;
	
	bool const validate = !arginfo.getValue<unsigned int>("disablevalidation",getDefaultDisableValidation());
	int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue("level",Z_DEFAULT_COMPRESSION));
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	bool const recompindetonly = arginfo.getValue<unsigned int>("recompindetonly",getDefaultRecompIndetOnly());
	bool const warnchange = arginfo.getValue<unsigned int>("warnchange",getDefaultWarnChange());

	libmaus2::aio::PosixFdOutputStream::unique_ptr_type Ppfos;
	libmaus2::bambam::BamWriter::unique_ptr_type Pout;
	libmaus2::bambam::BamWriter::stream_type * bamout = 0;
	
	libmaus2::bambam::MdNmRecalculation::unique_ptr_type Precalc;

	while ( !(bgzfinfo = infl.readAndInfo(B.begin(),B.size())).streameof )
	{
		uint8_t const * pa = reinterpret_cast<uint8_t *>(B.begin()); // buffer current pointer
		uint8_t const * pc = pa + bgzfinfo.uncompressed; // buffer end pointer
		
		if ( (! haveheader) && (pa != pc) )
		{			
			::libmaus2::util::GetObject<uint8_t const *> G(pa);
			std::pair<bool,uint64_t> const P = bamheaderparsestate.parseHeader(G,bgzfinfo.uncompressed);

			// header complete?
			if ( P.first )
			{
				header.init(bamheaderparsestate);
				haveheader = true;
				pa = reinterpret_cast<uint8_t *>(B.begin()) + P.second;
				
				::libmaus2::bambam::BamHeader::unique_ptr_type uphead(libmaus2::bambam::BamHeaderUpdate::updateHeader(arginfo,header,"bammdnm",std::string(PACKAGE_VERSION)));

				libmaus2::aio::PosixFdOutputStream::unique_ptr_type Tpfos(new libmaus2::aio::PosixFdOutputStream(STDOUT_FILENO,ioblocksize));
				Ppfos = UNIQUE_PTR_MOVE(Tpfos);

				libmaus2::bambam::BamWriter::unique_ptr_type Tout(new libmaus2::bambam::BamWriter(*Ppfos,*uphead,level,Pcbs));
				Pout = UNIQUE_PTR_MOVE(Tout);
				bamout = &(Pout->getStream());
				
				libmaus2::bambam::MdNmRecalculation::unique_ptr_type Trecalc(new libmaus2::bambam::MdNmRecalculation(reference,validate,recompindetonly,warnchange,ioblocksize));
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
						#if defined(LIBMAUS2_HAVE_i386)
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
										algn.D = ::libmaus2::bambam::BamAlignment::D_array_type(blocklen,false);
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
									algn.D = ::libmaus2::bambam::BamAlignment::D_array_type(blocklen,false);
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
											algn.D = ::libmaus2::bambam::BamAlignment::D_array_type(blocklen,false);
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
										algn.D = ::libmaus2::bambam::BamAlignment::D_array_type(blocklen,false);
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
			
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam2::Licensing::formatNumber(getDefaultLevel())+"]>", libmaus2::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam2::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( "disablevalidation=<["+::biobambam2::Licensing::formatNumber(getDefaultDisableValidation())+"]>", "disable input validation (default is 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA" ) );
				V.push_back ( std::pair<std::string,std::string> ( "recompindetonly=<["+::biobambam2::Licensing::formatNumber(getDefaultRecompIndetOnly())+"]>", "only compute MD/NM fields in the presence of indeterminate bases (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "warnchange=<["+::biobambam2::Licensing::formatNumber(getDefaultWarnChange())+"]>", "print a warning message when MD/NM field is present but different from the recomputed value (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "ioblocksize=<["+::biobambam2::Licensing::formatNumber(getDefaultIOBlockSize())+"]>", "block size for I/O operations" ) );

				::biobambam2::Licensing::printMap(std::cerr,V);

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
