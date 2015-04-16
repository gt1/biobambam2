/*
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
*/
#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/bambam/BamAlignmentDecoder.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamHeaderLowMem.hpp>
#include <libmaus/lz/BgzfDeflate.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/lz/BgzfInflateStream.hpp>

#include <biobambam/Licensing.hpp>

#include <config.h>

int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
int getDefaultVerbose() { return 1; }
int getDefaultMD5() { return 0; }

/*
 * compute bit vector of used sequences
 */
static ::libmaus::bitio::IndexedBitVector::unique_ptr_type getUsedSeqVector(libmaus::util::ArgInfo const & arginfo, std::istream & in)
{
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	libmaus::lz::BgzfInflateStream bgzfin(in);

	libmaus::bambam::BamHeaderLowMem::unique_ptr_type PBHLM ( libmaus::bambam::BamHeaderLowMem::constructFromBAM(bgzfin));

	::libmaus::bambam::BamAlignment algn;
	::libmaus::bitio::IndexedBitVector::unique_ptr_type Psqused(new ::libmaus::bitio::IndexedBitVector(PBHLM->getNumRef()));
	::libmaus::bitio::IndexedBitVector & sqused = *Psqused;
	uint64_t c = 0;
	while ( 
		libmaus::bambam::BamAlignmentDecoder::readAlignmentGz(bgzfin,algn)
	)
	{
		if ( algn.isMapped() )
		{
			int64_t const refid = algn.getRefID();
			assert ( refid >= 0 );
			sqused.set(refid,true);
		}
		if ( algn.isPaired() && algn.isMapped() )
		{
			int64_t const refid = algn.getNextRefID();
			assert ( refid >= 0 );
			sqused.set(refid,true);
		}
		
		if ( verbose && (((++c) & (1024*1024-1)) == 0) )
			std::cerr << "[V] " << c/(1024*1024) << std::endl;
	}
	
	sqused.setupIndex();
	
	return UNIQUE_PTR_MOVE(Psqused);
}

static void filterBamUsedSequences(
	libmaus::util::ArgInfo const & arginfo,
	std::istream & in,
	::libmaus::bitio::IndexedBitVector const & IBV,
	std::ostream & out
)
{
	libmaus::lz::BgzfInflateStream bgzfin(in);
	libmaus::bambam::BamHeaderLowMem::unique_ptr_type PBHLM ( libmaus::bambam::BamHeaderLowMem::constructFromBAM(bgzfin));

	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
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

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	libmaus::lz::BgzfDeflate<std::ostream>::unique_ptr_type Pbgzfout(
		new libmaus::lz::BgzfDeflate<std::ostream>(
			out,level
		)
	);
	libmaus::lz::BgzfDeflate<std::ostream> & bgzfout = *Pbgzfout;
	
	if ( verbose )
		std::cerr << "[V] writing filtered header...";
	PBHLM->serialiseSequenceSubset(bgzfout,IBV,"bamfilterheader2" /* id */,"bamfilterheader2" /* pn */,
		arginfo.commandline /* pgCL */, PACKAGE_VERSION /* pgVN */
	);
	if ( verbose )
		std::cerr << "done." << std::endl;

	::libmaus::bambam::BamAlignment algn;
	uint64_t c = 0;
	while ( libmaus::bambam::BamAlignmentDecoder::readAlignmentGz(bgzfin,algn) )
	{
		if ( algn.isMapped() )
		{
			int64_t const refid = algn.getRefID();
			assert ( refid >= 0 );
			assert ( IBV.get(refid) );
			algn.putRefId(IBV.rank1(refid)-1);
		}
		else
		{
			algn.putRefId(-1);
		}
		
		if ( algn.isPaired() && algn.isMapped() )
		{
			int64_t const refid = algn.getNextRefID();
			assert ( refid >= 0 );
			assert ( IBV.get(refid) );
			algn.putNextRefId(IBV.rank1(refid)-1);
		}
		else
		{
			algn.putNextRefId(-1);
		}
		
		algn.serialise(bgzfout);
		
		if ( verbose && ( ((++c) & (1024*1024-1)) == 0 ) )
			std::cerr << "[V] " << c/(1024*1024) << std::endl;
	}
	
	bgzfout.flush();
	bgzfout.addEOFBlock();	
		
	Pbgzfout.reset();

	if ( Pmd5cb )
		Pmd5cb->saveDigestAsFile(md5filename);
}

int bamfilterheader2(libmaus::util::ArgInfo const & arginfo)
{
	std::string const fn = arginfo.getUnparsedRestArg(0);
	
	::libmaus::bitio::IndexedBitVector::unique_ptr_type PIBV;

	// compute vector of used sequences
	{
		libmaus::aio::PosixFdInputStream in(fn);
		::libmaus::bitio::IndexedBitVector::unique_ptr_type TIBV(getUsedSeqVector(arginfo,in));
		PIBV = UNIQUE_PTR_MOVE(TIBV);
	}
	
	// filter file and remove all unused sequences from header
	{
		libmaus::aio::PosixFdInputStream in(fn);
		filterBamUsedSequences(arginfo,in,*PIBV,std::cout);
	}
	
	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::util::ArgInfo const arginfo(argc,argv);
		
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
				V.push_back ( std::pair<std::string,std::string> ( "I=<[input filename]>", "name of the input file (mandatory)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}

		return bamfilterheader2(arginfo);		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
