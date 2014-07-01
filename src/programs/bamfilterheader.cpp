/*
    libmaus
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
*/

#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bitio/BitVector.hpp>
#include <libmaus/util/ArgInfo.hpp>

#include <biobambam/Licensing.hpp>

#include "config.h"

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; };

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>


static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

void getUsedRefSeqs(
	libmaus::util::ArgInfo const & arginfo,
	libmaus::bitio::IndexedBitVector::unique_ptr_type & usedrefseq,
	libmaus::bitio::IndexedBitVector::unique_ptr_type & usedrg,
	libmaus::bambam::BamHeader::unique_ptr_type & uheader
)
{
	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamHeader const & header = dec.getHeader();
	::libmaus::bambam::BamAlignment const & algn = dec.getAlignment();
	uint64_t const numrefseq = header.getNumRef();	
	libmaus::bitio::IndexedBitVector::unique_ptr_type tusedrefseq(new libmaus::bitio::IndexedBitVector(numrefseq));
	libmaus::bitio::IndexedBitVector::unique_ptr_type tusedrg(new libmaus::bitio::IndexedBitVector(header.getNumReadGroups()));
	
	while ( dec.readAlignment() )
	{
		if ( (!algn.isPaired()) && algn.isMapped() )
		{
			assert ( algn.getRefID() >= 0 );
			assert ( algn.getRefID() < static_cast<int64_t>(tusedrefseq->size()) );
			tusedrefseq->set(algn.getRefID(),true);		
		}
		if ( algn.isPaired() && algn.isMapped() )
		{
			assert ( algn.getRefID() >= 0 );
			assert ( algn.getRefID() < static_cast<int64_t>(tusedrefseq->size()) );
			tusedrefseq->set(algn.getRefID(),true);
		}
		if ( algn.isPaired() && algn.isMateMapped() )
		{
			assert ( algn.getNextRefID() >= 0 );
			assert ( algn.getNextRefID() < static_cast<int64_t>(tusedrefseq->size()) );
			tusedrefseq->set(algn.getNextRefID(),true);
		}
			
		int64_t const rgid = header.getReadGroupId(algn.getReadGroup());
		
		if ( rgid >= 0 )
		{
			assert ( rgid < static_cast<int64_t>(header.getNumReadGroups()) );
			tusedrg->set(rgid,true);
		}
	}
	
	tusedrefseq->setupIndex();
	tusedrg->setupIndex();
	
	usedrefseq = UNIQUE_PTR_MOVE(tusedrefseq);
	usedrg = UNIQUE_PTR_MOVE(tusedrg);
	libmaus::bambam::BamHeader::unique_ptr_type tuheader(header.uclone());
	uheader= UNIQUE_PTR_MOVE(tuheader);
}

uint64_t bamheaderfilter(libmaus::util::ArgInfo const & arginfo)
{
	std::string const inputfilename = arginfo.getUnparsedValue("I","");

	if ( ! inputfilename.size() || inputfilename == "-" )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "No input filename given, please set the I key appropriately." << std::endl;
		se.finish();
		throw se;
	}

	libmaus::bitio::IndexedBitVector::unique_ptr_type usedrefseq;
	libmaus::bitio::IndexedBitVector::unique_ptr_type usedrg;
	libmaus::bambam::BamHeader::unique_ptr_type uheader;

	getUsedRefSeqs(arginfo,usedrefseq,usedrg,uheader);

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
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
	/*
	 * end md5/index callbacks
	 */

	std::string headertext(uheader->text);
	std::vector<libmaus::bambam::HeaderLine> hl = libmaus::bambam::HeaderLine::extractLines(headertext);
	
	std::ostringstream headertextostr;
	uint64_t rscnt = 0;
	uint64_t rgcnt = 0;
	for ( uint64_t i = 0; i < hl.size(); ++i )
	{
		if ( hl[i].type == "SQ" )
		{
			if ( usedrefseq->get(rscnt) )
				headertextostr << hl[i].line << std::endl;

			rscnt += 1;
		}
		else if ( hl[i].type == "RG" )
		{
			if ( usedrg->get(rgcnt) )
				headertextostr << hl[i].line << std::endl;
			
			rgcnt += 1;
		}
		else
		{
			headertextostr << hl[i].line << std::endl;
		}
	}
	headertext = headertextostr.str();

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamheaderfilter", // ID
		"bamheaderfilter", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);
	libmaus::bambam::BamBlockWriterBase::unique_ptr_type Pout ( libmaus::bambam::BamBlockWriterBaseFactory::construct(uphead, arginfo, &cbs) );

	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamAlignment & algn = dec.getAlignment();
	
	while ( dec.readAlignment() )
	{
		if ( (!algn.isPaired()) && algn.isMapped() )
		{
			assert ( algn.getRefID() >= 0 );
			assert ( algn.getRefID() < static_cast<int64_t>(usedrefseq->size()) );
			assert ( usedrefseq->get(algn.getRefID()) );
			assert ( usedrefseq->rank1(algn.getRefID())-1 < uphead.getNumRef() );
			algn.putRefId(usedrefseq->rank1(algn.getRefID())-1);
		}
		if ( algn.isPaired() && algn.isMapped() )
		{
			assert ( algn.getRefID() >= 0 );
			assert ( algn.getRefID() < static_cast<int64_t>(usedrefseq->size()) );
			assert ( usedrefseq->get(algn.getRefID()) );
			assert ( usedrefseq->rank1(algn.getRefID())-1 < uphead.getNumRef() );
			algn.putRefId(usedrefseq->rank1(algn.getRefID())-1);
		}
		if ( algn.isPaired() && algn.isMateMapped() )
		{
			assert ( algn.getNextRefID() >= 0 );
			assert ( algn.getNextRefID() < static_cast<int64_t>(usedrefseq->size()) );
			assert ( usedrefseq->get(algn.getNextRefID()) );
			assert ( usedrefseq->rank1(algn.getNextRefID())-1 < uphead.getNumRef() );
			algn.putNextRefId(usedrefseq->rank1(algn.getNextRefID())-1);
		}
		
		// erase unmapped refid and pos
		if ( algn.isUnmap() )
		{
			algn.putRefId(-1);
			algn.putPos(-1);
		}
		if ( algn.isMateUnmap() )
		{
			algn.putNextRefId(-1);
			algn.putNextPos(-1);
		}

		Pout->writeAlignment(algn);
	}


	Pout.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
	}
	
	return 0;
}

int main(int argc, char *argv[])
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
				V.push_back ( std::pair<std::string,std::string> ( "I=<[input filename]>", "name of the input file" ) );
				// V.push_back ( std::pair<std::string,std::string> ( "numthreads=<["+::biobambam::Licensing::formatNumber(getDefaultNumThreads())+"]>", "number of recoding threads" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamheaderfilter(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
