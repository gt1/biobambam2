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

#include <libmaus/bambam/AdapterFilter.hpp>
#include <libmaus/bambam/BamAlignment.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/fastx/acgtnMap.hpp>
#include <libmaus/rank/popcnt.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/Histogram.hpp>


#include <biobambam/Licensing.hpp>
#include <biobambam/Split12.hpp>
#include <biobambam/Strip12.hpp>
#include <biobambam/ClipReinsert.hpp>
#include <biobambam/zzToName.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static uint64_t getDefaultMod() { return 1024*1024; }
static uint64_t getDefaultRankSplit() { return 1; }
static uint64_t getDefaultRankStrip() { return 1; }
static uint64_t getDefaultClipReinsert() { return 1; }
static uint64_t getDefaultZZToName() { return 1; }

#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }

int bam12auxmerge(::libmaus::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	std::string const prefilename = arginfo.getRestArg<std::string>(0);
	libmaus::bambam::BamDecoder bampredec(prefilename);

	int const level = libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const ranksplit = arginfo.getValue<int>("ranksplit",getDefaultRankSplit());
	int const rankstrip = arginfo.getValue<int>("rankstrip",getDefaultRankSplit());
	int const clipreinsert = arginfo.getValue<int>("clipreinsert",getDefaultClipReinsert());
	int const zztoname = arginfo.getValue<int>("zztoname",getDefaultZZToName());
	uint64_t const mod = arginfo.getValue<int>("mod",getDefaultMod());
	uint64_t const bmod = libmaus::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;
	
	libmaus::autoarray::AutoArray<char> Aread;

	::libmaus::bambam::BamDecoder bamdec(std::cin,false);
	::libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	::libmaus::bambam::BamHeader const & preheader = bampredec.getHeader();

	std::string const headertext(header.text);
	std::string const preheadertext(libmaus::bambam::HeaderLine::removeSequenceLines(preheader.text));
	
	libmaus::bambam::ProgramHeaderLineSet headerlines(headertext);
	libmaus::bambam::ProgramHeaderLineSet preheaderlines(preheadertext);
	
	std::vector<libmaus::bambam::HeaderLine> allheaderlines = libmaus::bambam::HeaderLine::extractLines(headertext);
	
	std::string const lastid = preheaderlines.getLastIdInChain();
	
	std::stack < std::pair<uint64_t,std::string> > pgtodo;
	for ( uint64_t i = 0; i < headerlines.roots.size(); ++i )
		pgtodo.push(std::pair<uint64_t,std::string>(headerlines.roots[i],lastid));
	
	std::string upheadtext = preheadertext;
	while ( pgtodo.size() )
	{
		uint64_t const hid = pgtodo.top().first;
		std::string const PP = pgtodo.top().second;
		pgtodo.pop();
		libmaus::bambam::HeaderLine const & line = headerlines.lines[hid];
		
		// ID, PP, PN, CL, VN
		std::string       ID = (line.M.find("ID") != line.M.end()) ? line.M.find("ID")->second : "";
		std::string const PN = (line.M.find("PN") != line.M.end()) ? line.M.find("PN")->second : "";
		std::string const CL = (line.M.find("CL") != line.M.end()) ? line.M.find("CL")->second : "";
		std::string const VN = (line.M.find("VN") != line.M.end()) ? line.M.find("VN")->second : "";
		
		upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLineRef(
			upheadtext,
			ID,
			PN,
			CL,
			PP,
			VN
		);
	
		if ( headerlines.edges.find(hid) != headerlines.edges.end() )
		{
			std::vector<uint64_t> const & children = headerlines.edges.find(hid)->second;
			
			for ( uint64_t j = 0; j < children.size(); ++j )
				pgtodo.push(std::pair<uint64_t,std::string>(children[j],ID));
		}
	}
	
	/* copy SQ lines */
	std::ostringstream sqconcstr;
	sqconcstr << upheadtext;
	for ( uint64_t i = 0; i < allheaderlines.size(); ++i )
		if ( allheaderlines[i].type == "SQ" )
			sqconcstr << allheaderlines[i].line << "\n";
	upheadtext = sqconcstr.str();

	::libmaus::bambam::BamHeader uphead(upheadtext);
	uphead.changeSortOrder("unknown");

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
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5/index callbacks
	 */

	::libmaus::bambam::BamWriter::unique_ptr_type writer(new ::libmaus::bambam::BamWriter(std::cout,uphead,level,Pcbs));
	
	::libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	::libmaus::bambam::BamAlignment & prealgn = bampredec.getAlignment();
	int64_t curid = -1;
	
	libmaus::autoarray::AutoArray< std::pair<uint8_t,uint8_t> > auxpre;
	libmaus::autoarray::AutoArray< std::pair<uint8_t,uint8_t> > auxnew;
	
	libmaus::bambam::BamAuxFilterVector auxfilter;

	// helpers for clipReinsert
	libmaus::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > auxtags;
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	std::stack < libmaus::bambam::cigar_operation > hardstack;
	libmaus::bambam::BamAlignment::D_array_type Tcigar;
	libmaus::bambam::BamAuxFilterVector bafv;
	libmaus::bambam::BamAuxFilterVector auxfilterout;
	auxfilterout.set('q','s');
	auxfilterout.set('q','q');

	// helpers for zztoname	
 	libmaus::bambam::BamAuxFilterVector zzbafv;
 	zzbafv.set('z','z');

	// loop over aligned BAM file
	while ( bamdec.readAlignment() )
	{
		if ( ranksplit )
			split12(algn);
	
		// extract rank
		char const * name = algn.getName();
		char const * u1 = name;
		bool ok = true;
		uint64_t rank = 0;
		while ( *u1 && *u1 != '_' )
		{
			rank *= 10;
			rank += (*u1-'0');
			ok = ok && isdigit(*u1);
			++u1;
		}
		
		// unable to find rank?	write out as is and continue
		if ( ! ok )
		{
			algn.serialise(writer->getStream());
			continue;
		}
		
		// loop over unaligned BAM file
		while ( curid != static_cast<int64_t>(rank) )
		{
			bool const a_ok = bampredec.readAlignment();
			
			if ( ! a_ok )
			{
				libmaus::exception::LibMausException se;
				se.getStream() << "Found unexpected EOF on file " << prefilename << std::endl;
				se.finish();
				throw se;
			}
			assert ( a_ok );
			++curid;
			
			if ( verbose && (! (curid & bmask)) )
				std::cerr << "[V] " << (curid / bmod) << std::endl;
		}
			
		if ( verbose > 1 )
			std::cerr << "Merging:\n" << algn.formatAlignment(header) << "\n" << prealgn.formatAlignment(preheader) << std::endl;
		
		uint64_t pretagnum = prealgn.enumerateAuxTags(auxpre);
		uint64_t newtagnum = algn.enumerateAuxTags(auxnew);
		
		std::sort(auxpre.begin(),auxpre.begin()+pretagnum);
		std::sort(auxnew.begin(),auxnew.begin()+newtagnum);
		
		if ( verbose > 1 )
			std::cerr << "pretagnum=" << pretagnum << " newtagnum=" << newtagnum << std::endl;
		
		std::pair<uint8_t,uint8_t> * prec = auxpre.begin();
		std::pair<uint8_t,uint8_t> * pree = prec + pretagnum;
		std::pair<uint8_t,uint8_t> * preo = prec;

		std::pair<uint8_t,uint8_t> * newc = auxnew.begin();
		std::pair<uint8_t,uint8_t> * newe = newc + newtagnum;
		std::pair<uint8_t,uint8_t> * newo = newc;
		
		while ( prec != pree && newc != newe )
		{
			// pre which is not in new
			if ( *prec < *newc )
			{
				*(preo++) = *(prec++);
			}
			// tag in both, drop pre
			else if ( *prec == *newc )
			{
				*(newo++) = *(newc++);
				prec++;
			}
			// new not in pre
			else
			{
				*(newo++) = *(newc++);					
			}
		}
		
		while ( prec != pree )
			*(preo++) = *(prec++);
		while ( newc != newe )
			*(newo++) = *(newc++);
			
		pretagnum = preo-auxpre.begin();
		newtagnum = newo-auxnew.begin();
		
		for ( uint64_t i = 0; i < pretagnum; ++i )
			auxfilter.set(auxpre[i].first,auxpre[i].second);
		
		algn.copyAuxTags(prealgn, auxfilter);

		for ( uint64_t i = 0; i < pretagnum; ++i )
			auxfilter.clear(auxpre[i].first,auxpre[i].second);

		if ( verbose > 1 )
		{
			std::cerr << "pretagnum=" << pretagnum << " newtagnum=" << newtagnum << std::endl;	
			std::cerr << "result: " << algn.formatAlignment(header) << std::endl;
		}

		// copy QC fail flag from original file to aligner output		
		if ( prealgn.isQCFail() )
			algn.putFlags( algn.getFlags() | libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FQCFAIL );
		
		if ( rankstrip )
			strip12(algn);

		if ( clipreinsert )
			clipReinsert(algn,auxtags,bafv,cigop,Tcigar,hardstack,auxfilterout);
			
		if ( zztoname )
			zzToRank(algn,zzbafv);	
		
		algn.serialise(writer->getStream());
	}

	writer.reset();

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	if ( Pindex )
	{
		Pindex->flush(std::string(indexfilename));
	}
	
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
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for every mod'th alignment if verbose" ) );
				V.push_back ( std::pair<std::string,std::string> ( "ranksplit=<["+::biobambam::Licensing::formatNumber(getDefaultRankSplit())+"]>", "split rank pairs in names (see bam12split command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rankstrip=<["+::biobambam::Licensing::formatNumber(getDefaultRankStrip())+"]>", "strip ranks of names (see bam12strip command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "clipreinsert=<["+::biobambam::Licensing::formatNumber(getDefaultClipReinsert())+"]>", "reinsert clipped sequence fragments (see bamclipreinsert command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "zztoname=<["+::biobambam::Licensing::formatNumber(getDefaultZZToName())+"]>", "move rank from zz to name aux field to name (see bamzztoname command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bam12auxmerge(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

