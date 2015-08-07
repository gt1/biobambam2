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

#include <libmaus2/aio/OutputStreamInstance.hpp>

#include <libmaus2/bambam/AdapterFilter.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/rank/popcnt.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/Histogram.hpp>


#include <biobambam2/Licensing.hpp>
#include <biobambam2/Split12.hpp>
#include <biobambam2/Strip12.hpp>
#include <biobambam2/ClipReinsert.hpp>
#include <biobambam2/zzToName.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static int getDefaultVerbose() { return 1; }
static int getDefaultSanity() { return 0; }
static uint64_t getDefaultMod() { return 1024*1024; }
static uint64_t getDefaultRankSplit() { return 1; }
static uint64_t getDefaultRankStrip() { return 1; }
static uint64_t getDefaultClipReinsert() { return 1; }
static uint64_t getDefaultZZToName() { return 1; }

#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/bambam/BgzfDeflateOutputCallbackBamIndex.hpp>
static int getDefaultMD5() { return 0; }
static int getDefaultIndex() { return 0; }


bool is_suffix(const char *str, const char *suf)
{
	if ( !str || !suf )
    	    	return false;

	size_t len_str = strlen(str);
	size_t len_suf = strlen(suf);

	if ( len_str < len_suf )
    	    	return false;

	if ( strcmp(str + len_str - len_suf, suf) == 0 )
    	    	return true;
	else
    	    	return false;
} 

int bam12auxmerge(::libmaus2::util::ArgInfo const & arginfo)
{
	if ( isatty(STDIN_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing to read binary data from terminal, please redirect standard input to pipe or file." << std::endl;
		se.finish();
		throw se;
	}

	if ( isatty(STDOUT_FILENO) )
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "Refusing write binary data to terminal, please redirect standard output to pipe or file." << std::endl;
		se.finish();
		throw se;
	}
	
	std::string const prefilename = arginfo.getRestArg<std::string>(0);
	libmaus2::bambam::BamDecoder bampredec(prefilename);

	int const level = libmaus2::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	int const ranksplit = arginfo.getValue<int>("ranksplit",getDefaultRankSplit());
	int const rankstrip = arginfo.getValue<int>("rankstrip",getDefaultRankSplit());
	int const clipreinsert = arginfo.getValue<int>("clipreinsert",getDefaultClipReinsert());
	int const zztoname = arginfo.getValue<int>("zztoname",getDefaultZZToName());
	int const sanity = arginfo.getValue<int>("sanity",getDefaultSanity());
	uint64_t const mod = arginfo.getValue<int>("mod",getDefaultMod());
	uint64_t const bmod = libmaus2::math::nextTwoPow(mod);
	uint64_t const bmask = bmod-1;
	
	libmaus2::autoarray::AutoArray<char> Aread;

	::libmaus2::bambam::BamDecoder bamdec(std::cin,false);
	::libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
	::libmaus2::bambam::BamHeader const & preheader = bampredec.getHeader();

	std::string const headertext(header.text);
	std::string const preheadertext(libmaus2::bambam::HeaderLine::removeSequenceLines(preheader.text));
	
	libmaus2::bambam::ProgramHeaderLineSet headerlines(headertext);
	libmaus2::bambam::ProgramHeaderLineSet preheaderlines(preheadertext);
	
	std::vector<libmaus2::bambam::HeaderLine> allheaderlines = libmaus2::bambam::HeaderLine::extractLines(headertext);
	
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
		libmaus2::bambam::HeaderLine const & line = headerlines.lines[hid];
		
		// ID, PP, PN, CL, VN
		std::string       ID = (line.M.find("ID") != line.M.end()) ? line.M.find("ID")->second : "";
		std::string const PN = (line.M.find("PN") != line.M.end()) ? line.M.find("PN")->second : "";
		std::string const CL = (line.M.find("CL") != line.M.end()) ? line.M.find("CL")->second : "";
		std::string const VN = (line.M.find("VN") != line.M.end()) ? line.M.find("VN")->second : "";
		
		upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLineRef(
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

	::libmaus2::bambam::BamHeader uphead(upheadtext);
	uphead.changeSortOrder("unknown");

	/*
	 * start index/md5 callbacks
	 */
	std::string const tmpfilenamebase = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const tmpfileindex = tmpfilenamebase + "_index";
	::libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfileindex);

	std::string md5filename;
	std::string indexfilename;

	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus2::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
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
	/*
	 * end md5/index callbacks
	 */

	::libmaus2::bambam::BamWriter::unique_ptr_type writer(new ::libmaus2::bambam::BamWriter(std::cout,uphead,level,Pcbs));
	
	::libmaus2::bambam::BamAlignment & algn = bamdec.getAlignment();
	::libmaus2::bambam::BamAlignment & prealgn = bampredec.getAlignment();
	int64_t curid = -1;
	
	libmaus2::autoarray::AutoArray< std::pair<uint8_t,uint8_t> > auxpre;
	libmaus2::autoarray::AutoArray< std::pair<uint8_t,uint8_t> > auxnew;
	
	libmaus2::bambam::BamAuxFilterVector auxfilter;

	// helpers for clipReinsert
	libmaus2::autoarray::AutoArray < std::pair<uint8_t,uint8_t> > auxtags;
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;
	std::stack < libmaus2::bambam::cigar_operation > hardstack;
	libmaus2::bambam::BamAlignment::D_array_type Tcigar;
	libmaus2::bambam::BamAuxFilterVector bafv;
	libmaus2::bambam::BamAuxFilterVector auxfilterout;
	auxfilterout.set('q','s');
	auxfilterout.set('q','q');

	// helpers for zztoname	
 	libmaus2::bambam::BamAuxFilterVector zzbafv;
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
				libmaus2::exception::LibMausException se;
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
		
		// do some sanity checking
		if ( sanity )
		{
		    	// first do a name check
			char const * prename = prealgn.getName();
			u1++; // put on the first letter of readname
			
			if ( verbose > 1 )
    	    	    	    	std::cerr << "Sanity: comparing " << name << " and " << prename << std::endl;			    
			
			if ( !is_suffix(prename, u1) ) // names do not match
			{
			    	libmaus2::exception::LibMausException se;
			    	se.getStream() << "Sanity check failed on read names, found " << name << " and " << prename << std::endl;
				se.finish();
				throw se;
			}
			
			// now the names match so try the flags
			
			if ( !(algn.isPaired() == prealgn.isPaired() &&
			     algn.isRead1() == prealgn.isRead1() &&
			     algn.isRead2() == prealgn.isRead2()) )
			{
			    	libmaus2::exception::LibMausException se;
				se.getStream() << "Sanity check failed on flags, " << std::endl
				    	       << "Aligned " << name << " paired " << algn.isPaired() << " first " << algn.isRead1() << " last " << algn.isRead2() << std::endl
			    	    	       << "Unaligned " << prename << " paired " << prealgn.isPaired() << " first " << prealgn.isRead1() << " last " << prealgn.isRead2() << std::endl;
				se.finish();
				throw se;
			}
			
			if ( verbose > 1 )
			    std::cerr << "Sanity check on flags: " << std::endl
				    	       << "Aligned " << name << " paired " << algn.isPaired() << " first " << algn.isRead1() << " last " << algn.isRead2() << std::endl
			    	    	       << "Unaligned " << prename << " paired " << prealgn.isPaired() << " first " << prealgn.isRead1() << " last " << prealgn.isRead2() << std::endl;	
			
			
		}
		
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
			algn.putFlags( algn.getFlags() | libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FQCFAIL );
		
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
				V.push_back ( std::pair<std::string,std::string> ( "mod=<["+::biobambam2::Licensing::formatNumber(getDefaultMod())+"]>", "print progress for every mod'th alignment if verbose" ) );
				V.push_back ( std::pair<std::string,std::string> ( "ranksplit=<["+::biobambam2::Licensing::formatNumber(getDefaultRankSplit())+"]>", "split rank pairs in names (see bam12split command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "rankstrip=<["+::biobambam2::Licensing::formatNumber(getDefaultRankStrip())+"]>", "strip ranks of names (see bam12strip command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "clipreinsert=<["+::biobambam2::Licensing::formatNumber(getDefaultClipReinsert())+"]>", "reinsert clipped sequence fragments (see bamclipreinsert command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "zztoname=<["+::biobambam2::Licensing::formatNumber(getDefaultZZToName())+"]>", "move rank from zz to name aux field to name (see bamzztoname command)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam2::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "index=<["+::biobambam2::Licensing::formatNumber(getDefaultIndex())+"]>", "create BAM index (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "indexfilename=<filename>", "file name for BAM index file (default: extend output file name)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=<filename>", "prefix for temporary files, default: create files in current directory" ) );
    	    	    	    	V.push_back ( std::pair<std::string,std::string> ( "sanity=<["+::biobambam2::Licensing::formatNumber(getDefaultSanity())+"]>", "extra checking of reads" ) );
				
				::biobambam2::Licensing::printMap(std::cerr,V);

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

