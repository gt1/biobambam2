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
#include <libmaus/bambam/CollatingBamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>

#include <config.h>
#include <biobambam/Licensing.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
static unsigned int getDefaultColHashBits() { return 20; }
static uint64_t getDefaultColListSize() { return 512*1024; }
static unsigned int getDefaultPairsOnly() { return 0; }

int bamCollate(::libmaus::util::ArgInfo const & arginfo)
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
	
	std::string const tmpfile = arginfo.getValue<std::string>("tmpfile",arginfo.getDefaultTmpFileName());
	std::string const readgroups = arginfo.getValue<std::string>("readgroups",std::string());
	bool const pairsonly = arginfo.getValue<unsigned int>("pairsonly",false);
	::libmaus::util::TempFileRemovalContainer::setup();
	::libmaus::util::TempFileRemovalContainer::addTempFile(tmpfile);
	std::vector < std::string > vreadgroups;
	::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type LHTsnofailure;
	
	if ( readgroups.size() )
	{
		std::deque<std::string> qreadgroups = ::libmaus::util::stringFunctions::tokenize(readgroups,std::string(","));
		vreadgroups = std::vector<std::string>(qreadgroups.begin(),qreadgroups.end());
		::libmaus::trie::Trie<char> trienofailure;
		trienofailure.insertContainer(vreadgroups);
		::libmaus::trie::LinearHashTrie<char,uint32_t>::unique_ptr_type LHTnofailure = UNIQUE_PTR_MOVE(trienofailure.toLinearHashTrie<uint32_t>());
		LHTsnofailure = ::libmaus::trie::LinearHashTrie<char,uint32_t>::shared_ptr_type(LHTnofailure.release());
	}
	
	unsigned int const colhashbits = arginfo.getValue<unsigned int>("colhashbits",getDefaultColHashBits());
	unsigned int const collistsize = arginfo.getValue<unsigned int>("collistsize",getDefaultColListSize());
	int const level = arginfo.getValue<int>("level",getDefaultLevel());
	
	switch ( level )
	{
		case Z_NO_COMPRESSION:
		case Z_BEST_SPEED:
		case Z_BEST_COMPRESSION:
		case Z_DEFAULT_COMPRESSION:
			break;
		default:
		{
			::libmaus::exception::LibMausException se;
			se.getStream()
				<< "Unknown compression level, please use"
				<< " level=" << Z_DEFAULT_COMPRESSION << " (default) or"
				<< " level=" << Z_BEST_SPEED << " (fast) or"
				<< " level=" << Z_BEST_COMPRESSION << " (best) or"
				<< " level=" << Z_NO_COMPRESSION << " (no compression)" << std::endl;
			se.finish();
			throw se;
		}
			break;
	}

	::libmaus::bambam::CollatingBamDecoder CBD(std::cin,tmpfile,false/* add rank */,colhashbits,collistsize);
	::libmaus::bambam::BamHeader const & bamheader = CBD.bamdecoder.bamheader;

	std::string const headertext(bamheader.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bamcollate", // ID
		"bamcollate", // PN
		arginfo.commandline, // CL
		::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN			
	);
	// construct new header
	::libmaus::bambam::BamHeader uphead(upheadtext);
	// std::string const insortorder = uphead.getSortOrder();
	uphead.changeSortOrder("unknown");

	::libmaus::bambam::BamWriter writer(std::cout,uphead);
	// bool tryPair(std::pair<alignment_ptr_type,alignment_ptr_type> & P)
	// std::pair<alignment_ptr_type,alignment_ptr_type> P;
	typedef ::libmaus::bambam::CollatingBamDecoder::alignment_ptr_type alignment_ptr_type;
		
	alignment_ptr_type a;
	if ( pairsonly )
	{
		while ( a = CBD.getPair() )
		{
			alignment_ptr_type b = CBD.getPair();
			assert ( b );
			
			if ( !readgroups.size() 
				||
				( 
					LHTsnofailure->searchCompleteNoFailure(std::string(a->getReadGroup())) != -1 &&
					LHTsnofailure->searchCompleteNoFailure(std::string(b->getReadGroup())) != -1
				)	
			)
			{
				a->serialise(writer.bgzfos);
				b->serialise(writer.bgzfos);
			}
		}
	}
	else
	{
		while ( a = CBD.get() )
			if ( !readgroups.size() || LHTsnofailure->searchCompleteNoFailure(std::string(a->getReadGroup())) != -1 )
				a->serialise(writer.bgzfos);
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
				
				V.push_back ( std::pair<std::string,std::string> ( "tmpfile=[<filename>]", "prefix for temporary files, default: create files in current directory" ) );
				V.push_back ( std::pair<std::string,std::string> ( "readgroups=[]", "filter for read groups, default: do not filter" ) );
				V.push_back ( std::pair<std::string,std::string> ( "pairsonly=["+ ::biobambam::Licensing::formatNumber(getDefaultPairsOnly())+ "]", "output complete pairs only (1=yes,0=no)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "level=<["+::biobambam::Licensing::formatNumber(getDefaultLevel())+"]>", "compression settings for output bam file (0=uncompressed,1=fast,9=best,-1=zlib default)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "colhashbits=<["+::biobambam::Licensing::formatNumber(getDefaultColHashBits())+"]>", "log_2 of size of hash table used for collation" ) );
				V.push_back ( std::pair<std::string,std::string> ( "collistsize=<["+::biobambam::Licensing::formatNumber(getDefaultColListSize())+"]>", "output list size for collation" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamCollate(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

