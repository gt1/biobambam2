/*
    biobambam
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
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/util/SimpleCountingHash.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 1; }
static std::string getDefaultInputFormat() { return "bam"; }

int bammapdist(libmaus::util::ArgInfo const & arginfo)
{
	bool const verbose = arginfo.getValue<unsigned int>("verbose",getDefaultVerbose());
	// input decoder wrapper
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(
			arginfo,false // put rank
		)
	);
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamHeader const & header = dec.getHeader();
	
	// getNumRef();
	::libmaus::bambam::BamAlignment const & algn = dec.getAlignment();

	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> pairmap(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> pairmap5(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> properpairmap(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> singlemap(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> partial1map(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> partial2map(8);
	double const loadthres = 0.9;
	uint64_t pairsunmapped = 0;
	uint64_t singleunmapped = 0;
	uint64_t numpairs = 0;
	uint64_t numsingle = 0;
	uint64_t c = 0;
	bool warningprinted = false;
	uint64_t skipped = 0;
	
	while ( dec.readAlignment() )
	{
		if ( algn.isMapped() && algn.getRefID() < 0 )
		{
			if ( ! warningprinted )
			{
				std::ostringstream estr;
				estr << "[D] Alignment has read marked as mapped but read coordinates are undefined:\n";
				estr << "[D] " << algn.formatAlignment(header) << std::endl;
				estr << "[D] This file is broken, skipping over broken alignments, stats will be incorrect.\n";
				std::cerr << estr.str();
				warningprinted = true;
			}
			skipped += 1;
			continue;
		}
		
		if ( algn.isPaired() && (algn.isMateMapped() && algn.getNextRefID() < 0) )
		{
			if ( ! warningprinted )
			{
				std::ostringstream estr;
				estr << "Alignment has mate marked as mapped but mate coordinates are undefined.\n";
				estr << "[D] " << algn.formatAlignment(header) << std::endl;
				estr << "[D] This file is broken, skipping over broken alignments, stats will be incorrect.\n";
				std::cerr << estr.str();
				warningprinted = true;
			}
			skipped += 1;
			continue;
		}

		if ( algn.isPaired() )
		{
			if ( 
				algn.isRead1() 
				&&
				(!algn.isSupplementary())
				&&
				(!algn.isSecondary())
			)
			{
				if ( algn.isMapped() && algn.isMateMapped() )
				{
					pairmap.insertExtend(
						(static_cast<uint64_t>(algn.getRefID()) << 32)
						|
						(static_cast<uint64_t>(algn.getNextRefID()) << 0),
						1,
						loadthres
					);
					
					if ( algn.getMapQ() >= 5 )
					{
						pairmap5.insertExtend(
							(static_cast<uint64_t>(algn.getRefID()) << 32)
							|
							(static_cast<uint64_t>(algn.getNextRefID()) << 0),
							1,
							loadthres
						);
					}
					
					if ( algn.isProper() )
						properpairmap.insertExtend(algn.getRefID(),1,loadthres);
				}
				else if ( algn.isMapped() )
				{
					partial1map.insertExtend(algn.getRefID(),1,loadthres);
				}
				else if ( algn.isMateMapped() )
				{
					partial2map.insertExtend(algn.getNextRefID(),1,loadthres);				
				}
				else
				{
					pairsunmapped++;
				}
				
				numpairs++;
			}
		}
		// single
		else
		{

			if ( 
				(!algn.isSupplementary())
				&&
				(!algn.isSecondary())			
			)
			{
				if ( 
					algn.isMapped() 
				)
				{
					singlemap.insertExtend(algn.getRefID(),1,loadthres);
				}
				else
				{
					singleunmapped++;
				}
				
				numsingle++;
			}
		}
		
		if ( verbose && ((++c & (1024*1024-1)) == 0) )
		{
			std::cerr << "[V] " << c 
				<< " "
				<< numpairs
				<< " "
				<< numsingle
				<< " "
				<< libmaus::util::MemUsage()
				<< std::endl;	
		}
	}
	
	std::cout << "number of pairs\t" << numpairs << std::endl;
	std::cout << "number of unmapped pairs\t" << pairsunmapped << std::endl;
	std::cout << "number of single\t" << numsingle << std::endl;
	std::cout << "number of unmapped single\t" << singleunmapped << std::endl;
	
	if ( skipped )
		std::cout << "number of skipped broken alignments\t" << skipped << std::endl;

	uint64_t split = 0;
	uint64_t split5 = 0;
	uint64_t mappedpairs = 0;
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> splitfirst(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> splitsecond(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> splitfirst5(8);
	::libmaus::util::SimpleCountingHash<uint64_t,uint64_t> splitsecond5(8);
	for ( uint64_t const * k = pairmap.begin(); k != pairmap.end(); ++k )
		if ( *k != pairmap.unused() )
		{
			if ( (*k >> 32) != ((*k) & 0xFFFFFFFFull) )
			{
				split += pairmap.getCount(*k);
				splitfirst.insertExtend((*k >> 32),pairmap.getCount(*k),loadthres);
				splitsecond.insertExtend((*k & 0xFFFFFFFFul),pairmap.getCount(*k),loadthres);
			}
			
			mappedpairs += pairmap.getCount(*k);
		}
	for ( uint64_t const * k = pairmap5.begin(); k != pairmap5.end(); ++k )
		if ( *k != pairmap5.unused() )
		{
			if ( (*k >> 32) != ((*k) & 0xFFFFFFFFull) )
			{
				split5 += pairmap.getCount(*k);
				splitfirst5.insertExtend((*k >> 32),pairmap.getCount(*k),loadthres);
				splitsecond5.insertExtend((*k & 0xFFFFFFFFul),pairmap.getCount(*k),loadthres);
			}
		}

	if ( mappedpairs )
	{
		std::cout << "total number of mapped pairs\t" << mappedpairs << std::endl;

		if ( split )
		{
			std::cout << "number of split reads\t" << split << std::endl;
			std::cout << "fraction of split reads in all mapped pairs\t" << static_cast<double>(split)/mappedpairs << std::endl;
		}
		if ( split5 )
		{
			std::cout << "number of split reads (MapQ>=5)\t" << split5 << std::endl;
			std::cout << "fraction of split reads (MapQ>=5) in all mapped pairs\t" << static_cast<double>(split5)/mappedpairs << std::endl;
		}
	}

	for ( uint64_t r = 0; r < header.getNumRef(); ++r )
	{
		std::ostringstream ostr;
		if ( singlemap.contains(r) )
			ostr << "\tsingle\t" << singlemap.getCount(r) << std::endl;

		if ( pairmap.contains( (r << 32) | r ) )
		{
			ostr << "\tmapped pairs\t" << pairmap.getCount((r << 32) | r) << std::endl;
			if ( properpairmap.contains(r) )
				ostr << "\tproper pairs\t" << properpairmap.getCount(r) << std::endl;
		}
		if ( partial1map.contains(r) )
			ostr << "\tpartial first\t" << partial1map.getCount(r) << std::endl;
		if ( partial2map.contains(r) )
			ostr << "\tpartial second\t" << partial2map.getCount(r) << std::endl;
			
		uint64_t splitref = 0;
		if ( splitfirst.contains(r) )
		{
			ostr << "\tsplit first\t" << splitfirst.getCount(r) << std::endl;
			ostr << "\tsplit first fraction in all\t" << static_cast<double>(splitfirst.getCount(r)) / split << std::endl;
			splitref += splitfirst.getCount(r);
		}
		if ( splitsecond.contains(r) )
		{
			ostr << "\tsplit second\t" << splitsecond.getCount(r) << std::endl;
			ostr << "\tsplit second fraction in all\t" << static_cast<double>(splitsecond.getCount(r)) / split << std::endl;
			splitref += splitsecond.getCount(r);
		}
		if ( splitref )
			ostr << "\tsplit fraction in all\t" << static_cast<double>(splitref)/split << std::endl;

		uint64_t splitref5 = 0;
		if ( splitfirst5.contains(r) )
		{
			ostr << "\tsplit (MapQ>=5) first\t" << splitfirst5.getCount(r) << std::endl;
			ostr << "\tsplit (MapQ>=5) first fraction in all\t" << static_cast<double>(splitfirst5.getCount(r)) / split5 << std::endl;
			splitref5 += splitfirst5.getCount(r);
		}
		if ( splitsecond5.contains(r) )
		{
			ostr << "\tsplit (MapQ>=5) second\t" << splitsecond5.getCount(r) << std::endl;
			ostr << "\tsplit (MapQ>=5) second fraction in all\t" << static_cast<double>(splitsecond5.getCount(r)) / split5 << std::endl;
			splitref5 += splitsecond5.getCount(r);
		}
		if ( splitref5 )
			ostr << "\tsplit (MapQ>=5) fraction in all\t" << static_cast<double>(splitref5)/split5 << std::endl;

		if ( ostr.str().size() )
		{
			std::cout << "reference sequence\t" << header.getRefIDName(r) << std::endl;
			std::cout << ostr.str();	
		}
	}
	
	return skipped ? EXIT_FAILURE : EXIT_SUCCESS;	
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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bammapdist(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;	
	}
}
