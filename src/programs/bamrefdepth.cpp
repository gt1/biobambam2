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

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 0; }

int bamrefdepth(libmaus::util::ArgInfo const & arginfo)
{
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());

	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type pdec(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	libmaus::bambam::BamAlignmentDecoder & bamdec = pdec->getDecoder();
	libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	std::string const sortorder = libmaus::bambam::BamHeader::getSortOrderStatic(header.text);
	
	libmaus::bambam::BamAlignment prevalgn;
	bool hasprev = false;
	uint64_t c = 0;
	
	std::deque<uint64_t> Q;
	uint64_t leftpos = 0;
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	libmaus::autoarray::AutoArray<char> decread;
	
	std::vector < std::string > refnames;
	for ( uint64_t i = 0; i < header.getNumRef(); ++i )
		refnames.push_back(header.getRefIDName(i));
                        	
	while ( bamdec.readAlignment() )
	{
		bool const ok =
			(!hasprev)
			||
			(
				(static_cast<uint32_t>(    algn.getRefID()) >
				 static_cast<uint32_t>(prevalgn.getRefID())
				)
				||
				(
					(static_cast<uint32_t>(    algn.getRefID()) ==
					 static_cast<uint32_t>(prevalgn.getRefID())
					)
					&&
					(static_cast<uint32_t>(    algn.getPos()) >=
					 static_cast<uint32_t>(prevalgn.getPos())
					)
				)
			);
			
		if ( ! ok )
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "File is not ordered by coordinate:";
			se.getStream() << prevalgn.formatAlignment(header) << std::endl;
			se.getStream() <<     algn.formatAlignment(header) << std::endl;
			se.finish();
			throw se;
		}
		
		// next reference sequence
		if ( hasprev && (algn.getRefID() != prevalgn.getRefID()) )
		{
			while ( Q.size() )
			{
				if ( Q[0] )
					std::cout << refnames[prevalgn.getRefID()] << "\t" << leftpos << "\t" << Q[0] << std::endl;
				Q.pop_front();
				leftpos++;				
			}
		
			// Q.resize(0);
			leftpos = 0;
		}
		
		if ( algn.isMapped() )
		{
			uint32_t const numcigop = algn.getCigarOperations(cigop);
			uint64_t const readlen = algn.decodeRead(decread);
			int64_t pos = algn.getPos();
			uint64_t readpos = 0;
			
			// std::cerr << "Q.size()=" << Q.size() << std::endl;

			while ( Q.size() && pos > static_cast<int64_t>(leftpos) )
			{
				if ( Q[0] )
					std::cout << refnames[algn.getRefID()] << "\t" << leftpos << "\t" << Q[0] << std::endl;
				Q.pop_front();
				leftpos++;
			}
			
			// skip soft clipping at front
			uint64_t cidx = 0;
			bool frontskip = true;
			while ( cidx < numcigop && frontskip )
				switch ( cigop[cidx].first )
				{
					// padding, ref skip, hard clip, deletion
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CPAD:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL:
						cidx += 1;
						break;
					// insertion/soft clipping, advance on read
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
						readpos += cigop[cidx++].second;
						break;
					// match/mismatch
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
						frontskip = false;
						break;
				}

			for ( ; cidx < numcigop ; ++cidx )
				switch ( cigop[cidx].first )
				{
					// padding, hard clipping (ignore)
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CPAD:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP:
						break;
					// ref skip, deletion (advance on reference)
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL:
						pos += cigop[cidx].second;
						break;
					// insertion/soft clipping (advance on read)
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
						readpos += cigop[cidx].second;
						break;
					// match/mismatch
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
						for ( uint64_t i = 0; static_cast<int64_t>(i) < cigop[cidx].second; ++i, ++pos, ++readpos )
						{
							if ( ! Q.size() )
							{
								Q.push_back(1);
								leftpos = pos;
							}
							else
							{
								while ( pos-leftpos >= Q.size() )
									Q.push_back(0);
							
								assert ( pos-leftpos < Q.size() );		
								Q[pos-leftpos] += 1;
							}
						}
						break;
				}
			
			#if 0
			if ( readpos != readlen )
			{
				std::cerr << "readpos=" << readpos << " readlen=" << readlen << std::endl;			
				std::cerr << algn.formatAlignment(header) << std::endl;
			}	
			#endif
			assert (readpos == readlen);
		}
			
		prevalgn.swap(algn);
		hasprev = true;

		if ( verbose && ( ((++c) & ((1ull<<20)-1)) == 0 ) )
			std::cerr << "[V] " << c << std::endl;
	}

	while ( Q.size() )
	{
		if ( Q[0] )
			std::cout << refnames[prevalgn.getRefID()] << "\t" << leftpos << "\t" << Q[0] << std::endl;
		Q.pop_front();
		leftpos++;				
	}
		
	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

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
			
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<["+::biobambam::Licensing::formatNumber(getDefaultVerbose())+"]>", "print progress report" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamrefdepth(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}

