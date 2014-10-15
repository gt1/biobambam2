/**
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
**/
#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/timing/RealTimeClock.hpp>
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/regex/PosixRegex.hpp>

static std::string getDefaultInputFormat()
{
	return "bam";
}

void bamalignfrac(::libmaus::util::ArgInfo const & arginfo)
{
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
		libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	::libmaus::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
	::libmaus::bambam::BamAlignmentDecoder & dec = *ppdec;
	::libmaus::bambam::BamAlignment const & algn = dec.getAlignment();
        libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;

        uint64_t basealgn = 0;
        uint64_t clip = 0;
        uint64_t totalbases = 0;

        #if defined(LIBMAUS_HAVE_REGEX_H)
        std::string const regexs = arginfo.getUnparsedValue("name","");
        libmaus::util::unique_ptr<libmaus::regex::PosixRegex>::type regex_ptr;
        if ( regexs.size() )
	{
	        libmaus::util::unique_ptr<libmaus::regex::PosixRegex>::type tregex_ptr(new libmaus::regex::PosixRegex(regexs));
	        regex_ptr = UNIQUE_PTR_MOVE(tregex_ptr);
	}
	#endif

	while ( dec.readAlignment() )
	{
		if ( 
			algn.isMapped()
			#if defined(LIBMAUS_HAVE_REGEX_H)
			&&
			(
				(!regex_ptr)
				||
				(regex_ptr->findFirstMatch(algn.getName()) != -1)
			)
			#endif
		)
	        {
		        uint32_t const numcig = algn.getCigarOperations(cigop);
		        
		        totalbases += algn.getLseq();
		        
		        for ( uint64_t i = 0; i < numcig; ++i )
		        {
		        	switch ( cigop[i].first )
		        	{
		        		case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CMATCH:
		        		case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CINS:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CEQUAL:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDIFF:
						basealgn += cigop[i].second;
						break;
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CSOFT_CLIP:
						clip += cigop[i].second;
						break;
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CHARD_CLIP:
						totalbases += cigop[i].second;
						clip += cigop[i].second;
						break;
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CDEL:
					case libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_CREF_SKIP:
						break;
		        	}
		        }
	        }                                                                        
	}
	
	std::cerr << "total bases in mapped reads\t" << totalbases << std::endl;
	std::cerr << "clipped (hard and soft) bases in mapped reads\t" << clip << std::endl;
	std::cerr << "aligned bases in mapped reads\t" << basealgn << std::endl;
}

int main(int argc, char *argv[])
{
	try
	{
		libmaus::timing::RealTimeClock rtc; rtc.start();
		
		::libmaus::util::ArgInfo arginfo(argc,argv);
		
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
				std::cerr << ::biobambam::Licensing::license() << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;
								
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (default: read file from standard input)" ) );
				#if defined(BIOBAMBAM_LIBMAUS_HAVE_IO_LIB)
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", "input format: cram, bam or sam" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<[]>", "name of reference FastA in case of inputformat=cram" ) );
				#else
				V.push_back ( std::pair<std::string,std::string> ( "inputformat=<[bam]>", "input format: bam" ) );
				#endif
				
				#if defined(LIBMAUS_HAVE_REGEX_H)
				V.push_back ( std::pair<std::string,std::string> ( "name=<[]>", "consider only reads with names matching the given regualr expression (default: use all reads)" ) );
				#endif
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;

				return EXIT_SUCCESS;
			}
		
		bamalignfrac(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;		
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
