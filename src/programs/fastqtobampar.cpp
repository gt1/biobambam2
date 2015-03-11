/**
    biobambam
    Copyright (C) 2009-2015 German Tischler
    Copyright (C) 2011-2015 Genome Research Limited

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

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/aio/PosixFdOutputStream.hpp>
#include <libmaus/bambam/RgInfo.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/parallel/FastqToBamControl.hpp>
#include <libmaus/lz/BgzfDeflate.hpp>
#include <libmaus/parallel/NumCpus.hpp>
#include <libmaus/util/MemUsage.hpp>

static std::string writeHeader(libmaus::util::ArgInfo const & arginfo, std::ostream & out)
{
	libmaus::bambam::RgInfo const rginfo(arginfo);

	std::ostringstream headerostr;
	headerostr << "@HD\tVN:1.4\tSO:unknown\n";
	headerostr 
		<< "@PG"<< "\t" 
		<< "ID:" << "fastqtobam" << "\t" 
		<< "PN:" << "fastqtobam" << "\t"
		<< "CL:" << arginfo.commandline << "\t"
		<< "VN:" << std::string(PACKAGE_VERSION)
		<< std::endl;
	headerostr << rginfo.toString();
	::libmaus::bambam::BamHeader bamheader;
	bamheader.text = headerostr.str();		

	libmaus::lz::BgzfOutputStream bgzf(out);
	bamheader.serialise(bgzf);
	bgzf.flush();
	
	return rginfo.ID;
}

static int fastqtobampar(libmaus::util::ArgInfo const & arginfo)
{
	std::ostream & out = std::cout;
	uint64_t const numlogcpus = arginfo.getValue<int>("threads", 1 /* libmaus::parallel::NumCpus::getNumLogicalProcessors() */);
	int const level = arginfo.getValue<int>("level",Z_DEFAULT_COMPRESSION);
		
	std::string const rgid = writeHeader(arginfo,out);

	libmaus::parallel::SimpleThreadPool STP(numlogcpus);
	uint64_t const outblocks = 1024;
	uint64_t const inputblocksize = 1024*1024*64;
	uint64_t const inblocks = 64;
	libmaus::aio::PosixFdInputStream PFIS(STDIN_FILENO);
	libmaus::bambam::parallel::FastqToBamControl FTBC(PFIS,out,STP,inblocks,inputblocksize,outblocks,level,rgid);

	FTBC.enqueReadPackage();
	FTBC.waitCompressionFinished();
		
	STP.terminate();		
	STP.join();

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::timing::RealTimeClock rtc; rtc.start();	
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
				std::cerr << ::biobambam::Licensing::license() << std::endl;
				std::cerr << "Key=Value pairs:" << std::endl;
				std::cerr << std::endl;
				
				std::vector< std::pair<std::string,std::string> > V;

				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
				#if 0
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<[0]>", "print progress report (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[input file name]>", "input file names (standard input if unset)" ) );
				#endif
				V.push_back ( std::pair<std::string,std::string> ( "threads=<[1]>", "number of threads used (default: 1)" ) );

				V.push_back ( std::pair<std::string,std::string> ( "RGID=<>", "read group id for reads (default: do not set a read group id)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGCN=<>", "CN field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGDS=<>", "DS field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGDT=<>", "DT field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGFO=<>", "FO field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGKS=<>", "KS field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGLB=<>", "LB field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGPG=<fastqtobam>", "CN field of RG header line if RGID is set (default: fastqtobam)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGPI=<>", "PI field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGPL=<>", "PL field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGPU=<>", "PU field of RG header line if RGID is set (default: not present)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "RGSM=<>", "SM field of RG header line if RGID is set (default: not present)" ) );
				
				#if 0
				V.push_back ( std::pair<std::string,std::string> ( "qualityoffset=<["+::biobambam::Licensing::formatNumber(getDefaultQualityOffset())+"]>", "quality offset (default: 33)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "qualitymax=<["+::biobambam::Licensing::formatNumber(getDefaultQualityMaximum())+"]>", "maximum valid quality value (default: 41)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "checkquality=<["+::biobambam::Licensing::formatNumber(getDefaultCheckQuality())+"]>", "check quality (default: true)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("namescheme=<[")+(getDefaultNameScheme())+"]>", "read name scheme (generic, c18s, c18pe, pairedfiles)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "qualityhist=<["+::biobambam::Licensing::formatNumber(getDefaultQualityHist())+"]>", "compute quality histogram and print it on standard error" ) );
				#endif
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				
				std::cerr << "The I key can be given twice for a pair of synced FastQ files." << std::endl;
				std::cerr << "Any none key=value arguments will be considered as input file names." << std::endl;

				return EXIT_SUCCESS;
			}
		
			
		fastqtobampar(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
