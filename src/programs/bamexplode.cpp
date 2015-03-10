/*
    libmaus
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
*/
#include <libmaus/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }
int getDefaultVerbose() { return 0; }
std::string getDefaultInputFormat() { return "bam"; }
uint64_t getDefaultSizeThres() { return 32*1024*1024; }
std::string getDefaultPrefix() { return "split_"; }

int bamexplode(libmaus::util::ArgInfo const & arginfo)
{
	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type Preader(libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));

	libmaus::bambam::BamBlockWriterBase::unique_ptr_type Pwriter;
	
	libmaus::bambam::BamAlignmentDecoder & decoder = Preader->getDecoder();
	libmaus::bambam::BamHeader const & header = decoder.getHeader();
	libmaus::bambam::BamAlignment const & algn = decoder.getAlignment();
	uint64_t nextfn = 0;
	uint64_t written = std::numeric_limits<uint64_t>::max();
	int32_t prevrefid = std::numeric_limits<int32_t>::max();
	std::string const outputformat = arginfo.getUnparsedValue("outputformat",libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat());
	std::string const prefix = arginfo.getUnparsedValue("prefix",getDefaultPrefix());
	uint64_t const thres = arginfo.getValueUnsignedNumeric("sizethres",getDefaultSizeThres());
	
	while ( decoder.readAlignment() )
	{
		int32_t const refid = algn.getRefID();

		if ( refid != prevrefid && written > thres )
		{
			Pwriter.reset();
			libmaus::util::ArgInfo argcopy(arginfo);
			std::ostringstream fnostr;
			fnostr << prefix << std::setw(6) << std::setfill('0') << nextfn++ << std::setw(0) << "." << outputformat;
			argcopy.replaceKey("O",fnostr.str());
			libmaus::bambam::BamBlockWriterBase::unique_ptr_type Twriter(libmaus::bambam::BamBlockWriterBaseFactory::construct(header,argcopy));
			Pwriter = UNIQUE_PTR_MOVE(Twriter);
			written = 0;
		}
		
		Pwriter->writeAlignment(algn);
		
		prevrefid = refid;
		written ++;
	}
	
	Pwriter.reset();

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
				V.push_back ( std::pair<std::string,std::string> ( std::string("inputformat=<[")+getDefaultInputFormat()+"]>", std::string("input format (") + libmaus::bambam::BamMultiAlignmentDecoderFactory::getValidInputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("outputformat=<[")+libmaus::bambam::BamBlockWriterBaseFactory::getDefaultOutputFormat()+"]>", std::string("output format (") + libmaus::bambam::BamBlockWriterBaseFactory::getValidOutputFormats() + ")" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[stdin]>", "input filename (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "inputthreads=<[1]>", "input helper threads (for inputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "reference=<>", "reference FastA (.fai file required, for cram i/o only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "range=<>", "coordinate range to be processed (for coordinate sorted indexed BAM input only)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "outputthreads=<[1]>", "output helper threads (for outputformat=bam only, default: 1)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "O=<[stdout]>", "output filename (standard output if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("prefix=<[")+getDefaultPrefix()+"]>", "prefix of output file names" ) );
				V.push_back ( std::pair<std::string,std::string> ( "thres=<["+::biobambam::Licensing::formatNumber(getDefaultSizeThres())+"]>", "size threshold for the creation of next file" ) );

				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				return EXIT_SUCCESS;
			}
			
		return bamexplode(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
