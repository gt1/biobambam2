/*
    biobambam2
    Copyright (C) 2016 German Tischler

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

#include <iostream>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>

#include <biobambam2/BamBamConfig.hpp>
#include <biobambam2/Licensing.hpp>

static std::string getDefaultPrefix()
{
	return "fasta_";
}

static uint64_t getDefaultLineLength()
{
	return 80;
}

int fastaexplode(libmaus2::util::ArgParser const & arg)
{
	libmaus2::fastx::StreamFastAReaderWrapper SFA(std::cin);
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

	uint64_t cnt = 0;
	std::string const prefixarg =
		arg.uniqueArgPresent("p") ?
			"p" : ( arg.uniqueArgPresent("prefix") ? "prefix" : "" );
	std::string const prefix = prefixarg.size() ? arg[prefixarg] : getDefaultPrefix();
	bool const singleline = arg.argPresent("s") || arg.argPresent("singleline");
	bool const longname = arg.argPresent("L") || arg.argPresent("longname");
	bool const dataonly = arg.argPresent("d") || arg.argPresent("dataonly");
	uint64_t const linelength = arg.uniqueArgPresent("c") ? arg.getUnsignedNumericArg<uint64_t>("c") : getDefaultLineLength();

	while ( SFA.getNextPatternUnlocked(pattern) )
	{
		std::ostringstream fnostr;
		fnostr << prefix << std::setw(6) << std::setfill('0') << cnt++ << std::setw(0) << (dataonly ? "" : ".fasta");
		std::string const fn = fnostr.str();
		libmaus2::aio::OutputStreamInstance OSI(fn);

		std::string & spat = pattern.spattern;

		for ( uint64_t i = 0; i < spat.size(); ++i )
			spat[i] = toupper(spat[i]);

		if ( !longname )
			pattern.sid = pattern.getShortStringId();

		if ( dataonly )
			OSI.write(pattern.spattern.c_str(),pattern.spattern.size());
		else if ( singleline )
			OSI << pattern;
		else
			pattern.printMultiLine(OSI,linelength);

		std::cout << fn << "\t" << pattern.getStringId() << std::endl;
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if (
			arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version")
		)
		{
			std::cerr << ::biobambam2::Licensing::license();
			return EXIT_SUCCESS;
		}
		else if (
			arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help")
		)
		{
			std::cerr << ::biobambam2::Licensing::license();
			std::cerr << std::endl;
			std::cerr << "options:" << std::endl;
			std::cerr << std::endl;

			std::vector< std::pair<std::string,std::string> > V;

			V.push_back ( std::pair<std::string,std::string> ( "-v/--version", "print version number and quit" ) );
			V.push_back ( std::pair<std::string,std::string> ( "-h/--help", "print help message and quit" ) );
			V.push_back ( std::pair<std::string,std::string> ( "-s/--singleline", "do not wrap sequence data lines" ) );
			V.push_back ( std::pair<std::string,std::string> ( "-L/--longname", "do not shorten name" ) );
			V.push_back ( std::pair<std::string,std::string> ( "-l<cols>", "line length (default: "+libmaus2::util::NumberSerialisation::formatNumber(getDefaultLineLength(),0)+")" ) );
			V.push_back ( std::pair<std::string,std::string> ( "-d/--dataonly", "do not print FastA header (data only)" ) );
			V.push_back ( std::pair<std::string,std::string> ( "-p", std::string("prefix for FastA output files (default: ")+getDefaultPrefix()+")" ) );

			::biobambam2::Licensing::printMap(std::cerr,V);

			std::cerr << std::endl;
			return EXIT_SUCCESS;

		}

		return fastaexplode(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
