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
#include <config.h>

#include <libmaus/bambam/BamDecoder.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeader.hpp>
#include <libmaus/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus/util/ArgInfo.hpp>

#include <libmaus/timing/RealTimeClock.hpp>

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		bool const verbose = arginfo.getValue<unsigned int>("verbose",1);
		
		::libmaus::timing::RealTimeClock rtc; rtc.start();
		
		// gzip compression level for output
		int const level = arginfo.getValue<int>("level",1);
		
		::libmaus::bambam::BamDecoder bamfile(std::cin);
		std::string const headertext(bamfile.getHeader().text);
	
		// "reconstruct" command line
		std::string cl;
		for ( int i = 0; i < argc; ++i )
		{
			cl += argv[i];
			if ( i+1 < argc )
				cl += " ";
		}

		// add PG line to header
		std::string const upheadtext = ::libmaus::bambam::ProgramHeaderLineSet::addProgramLine(
			headertext,
			"bamfixmatecoordinatesnamesorted", // ID
			"bamfixmatecoordinatesnamesorted", // PN
			cl, // CL
			::libmaus::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
			std::string(PACKAGE_VERSION) // VN			
		);
		// construct new header
		::libmaus::bambam::BamHeader uphead(upheadtext);
		
		if ( uphead.getSortOrder() != "queryname" )
			uphead.changeSortOrder("unknown");
			
		std::string const & finalheadtext = uphead.text;
		::libmaus::bambam::BamHeader finalheader(finalheadtext);

		::libmaus::bambam::BamWriter writer(std::cout,finalheader,level);
		std::pair< std::pair< ::libmaus::bambam::BamAlignment::shared_ptr_type, bool> , std::pair< ::libmaus::bambam::BamAlignment::shared_ptr_type, bool> > 
			P(std::pair< ::libmaus::bambam::BamAlignment::shared_ptr_type, bool>(::libmaus::bambam::BamAlignment::shared_ptr_type(),false),std::pair< ::libmaus::bambam::BamAlignment::shared_ptr_type, bool>(::libmaus::bambam::BamAlignment::shared_ptr_type(),false));
		
		// try to read two alignments	
		P.first.second  = bamfile.readAlignment();
		if ( P.first.second )
		{
			P.first.first   = bamfile.salignment();
			P.second.second = P.first.second && bamfile.readAlignment();
			P.second.first  = bamfile.salignment();
		}
		
		uint64_t single = 0, pairs = 0;
		uint64_t proc = 0;
		uint64_t lastproc = 0;
		uint64_t const mod = 1024*1024;
		
		// while we have two alignments
		while ( P.first.second && P.second.second )
		{
			uint32_t const aflags = P.first.first->getFlags();
			uint32_t const bflags = P.second.first->getFlags();
		
			// same name?
			if ( 
				(aflags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED)
				&&
				(bflags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED)
				&&
				(! strcmp(P.first.first->getName(),P.second.first->getName()))
			)
			{			
				unsigned int const amap = (aflags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP) ? 0 : 1;
				unsigned int const bmap = (bflags & ::libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP) ? 0 : 1;

				// std::cerr << "Pair " << bam1_qname(P.first.first->alignment) << " amap=" << amap << " bmap=" << bmap << std::endl;
				
				// if exactly one of the two is mapped
				if ( amap + bmap == 1 )
				{
					::libmaus::bambam::BamAlignment::shared_ptr_type mapped = amap ? P.first.first : P.second.first;
					int64_t const tid = mapped->getRefID();
					int64_t const pos = mapped->getPos();
					
					// std::cerr << "tid=" << tid << " pos=" << pos << std::endl;
					
					// set all tid and pos values
					P.first.first->putRefId(tid);
					P.first.first->putPos(pos);
					P.first.first->putNextRefId(tid);
					P.first.first->putNextPos(pos);
					P.second.first->putRefId(tid);
					P.second.first->putPos(pos);
					P.second.first->putNextRefId(tid);
					P.second.first->putNextPos(pos);
				}
			
				// write alignments
				P.first.first->serialise(writer.getStream());
				P.second.first->serialise(writer.getStream());
				// read new alignments
				P.first.second = bamfile.readAlignment();
				if ( P.first.second )
				{
					P.first.first = bamfile.salignment();
					P.second.second = bamfile.readAlignment();
					P.second.first = bamfile.salignment();
				}
				
				pairs++;
				proc += 2;
			}
			// different names
			else
			{
				// write first alignment
				P.first.first->serialise(writer.getStream());
				// move second to first
				std::swap(P.first,P.second);
				// read new second
				P.second.second = P.first.second && bamfile.readAlignment();
				if ( P.second.second )
					P.second.first = bamfile.salignment();
				
				single++;
				proc += 1;
			}
			
			if ( verbose && (proc/mod != lastproc/mod) )
			{
				std::cerr << proc << "\t" << single << "\t" << pairs << "\t" <<
					proc/rtc.getElapsedSeconds() << "al/s"
					<< std::endl;
				lastproc = proc;
			}
		}
		
		if ( P.first.second )
		{
			P.first.first->serialise(writer.getStream());
			single++;
			proc += 1;
		}

		if ( verbose )
			std::cerr << proc << "\t" << single << "\t" << pairs << "\t" <<
				proc/rtc.getElapsedSeconds() << "al/s"
				<< std::endl;
			
		assert ( ! P.second.second );
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
