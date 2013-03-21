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

#if ! defined(BAMFILE_HPP)
#define BAMFILE_HPP

#include <libmaus/util/stringFunctions.hpp>
#include <libmaus/util/unique_ptr.hpp>
#include <libmaus/util/shared_ptr.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/exception/LibMausException.hpp>

#include <bambam/BamBamConfig.hpp>

#if defined(BAMBAM_HAVE_SAMTOOLS)

#if defined(HAVE_BAM_H)
#include <bam.h>
#include <sam.h>
#elif defined(HAVE_SAMTOOLS_BAM_H)
#include <samtools/bam.h>
#include <samtools/sam.h>
#else
#error "Required bam.h header not available."
#endif


namespace bambam
{
	struct BamAlignment
	{
		typedef BamAlignment this_type;
		typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
		typedef ::libmaus::util::shared_ptr<this_type>::type shared_ptr_type;
	
		bam1_t * alignment;
	
		BamAlignment()
		{
			alignment = bam_init1();

			if ( ! alignment )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "Failed to init alignment by bam_init1.";
				se.finish();
				throw se;
			}
		}
		~BamAlignment()
		{
			bam_destroy1(alignment);		
		}
	};	

	struct BamFile // : public BamFileBase
	{
		typedef BamFile this_type;
		typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
		typedef ::libmaus::util::shared_ptr<this_type>::type shared_ptr_type;
	
		samfile_t * bamfile;
		bam1_t * alignment;
		bam_header_t * header;
				
		bool readAlignment(bam1_t * align)
		{
			return (samread(bamfile,align) >= 0);
		}

		bool readAlignment(BamAlignment & align)
		{
			return (samread(bamfile,align.alignment) >= 0);
		}

		bam_header_t * getBamHeader()
		{
			return header;
		}
		
		std::string getHeaderAsString() const
		{
			if ( ! header || ! header->text || ! header->l_text )
				return std::string();
			else
				return std::string(header->text,header->text+header->l_text);
		}
		
		void init(std::string const inputformat, std::string const inputfilename)
		{
			
			if ( inputformat != "bam" && inputformat != "sam" )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "Cannot handle input format " << inputformat << std::endl;
				se.finish();
				throw se;	
			}
			
			std::string const mode = (inputformat == "bam") ? "rb" : "r";
			
			// std::cerr << "Mode " << mode << std::endl;
			
			// open bam file
			bamfile = samopen(inputfilename.c_str(),mode.c_str(),0);

			if ( ! bamfile )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "Failed to parse standard input as bam.";
				se.finish();
				throw se;
			}
			
			// read bam header
			header = bamfile->header;

			if ( ! header )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "Failed to parse standard input as bam (failed to read header).";
				se.finish();
				throw se;		
			}
			
			alignment = bam_init1();

			if ( ! alignment )
			{
				samclose(bamfile);
				::libmaus::exception::LibMausException se;
				se.getStream() << "Failed to init alignment.";
				se.finish();
				throw se;
			}
		
		}

		BamFile(::libmaus::util::ArgInfo const & arginfo)
		: bamfile(0), alignment(0), header(0)
		{
			std::string const inputformat = arginfo.getValue<std::string>("inputformat",std::string("bam"));	
			std::string const inputfilename = arginfo.getValue<std::string>("filename",std::string("-"));	
			init(inputformat,inputfilename);
		}
		
		BamFile(std::string const inputformat, std::string const inputfilename)
		: bamfile(0), alignment(0), header(0)
		{
			init(inputformat,inputfilename);
		}

		
		~BamFile()
		{
			bam_destroy1(alignment);
			samclose(bamfile);
		}
	};
}
#endif

#endif
