/**
    biobambam
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
#include <libmaus/fastx/StreamFastQReader.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/util/MemUsage.hpp>

#include <libmaus/lz/BufferedGzipStream.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <iomanip>

#include <config.h>

int getDefaultQualityOffset()
{
	return 33;
}

bool getDefaultCheckQuality()
{
	return true;
}

int getDefaultMD5()
{
	return 0;
}

static int getDefaultLevel() 
{
	return Z_DEFAULT_COMPRESSION;
}
static int getDefaultVerbose() 
{
	return 0;
}

static int getDefaultGz() 
{
	return 0;
}
static std::string getDefaultNameScheme() 
{
	return "generic";
}

static int getLevel(libmaus::util::ArgInfo const & arginfo)
{
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
	
	return level;
}

struct RgInfo
{
	std::string ID;
	std::string CN;
	std::string DS;
	std::string DT;
	std::string FO;
	std::string KS;
	std::string LB;
	std::string PG; // = fastqtobam
	std::string PI;
	std::string PL;
	std::string PU;
	std::string SM;
	
	RgInfo() {}
	RgInfo(libmaus::util::ArgInfo const & arginfo)
	:
		ID(arginfo.getUnparsedValue("RGID","")),
		CN(arginfo.getUnparsedValue("RGCN","")),
		DS(arginfo.getUnparsedValue("RGDS","")),
		DT(arginfo.getUnparsedValue("RGDT","")),
		FO(arginfo.getUnparsedValue("RGFO","")),
		KS(arginfo.getUnparsedValue("RGKS","")),
		LB(arginfo.getUnparsedValue("RGLB","")),
		PG(arginfo.getUnparsedValue("RGPG","fastqtobam")),
		PI(arginfo.getUnparsedValue("RGPI","")),
		PL(arginfo.getUnparsedValue("RGPL","")),
		PU(arginfo.getUnparsedValue("RGPU","")),
		SM(arginfo.getUnparsedValue("RGSM",""))
	{
		
	}
	
	std::string toString() const
	{
		std::ostringstream ostr;
		
		if ( ID.size() )
		{
			ostr << "@RG\tID:" << ID;
			
			if ( CN.size() ) ostr << "\tCN:" << CN;
			if ( DS.size() ) ostr << "\tDS:" << DS;
			if ( DT.size() ) ostr << "\tDT:" << DT;
			if ( FO.size() ) ostr << "\tFO:" << FO;
			if ( KS.size() ) ostr << "\tKS:" << KS;
			if ( LB.size() ) ostr << "\tLB:" << LB;
			if ( PG.size() ) ostr << "\tPG:" << PG;
			if ( PI.size() ) ostr << "\tPI:" << PI;
			if ( PL.size() ) ostr << "\tPL:" << PL;
			if ( PU.size() ) ostr << "\tPU:" << PU;
			if ( SM.size() ) ostr << "\tSM:" << SM;
			
			ostr << "\n";
		}
		
		return ostr.str();
	}
};

enum fastq_name_scheme_type { 
	fastq_name_scheme_generic,
	fastq_name_scheme_casava18_single,
	fastq_name_scheme_casava18_paired_end
};

static fastq_name_scheme_type parseNameScheme(std::string const & schemename)
{
	if ( schemename == "generic" )
		return fastq_name_scheme_generic;
	else if ( schemename == "c18s" )
		return fastq_name_scheme_casava18_single;
	else if ( schemename == "c18pe" )
		return fastq_name_scheme_casava18_paired_end;
	else
	{
		libmaus::exception::LibMausException lme;
		lme.getStream() << "unknown read name scheme " << schemename << std::endl;
		lme.finish();
		throw lme;	
	}
}

struct NameInfo
{
	std::string name;
	bool ispair;
	bool isfirst;
	uint64_t gl;
	uint64_t gr;
	fastq_name_scheme_type const namescheme;
		
	NameInfo() : name(), ispair(false), isfirst(false), gl(0), gr(0), namescheme(fastq_name_scheme_generic) {}
	NameInfo(
		std::string const & rname, 
		libmaus::fastx::SpaceTable const & ST,
		fastq_name_scheme_type const rnamescheme
	)
	: name(rname), ispair(false), isfirst(true), gl(0), gr(name.size()), namescheme(rnamescheme)
	{
		switch ( namescheme )
		{
			case fastq_name_scheme_generic:
			{
				uint64_t l = 0;
				// skip space at start
				while ( l != name.size() && ST.spacetable[name[l]] )
					++l;
				if ( l != name.size() )
				{
					// skip non space symbols
					uint64_t r = l;
					while ( (r != name.size()) && (!ST.spacetable[name[r]]) )
						++r;
					
					// check whether this part of the name ends on /1 or /2
					if ( r-l >= 2 && name[r-2] == '/' && (name[r-1] == '1' || name[r-1] == '2') )
					{
						if ( name[r-1] == '2' )
							isfirst = false;
							
						gl = l;
						gr = r;
						ispair = true;
					}
					else
					{
						gl = l;
						gr = r;
						ispair = false;
					}
					
					l = r;
				}	
				
				break;
			}
			case fastq_name_scheme_casava18_single:
			case fastq_name_scheme_casava18_paired_end:
			{
				uint64_t l = 0;
				
				while ( l != name.size() && ST.spacetable[name[l]] )
					++l;
					
				uint64_t const l0 = l;
				
				while ( l != name.size() && ! ST.spacetable[name[l]] )
					++l;
				
				uint64_t const l0l = l-l0;	

				while ( l != name.size() && ST.spacetable[name[l]] )
					++l;
					
				uint64_t const l1 = l;

				while ( l != name.size() && ! ST.spacetable[name[l]] )
					++l;
					
				uint64_t const l1l = l-l1;
			
				// count number of colons	
				uint64_t l0c = 0;
				for ( uint64_t i = l0; i != l0+l0l; ++i )
					l0c += (name[i] == ':');
				uint64_t l1c = 0;
				for ( uint64_t i = l1; i != l1+l1l; ++i )
					l1c += (name[i] == ':');
			
				if ( l0c != 6 || l1c != 3 )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "malformed read name " << name << " (wrong number of colon separated fields)" << std::endl;
					se.finish();
					throw se;
				}
				
				gl = l0;
				gr = l0+l0l;
				
				if ( namescheme == fastq_name_scheme_casava18_single )
				{
					ispair = false;
				}
				else
				{
					ispair = true;
					
					uint64_t fragid = 0;
					unsigned int fragidlen = 0;
					uint64_t p = l1;
					while ( isdigit(name[p]) )
					{
						fragid *= 10;
						fragid += name[p++]-'0';
						fragidlen++;
					}
					
					if ( (! fragidlen) || (fragid<1) || (fragid>2) || name[p] != ':' )
					{
						::libmaus::exception::LibMausException se;
						se.getStream() << "malformed read name " << name << " (malformed fragment id)" << std::endl;
						se.finish();
						throw se;	
					}
					
					isfirst = (fragid==1);
				}
				
				break;
			}
		}
	}
	
	std::string getName() const
	{
		switch ( namescheme )
		{
			case fastq_name_scheme_generic:
				if ( ispair )
					return name.substr(gl,gr-gl-2);		
				else
					return name.substr(gl,gr-gl);
				break;
			case fastq_name_scheme_casava18_single:
			case fastq_name_scheme_casava18_paired_end:
				return name.substr(gl,gr-gl);
			default:
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "NameInfo::getName(): invalid name scheme" << std::endl;
				se.finish();
				throw se;			
			}
		}
	}
};

template<typename element_type>
void checkFastQElement(element_type const & element, int const qualityoffset)
{
	for ( uint64_t i = 0; i < element.quality.size(); ++i )
		if ( 
			static_cast<int>(element.quality[i]) - qualityoffset < 0
			||
			static_cast<int>(element.quality[i]) - qualityoffset > 41
		)
	{
		libmaus::exception::LibMausException lme;
		lme.getStream() << "quality string " << element.quality << " of read " << element.sid << " is invalid" << std::endl;
		lme.finish();
		throw lme;
	}	
}

template<typename writer_type>
void fastqtobamSingle(
	std::istream & in, writer_type & bamwr, int const verbose, std::string const & rgid, 
	int const qualityoffset,
	bool const checkquality,
	fastq_name_scheme_type const namescheme
	)
{
	typedef ::libmaus::fastx::StreamFastQReaderWrapper reader_type;
	typedef reader_type::pattern_type pattern_type;
	::libmaus::fastx::StreamFastQReaderWrapper fqin(in);	
	pattern_type element;
	libmaus::fastx::SpaceTable ST;
	uint64_t proc = 0;
		
	while ( fqin.getNextPatternUnlocked(element) )
	{
		std::string const & name = element.sid;
		NameInfo const NI(name,ST,namescheme);
		
		if ( checkquality )
			checkFastQElement(element,qualityoffset);

		if ( NI.ispair )
		{
			bamwr.encodeAlignment(
				NI.getName(),
				-1,
				-1,
				0,
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
				(NI.isfirst ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
				std::string(),
				-1,
				-1,
				0,
				element.spattern,
				element.quality,
				qualityoffset
			);
			if ( rgid.size() )
				bamwr.putAuxString("RG",rgid.c_str());
			bamwr.commit();
		}
		else
		{
			bamwr.encodeAlignment(
				NI.getName(),
				-1,
				-1,
				0,
				libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP,
				std::string(),
				-1,
				-1,
				0,
				element.spattern,
				element.quality,
				qualityoffset
			);
			if ( rgid.size() )
				bamwr.putAuxString("RG",rgid.c_str());
			bamwr.commit();		
		}

		proc += 1;
		if ( verbose && ((proc & (1024*1024-1)) == 0) )
			std::cerr << "[V] " << proc << std::endl;
	}
	if ( verbose )
		std::cerr << "[V] " << proc << std::endl;
}

template<typename writer_type>
void fastqtobamPair(
	std::istream & in_1, std::istream & in_2, 
	writer_type & bamwr, int const verbose, std::string const & rgid, 
	int const qualityoffset, bool const checkquality,
	fastq_name_scheme_type const namescheme
)
{
	typedef ::libmaus::fastx::StreamFastQReaderWrapper reader_type;
	typedef reader_type::pattern_type pattern_type;
	::libmaus::fastx::StreamFastQReaderWrapper fqin_1(in_1);	
	::libmaus::fastx::StreamFastQReaderWrapper fqin_2(in_2);
	pattern_type element_1, element_2;
	libmaus::fastx::SpaceTable ST;
	uint64_t proc = 0;
	
	while ( fqin_1.getNextPatternUnlocked(element_1) )
	{
		bool const ok_2 = fqin_2.getNextPatternUnlocked(element_2);
		
		if ( ! ok_2 )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "fastq input files contain a different number of reads" << std::endl;
			ex.finish();
			throw ex;
		}

		if ( checkquality )
		{
			checkFastQElement(element_1,qualityoffset);
			checkFastQElement(element_2,qualityoffset);
		}

		std::string const & name_1 = element_1.sid;
		std::string const & name_2 = element_2.sid;
		NameInfo const NI_1(name_1,ST,namescheme);
		NameInfo const NI_2(name_2,ST,namescheme);
		
		if ( (!NI_1.ispair) || (!NI_1.isfirst) )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "name " << name_1 << " does not look like a first mate read name" << std::endl;
			ex.finish();
			throw ex;		
		}
		if ( NI_2.isfirst )
		{
			libmaus::exception::LibMausException ex;
			ex.getStream() << "name " << name_2 << " does not look like a second mate read name" << std::endl;
			ex.finish();
			throw ex;		
		}

		bamwr.encodeAlignment(
			NI_1.getName(),
			-1,
			-1,
			0,
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
			(NI_1.isfirst ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
			std::string(),
			-1,
			-1,
			0,
			element_1.spattern,
			element_1.quality,
			qualityoffset
		);
		if ( rgid.size() )
			bamwr.putAuxString("RG",rgid.c_str());
		bamwr.commit();

		bamwr.encodeAlignment(
			NI_2.getName(),
			-1,
			-1,
			0,
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
			libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
			(NI_2.isfirst ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
			std::string(),
			-1,
			-1,
			0,
			element_2.spattern,
			element_2.quality,
			qualityoffset
		);
		if ( rgid.size() )
			bamwr.putAuxString("RG",rgid.c_str());
		bamwr.commit();
		
		proc += 2;
		if ( verbose && ((proc & (1024*1024-1)) == 0) )
			std::cerr << "[V] " << proc << std::endl;
	}

	if ( verbose )
		std::cerr << "[V] " << proc << std::endl;
}


void fastqtobam(libmaus::util::ArgInfo const & arginfo)
{
	std::vector<std::string> filenames = arginfo.getPairValues("I");
	for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
		filenames.push_back(arginfo.restargs[i]);

	int const level = getLevel(arginfo);
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const gz = arginfo.getValue<int>("gz",getDefaultGz());
	unsigned int const threads = arginfo.getValue<unsigned int>("threads",1);
	int const qualityoffset = arginfo.getValue<int>("qualityoffset",getDefaultQualityOffset());
	bool const checkquality = arginfo.getValue<int>("checkquality",getDefaultCheckQuality());
	fastq_name_scheme_type const namescheme = parseNameScheme(
		arginfo.getValue<std::string>("namescheme",getDefaultNameScheme())
	);
	RgInfo const rginfo(arginfo);
	std::string const rgid = rginfo.ID; // ;arginfo.getUnparsedValue("RG","");

	::libmaus::bambam::BamHeader bamheader;
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

	bamheader.text = headerostr.str();		

	/*
	 * start md5 callbacks
	 */
	std::string md5filename;

	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > cbs;
	::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Pmd5cb;
	if ( arginfo.getValue<unsigned int>("md5",getDefaultMD5()) )
	{
		if ( arginfo.hasArg("md5filename") &&  arginfo.getUnparsedValue("md5filename","") != "" )
			md5filename = arginfo.getUnparsedValue("md5filename","");
		else
			std::cerr << "[V] no filename for md5 given, not creating hash" << std::endl;

		if ( md5filename.size() )
		{
			::libmaus::lz::BgzfDeflateOutputCallbackMD5::unique_ptr_type Tmd5cb(new ::libmaus::lz::BgzfDeflateOutputCallbackMD5);
			Pmd5cb = UNIQUE_PTR_MOVE(Tmd5cb);
			cbs.push_back(Pmd5cb.get());
		}
	}
	std::vector< ::libmaus::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
	if ( cbs.size() )
		Pcbs = &cbs;
	/*
	 * end md5 callbacks
	 */

	if ( threads <= 1 )
	{
		::libmaus::bambam::BamWriter::unique_ptr_type bamwr(
			new ::libmaus::bambam::BamWriter(
				std::cout,bamheader,level,Pcbs
			)
		);

		if ( filenames.size() == 0 )
		{
			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS(std::cin);
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
			else
			{
				fastqtobamSingle(std::cin,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
		}
		else if ( filenames.size() == 1 )
		{
			libmaus::aio::CheckedInputStream CIS(filenames[0]);		
			
			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS(CIS);
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
			else
			{
				fastqtobamSingle(CIS,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
		}
		else
		{
			if ( filenames.size() > 2 )
				std::cerr << "[D] warning, ignoring additional input files past first two" << std::endl;
		
			libmaus::aio::CheckedInputStream CIS_1(filenames[0]);
			libmaus::aio::CheckedInputStream CIS_2(filenames[1]);

			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS_1(CIS_1);
				libmaus::lz::BufferedGzipStream BGS_2(CIS_2);
				fastqtobamPair(BGS_1,BGS_2,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
			else
			{
				fastqtobamPair(CIS_1,CIS_2,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
		}

		bamwr.reset();
	}
	else
	{
		::libmaus::bambam::BamParallelWriter::unique_ptr_type bamwr(
			new ::libmaus::bambam::BamParallelWriter(
				std::cout,threads,bamheader,level,Pcbs
			)
		);

		if ( filenames.size() == 0 )
		{
			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS(std::cin);
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
			else
			{
				fastqtobamSingle(std::cin,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
		}
		else if ( filenames.size() == 1 )
		{
			libmaus::aio::CheckedInputStream CIS(filenames[0]);		
			
			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS(CIS);
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);
			}
			else
			{
				fastqtobamSingle(CIS,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);		
			}
		}
		else
		{
			if ( filenames.size() > 2 )
				std::cerr << "[D] warning, ignoring additional input files past first two" << std::endl;
		
			libmaus::aio::CheckedInputStream CIS_1(filenames[0]);
			libmaus::aio::CheckedInputStream CIS_2(filenames[1]);

			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS_1(CIS_1);
				libmaus::lz::BufferedGzipStream BGS_2(CIS_2);
				fastqtobamPair(BGS_1,BGS_2,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);			
			}
			else
			{
				fastqtobamPair(CIS_1,CIS_2,*bamwr,verbose,rgid,qualityoffset,checkquality,namescheme);		
			}
		}

		bamwr.reset();
	
	}

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
}

#if defined(FASTQTOBAMDEBUG)
static void testNameInfo()
{
	{
		std::string const name = "read1";
		libmaus::fastx::SpaceTable const ST;
		NameInfo const NI(name,ST,fastq_name_scheme_generic);
		std::cerr << NI.getName() << " " << NI.ispair << " " << NI.isfirst << std::endl;
	}
	{
		std::string const name = "read1/1";
		libmaus::fastx::SpaceTable const ST;
		NameInfo const NI(name,ST,fastq_name_scheme_generic);
		std::cerr << NI.getName() << " " << NI.ispair << " " << NI.isfirst << std::endl;
	}
	{
		std::string const name = "read1/2";
		libmaus::fastx::SpaceTable const ST;
		NameInfo const NI(name,ST,fastq_name_scheme_generic);
		std::cerr << NI.getName() << " " << NI.ispair << " " << NI.isfirst << std::endl;
	}
	{
		std::string const name = "HWI-ST1289:72:C18CTACXX:8:1101:1116:2072 1:N:0:";
		libmaus::fastx::SpaceTable const ST;
		NameInfo const NI(name,ST,fastq_name_scheme_casava18_single);
		std::cerr << NI.getName() << " " << NI.ispair << " " << NI.isfirst << std::endl;
	}
	{
		std::string const name = "HWI-ST1289:72:C18CTACXX:8:1101:1116:2072 1:N:0:";
		libmaus::fastx::SpaceTable const ST;
		NameInfo const NI(name,ST,fastq_name_scheme_casava18_paired_end);
		std::cerr << NI.getName() << " " << NI.ispair << " " << NI.isfirst << std::endl;
	}
	{
		std::string const name = "HWI-ST1289:72:C18CTACXX:8:1101:1116:2072 2:N:0:";
		libmaus::fastx::SpaceTable const ST;
		NameInfo const NI(name,ST,fastq_name_scheme_casava18_paired_end);
		std::cerr << NI.getName() << " " << NI.ispair << " " << NI.isfirst << std::endl;
	}
}
#endif

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

				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", "compression setting for output BAM file (default: -1, zlib default settings)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "gz=<[0]>", "input is gzip compressed FastQ (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "verbose=<[0]>", "print progress report (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5=<["+::biobambam::Licensing::formatNumber(getDefaultMD5())+"]>", "create md5 check sum (default: 0)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "md5filename=<filename>", "file name for md5 check sum" ) );
				V.push_back ( std::pair<std::string,std::string> ( "I=<[input file name]>", "input file names (standard input if unset)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "threads=<[1]>", "additional BAM encoding helper threads (default: serial encoding)" ) );

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
				V.push_back ( std::pair<std::string,std::string> ( "qualityoffset=<["+::biobambam::Licensing::formatNumber(getDefaultQualityOffset())+"]>", "quality offset (default: 33)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "checkquality=<["+::biobambam::Licensing::formatNumber(getDefaultCheckQuality())+"]>", "check quality (default: true)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("namescheme=<[")+(getDefaultNameScheme())+"]>", "read name scheme (generic, c18s, c18pe)" ) );
				
				::biobambam::Licensing::printMap(std::cerr,V);

				std::cerr << std::endl;
				
				std::cerr << "The I key can be given twice for a pair of synced FastQ files." << std::endl;
				std::cerr << "Any none key=value arguments will be considered as input file names." << std::endl;

				return EXIT_SUCCESS;
			}
		
			
		fastqtobam(arginfo);
		
		std::cerr << "[V] " << libmaus::util::MemUsage() << " wall clock time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
