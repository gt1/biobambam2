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
#include <libmaus/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/BamHeaderUpdate.hpp>
#include <libmaus/bambam/BamFlagBase.hpp>
#include <libmaus/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/util/ToUpperTable.hpp>

#include <libmaus/lz/BufferedGzipStream.hpp>

#include <biobambam/BamBamConfig.hpp>
#include <biobambam/Licensing.hpp>

#include <iomanip>

#include <config.h>

int getDefaultQualityHist()
{
	return 0;
}
int getDefaultQualityOffset()
{
	return 33;
}
int getDefaultQualityMaximum()
{
	return 41;
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
	return libmaus::bambam::BamBlockWriterBaseFactory::checkCompressionLevel(arginfo.getValue<int>("level",getDefaultLevel()));
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

#include <cctype>

struct FastQSeqMapTable
{
	typedef FastQSeqMapTable this_type;
	typedef ::libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef ::libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	::libmaus::autoarray::AutoArray<uint8_t> touppertable;

	FastQSeqMapTable() : touppertable(256)
	{
		for ( unsigned int i = 0; i < 256; ++i )
			if ( isalpha(i) )
				touppertable[i] = ::toupper(i);
			else
				touppertable[i] = i;
		
		// map . to N
		touppertable['.'] = 'N';
	}
	
	char operator()(char const a) const
	{
		return static_cast<char>(touppertable[static_cast<uint8_t>(a)]);
	}

	uint8_t operator()(uint8_t const a) const
	{
		return touppertable[a];
	}
	
	void toupper(std::string & s) const
	{
		for ( std::string::size_type i = 0; i < s.size(); ++i )
			s[i] = static_cast<char>(touppertable[static_cast<uint8_t>(s[i])]);
	}

	std::string operator()(std::string const & s) const
	{
		std::string t = s;
		toupper(t);
		return t;
	}
};

enum fastq_name_scheme_type { 
	fastq_name_scheme_generic,
	fastq_name_scheme_casava18_single,
	fastq_name_scheme_casava18_paired_end,
	fastq_name_scheme_pairedfiles
};

std::ostream & operator<<(std::ostream & out, fastq_name_scheme_type const namescheme)
{
	switch ( namescheme )
	{
		case fastq_name_scheme_generic:
			out << "generic";
			break;
		case fastq_name_scheme_casava18_single:
			out << "c18s";
			break;
		case fastq_name_scheme_casava18_paired_end:
			out << "c18pe";
			break;
		case fastq_name_scheme_pairedfiles:
			out << "pairedfiles";
			break;
		default:
			out << "unknown";
			break;
	}
	return out;
}

static fastq_name_scheme_type parseNameScheme(std::string const & schemename)
{
	if ( schemename == "generic" )
		return fastq_name_scheme_generic;
	else if ( schemename == "c18s" )
		return fastq_name_scheme_casava18_single;
	else if ( schemename == "c18pe" )
		return fastq_name_scheme_casava18_paired_end;
	else if ( schemename == "pairedfiles" )
		return fastq_name_scheme_pairedfiles;
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
					se.getStream() << "malformed read name " << name << " (wrong number of colon separated fields) for name scheme " << namescheme << std::endl;
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
						se.getStream() << "malformed read name " << name << " (malformed fragment id) for name scheme" << namescheme << std::endl;
						se.finish();
						throw se;	
					}
					
					isfirst = (fragid==1);
				}
				
				break;
			}
			case fastq_name_scheme_pairedfiles:
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "pairedfiles name scheme is not supported in NameInfo" << std::endl;
				se.finish();
				throw se;	
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
void checkFastQElement(
	element_type const & element, 
	int const qualityoffset, 
	int const maxvalid = 41
)
{
	for ( uint64_t i = 0; i < element.quality.size(); ++i )
		if ( 
			static_cast<int>(static_cast<uint8_t>(element.quality[i])) - qualityoffset < 0
			||
			static_cast<int>(static_cast<uint8_t>(element.quality[i])) - qualityoffset > maxvalid
		)
	{
		libmaus::exception::LibMausException lme;
		lme.getStream() 
			<< "quality string " << element.quality << " of read " << element.sid 
			<< " is invalid, first invalid symbol is " << element.quality[i]
			<< " deduced value is " << static_cast<int>(static_cast<uint8_t>(element.quality[i])) - qualityoffset;

		if ( static_cast<int>(static_cast<uint8_t>(element.quality[i])) - qualityoffset < 0 )
		{
			lme.getStream() << " (negative values are invalid)\n";
		}
		else
		{
			lme.getStream() << " which is " << ( static_cast<int>(static_cast<uint8_t>(element.quality[i])) - qualityoffset - maxvalid )
				<< " above the maximum of " << maxvalid << "\n";
			lme.getStream() << "If this value is valid, then please set qualitymax to " << static_cast<int>(static_cast<uint8_t>(element.quality[i])) - qualityoffset << " (or above)" << std::endl;
		}
		
		lme.finish();
		throw lme;
	}	
}

template<typename element_type>
void updateQualityHistogram(element_type const & element, int const qualityoffset, uint64_t * const H)
{
	for ( uint64_t i = 0; i < element.quality.size(); ++i )
		H [ static_cast<uint8_t>(element.quality[i]-qualityoffset) ] ++;
}

template<typename writer_type, fastq_name_scheme_type namescheme>
void fastqtobamSingleTemplate(
	std::istream & in, writer_type & bamwr, int const verbose, std::string const & rgid, 
	int const qualityoffset,
	int const maxvalid,
	bool const checkquality,
	uint64_t * const H
	)
{
	typedef ::libmaus::fastx::StreamFastQReaderWrapper reader_type;
	typedef reader_type::pattern_type pattern_type;
	::libmaus::fastx::StreamFastQReaderWrapper fqin(in);	
	libmaus::fastx::SpaceTable ST;
	uint64_t proc = 0;
	FastQSeqMapTable const toup;
	
	switch ( namescheme )
	{
		case fastq_name_scheme_generic:
		case fastq_name_scheme_casava18_single:
		case fastq_name_scheme_casava18_paired_end:
		{
			pattern_type element;
			
			while ( fqin.getNextPatternUnlocked(element) )
			{
				toup.toupper(element.spattern);
				std::string const & name = element.sid;
				NameInfo const NI(name,ST,namescheme);
				
				if ( checkquality )
					checkFastQElement(element,qualityoffset,maxvalid);
				if ( H )
					updateQualityHistogram(element,qualityoffset,H);

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
			break;
		}
		case fastq_name_scheme_pairedfiles:
		{
			pattern_type element_1;
			pattern_type element_2;

			while ( fqin.getNextPatternUnlocked(element_1) )
			{
				bool const ok_2 = fqin.getNextPatternUnlocked(element_2);
				if ( ! ok_2 )
				{
					libmaus::exception::LibMausException ex;
					ex.getStream() << "number of fragments is not even" << std::endl;
					ex.finish();
					throw ex;
				}

				toup.toupper(element_1.spattern);
				toup.toupper(element_2.spattern);

				if ( checkquality )
				{
					checkFastQElement(element_1,qualityoffset,maxvalid);
					checkFastQElement(element_2,qualityoffset,maxvalid);
				}
				if ( H )
				{
					updateQualityHistogram(element_1,qualityoffset,H);
					updateQualityHistogram(element_2,qualityoffset,H);
				}

				// clip read names after first white space
				uint64_t l1 = element_1.sid.size();
				for ( uint64_t i = l1; i; )
					if ( ST.spacetable[static_cast<uint8_t>(element_1.sid[--i])] )
						l1 = i;
				if ( l1 != element_1.sid.size() )
					element_1.sid = element_1.sid.substr(0,l1);

				uint64_t l2 = element_2.sid.size();
				for ( uint64_t i = l2; i; )
					if ( ST.spacetable[static_cast<uint8_t>(element_2.sid[--i])] )
						l2 = i;
				if ( l2 != element_2.sid.size() )
					element_2.sid = element_2.sid.substr(0,l2);
				
				// if read names end in /1 and /2 then chop those off
				if ( 
					element_1.sid.size() >= 2 && element_2.sid.size() >= 2 &&
					element_1.sid[element_1.sid.size()-2] == '/' &&
					element_1.sid[element_1.sid.size()-1] == '1' &&
					element_2.sid[element_2.sid.size()-2] == '/' &&
					element_2.sid[element_2.sid.size()-1] == '2'
				)
				{
					element_1.sid = element_1.sid.substr(0,element_1.sid.size()-2);
					element_2.sid = element_2.sid.substr(0,element_2.sid.size()-2);
				}
				
				if ( element_1.sid != element_2.sid )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "Pair names " << element_1.sid << " and " << element_2.sid << " are not in sync." << std::endl;
					se.finish();
					throw se;
				}
			
				bamwr.encodeAlignment(
					element_1.sid,
					-1,
					-1,
					0,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
					(true ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
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
					element_2.sid,
					-1,
					-1,
					0,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
					(false ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
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
			}
			
			break;
		}
	}
	if ( verbose )
		std::cerr << "[V] " << proc << std::endl;
}

template<typename writer_type>
void fastqtobamSingle(
	std::istream & in, writer_type & bamwr, int const verbose, std::string const & rgid, 
	int const qualityoffset,
	int const maxvalid,
	bool const checkquality,
	uint64_t * const H,
	fastq_name_scheme_type const namescheme
	)
{
	switch ( namescheme )
	{
		case fastq_name_scheme_generic:
			fastqtobamSingleTemplate<writer_type,fastq_name_scheme_generic>(in,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
		case fastq_name_scheme_casava18_single:
			fastqtobamSingleTemplate<writer_type,fastq_name_scheme_casava18_single>(in,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
		case fastq_name_scheme_casava18_paired_end:
			fastqtobamSingleTemplate<writer_type,fastq_name_scheme_casava18_paired_end>(in,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
		case fastq_name_scheme_pairedfiles:
			fastqtobamSingleTemplate<writer_type,fastq_name_scheme_pairedfiles>(in,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
	}
}

template<typename writer_type, fastq_name_scheme_type namescheme>
void fastqtobamPairTemplate(
	std::istream & in_1, std::istream & in_2, 
	writer_type & bamwr, int const verbose, std::string const & rgid, 
	int const qualityoffset, int const maxvalid, bool const checkquality,
	uint64_t * const H
)
{
	typedef ::libmaus::fastx::StreamFastQReaderWrapper reader_type;
	typedef reader_type::pattern_type pattern_type;
	::libmaus::fastx::StreamFastQReaderWrapper fqin_1(in_1);	
	::libmaus::fastx::StreamFastQReaderWrapper fqin_2(in_2);
	pattern_type element_1, element_2;
	libmaus::fastx::SpaceTable ST;
	uint64_t proc = 0;
	FastQSeqMapTable const toup;
	
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
			checkFastQElement(element_1,qualityoffset,maxvalid);
			checkFastQElement(element_2,qualityoffset,maxvalid);
		}
		if ( H )
		{
			updateQualityHistogram(element_1,qualityoffset,H);
			updateQualityHistogram(element_2,qualityoffset,H);
		}

		toup.toupper(element_1.spattern);
		toup.toupper(element_2.spattern);

		std::string const & name_1 = element_1.sid;
		std::string const & name_2 = element_2.sid;
		
		switch ( namescheme )
		{
			case fastq_name_scheme_generic:
			case fastq_name_scheme_casava18_single:
			case fastq_name_scheme_casava18_paired_end:
			{
				NameInfo const NI_1(name_1,ST,namescheme);
				NameInfo const NI_2(name_2,ST,namescheme);
				
				if ( (!NI_1.ispair) || (!NI_1.isfirst) )
				{
					libmaus::exception::LibMausException ex;
					ex.getStream() << "name " << name_1 << " does not look like a first mate read name in name scheme " << namescheme << ". Please check the namescheme parameter." << std::endl;
					ex.finish();
					throw ex;		
				}
				if ( NI_2.isfirst )
				{
					libmaus::exception::LibMausException ex;
					ex.getStream() << "name " << name_2 << " does not look like a second mate read name in name scheme " << namescheme << ". Please check the namescheme parameter." << std::endl;
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

				break;
			}
			case fastq_name_scheme_pairedfiles:
			{
				// clip read names after first white space
				uint64_t l1 = element_1.sid.size();
				for ( uint64_t i = l1; i; )
					if ( ST.spacetable[static_cast<uint8_t>(element_1.sid[--i])] )
						l1 = i;
				if ( l1 != element_1.sid.size() )
					element_1.sid = element_1.sid.substr(0,l1);

				uint64_t l2 = element_2.sid.size();
				for ( uint64_t i = l2; i; )
					if ( ST.spacetable[static_cast<uint8_t>(element_2.sid[--i])] )
						l2 = i;
				if ( l2 != element_2.sid.size() )
					element_2.sid = element_2.sid.substr(0,l2);
				
				// if read names end in /1 and /2 then chop those off
				if ( 
					element_1.sid.size() >= 2 && element_2.sid.size() >= 2 &&
					element_1.sid[element_1.sid.size()-2] == '/' &&
					element_1.sid[element_1.sid.size()-1] == '1' &&
					element_2.sid[element_2.sid.size()-2] == '/' &&
					element_2.sid[element_2.sid.size()-1] == '2'
				)
				{
					element_1.sid = element_1.sid.substr(0,element_1.sid.size()-2);
					element_2.sid = element_2.sid.substr(0,element_2.sid.size()-2);
				}
				
				if ( element_1.sid != element_2.sid )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "Pair names " << element_1.sid << " and " << element_2.sid << " are not in sync." << std::endl;
					se.finish();
					throw se;
				}
			
				bamwr.encodeAlignment(
					element_1.sid,
					-1,
					-1,
					0,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
					(true ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
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
					element_2.sid,
					-1,
					-1,
					0,
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FPAIRED |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FUNMAP |
					libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FMUNMAP |
					(false ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD1 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREAD2),
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
			}
		}
		
		proc += 2;
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
	int const qualityoffset, int const maxvalid, bool const checkquality,
	uint64_t * const H,
	fastq_name_scheme_type const namescheme
)
{
	switch ( namescheme )
	{
		case fastq_name_scheme_generic:
			fastqtobamPairTemplate<writer_type,fastq_name_scheme_generic>(in_1,in_2,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
		case fastq_name_scheme_casava18_single:
			fastqtobamPairTemplate<writer_type,fastq_name_scheme_casava18_single>(in_1,in_2,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
		case fastq_name_scheme_casava18_paired_end:
			fastqtobamPairTemplate<writer_type,fastq_name_scheme_casava18_paired_end>(in_1,in_2,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
		case fastq_name_scheme_pairedfiles:
			fastqtobamPairTemplate<writer_type,fastq_name_scheme_pairedfiles>(in_1,in_2,bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H);
			break;
	}
}

void fastqtobam(libmaus::util::ArgInfo const & arginfo)
{
	std::vector<std::string> filenames = arginfo.getPairValues("I");
	for ( uint64_t i = 0; i < arginfo.restargs.size(); ++i )
		filenames.push_back(arginfo.restargs[i]);

	int const qualityhistflag = arginfo.getValue<int>("qualityhist",getDefaultQualityHist());

	uint64_t qualityhist[256];
	std::fill(&qualityhist[0],&qualityhist[sizeof(qualityhist)/sizeof(qualityhist[0])],0ull);
	uint64_t * const H = qualityhistflag ? (&qualityhist[0]) : 0;

	int const level = getLevel(arginfo);
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());
	bool const gz = arginfo.getValue<int>("gz",getDefaultGz());
	unsigned int const threads = arginfo.getValue<unsigned int>("threads",1);
	int const qualityoffset = arginfo.getValue<int>("qualityoffset",getDefaultQualityOffset());
	int const maxvalid = arginfo.getValue<int>("qualitymax",getDefaultQualityMaximum());
	
	uint64_t const bamqualmax = ('~'-'!'+1);
	if ( maxvalid > static_cast<int>(bamqualmax) )
	{
		libmaus::exception::LibMausException lme;
		lme.getStream() << "Quality values of more than " << bamqualmax << " cannot be represented in BAM files\n";
		lme.finish();
		throw lme;
	}
	
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
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
			else
			{
				fastqtobamSingle(std::cin,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
		}
		else if ( filenames.size() == 1 )
		{
			libmaus::aio::CheckedInputStream CIS(filenames[0]);		
			
			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS(CIS);
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
			else
			{
				fastqtobamSingle(CIS,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
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
				fastqtobamPair(BGS_1,BGS_2,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
			else
			{
				fastqtobamPair(CIS_1,CIS_2,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
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
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
			else
			{
				fastqtobamSingle(std::cin,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
		}
		else if ( filenames.size() == 1 )
		{
			libmaus::aio::CheckedInputStream CIS(filenames[0]);		
			
			if ( gz )
			{
				libmaus::lz::BufferedGzipStream BGS(CIS);
				fastqtobamSingle(BGS,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);
			}
			else
			{
				fastqtobamSingle(CIS,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);		
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
				fastqtobamPair(BGS_1,BGS_2,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);			
			}
			else
			{
				fastqtobamPair(CIS_1,CIS_2,*bamwr,verbose,rgid,qualityoffset,maxvalid,checkquality,H,namescheme);		
			}
		}

		bamwr.reset();
	
	}

	if ( Pmd5cb )
	{
		Pmd5cb->saveDigestAsFile(md5filename);
	}
	
	if ( H )
	{
		uint64_t const allbases = std::accumulate(H,H+256,0ull);
		double const dallbases = allbases;
		double acc = 0;
		
		for ( uint64_t ii = 0; ii < 256; ++ii )
		{
			uint64_t i = 256-ii-1;
			
			if ( H[i] )
			{
				acc += H[i]/dallbases;
				std::cerr << "[H]\t" << static_cast<char>(i+33) << "\t" << H[i] << "\t" << H[i]/dallbases << "\t" << acc << std::endl;
			}
		}
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

				V.push_back ( std::pair<std::string,std::string> ( "level=<[-1]>", libmaus::bambam::BamBlockWriterBaseFactory::getBamOutputLevelHelpText() ) );
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
				V.push_back ( std::pair<std::string,std::string> ( "qualitymax=<["+::biobambam::Licensing::formatNumber(getDefaultQualityMaximum())+"]>", "maximum valid quality value (default: 41)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "checkquality=<["+::biobambam::Licensing::formatNumber(getDefaultCheckQuality())+"]>", "check quality (default: true)" ) );
				V.push_back ( std::pair<std::string,std::string> ( std::string("namescheme=<[")+(getDefaultNameScheme())+"]>", "read name scheme (generic, c18s, c18pe, pairedfiles)" ) );
				V.push_back ( std::pair<std::string,std::string> ( "qualityhist=<["+::biobambam::Licensing::formatNumber(getDefaultQualityHist())+"]>", "compute quality histogram and print it on standard error" ) );
				
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
