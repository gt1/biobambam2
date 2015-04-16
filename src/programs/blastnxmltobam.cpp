/**
    biobambam
    Copyright (C) 2004-2014 German Tischler
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
#include <xercesc/framework/Wrapper4InputSource.hpp>
#include <xercesc/sax/AttributeList.hpp>
#include <xercesc/sax/DocumentHandler.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/InputSource.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/util/BinInputStream.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/util/XMLUniDefs.hpp>

#include <libmaus/fastx/acgtnMap.hpp>
#include <libmaus/bambam/BamWriter.hpp>
#include <libmaus/bambam/CramRange.hpp>
#include <libmaus/util/ArgInfo.hpp>

#include <iostream>
#include <map>

#include <config.h>

struct XercesUtf8Transcoder
{
	typedef XercesUtf8Transcoder this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	xercesc::XMLTransService * ts;
	xercesc::XMLTranscoder * utf8transcoder;
	
	XercesUtf8Transcoder()
	: ts(0), utf8transcoder(0)
	{
		ts = (xercesc::XMLPlatformUtils::fgTransService);

		xercesc::XMLTransService::Codes transretcode;

		char const * ucharset = "UTF-8";
		utf8transcoder = ts->makeNewTranscoderFor("UTF-8",transretcode,8192u,xercesc::XMLPlatformUtils::fgMemoryManager);
		if ( transretcode != xercesc::XMLTransService::Ok )
		{
			libmaus::exception::LibMausException lme;
			lme.getStream() << "Failed to instantiate xerces transcoder for code " << ucharset << std::endl;
			lme.finish();
			throw lme;
		}
	}
	~XercesUtf8Transcoder()
	{
		delete utf8transcoder;
	}

	std::string transcodeStringToUtf8(XMLCh const * text)
	{
		unsigned int len = xercesc::XMLString::stringLen(text);
		std::ostringstream out;
		XMLByte buffer[1024];
		while ( len )
		{
			#if XERCES_VERSION_MAJOR >= 3
			XMLSize_t eaten = 0;
			#else
			unsigned int eaten = 0;
			#endif

			unsigned int produced = 
				utf8transcoder->transcodeTo( 
				text,len,&buffer[0], 
				sizeof(buffer),eaten,
				xercesc::XMLTranscoder::UnRep_RepChar);
				len -= eaten; 
				text += eaten;
			out.write ( reinterpret_cast<char const *>(buffer), produced );
		}
		return out.str();
	}
};

//! A class that wraps std::istream for Xerces.
class StdInputBinStream : public xercesc::BinInputStream
{
	private:
	//! Input stream.
	std::istream & in;
	//! Position in stream.
	unsigned int pos;

	public:
	//! Constructor.
	/**
	 * Construct with an input stream.
	 *
	 * @param rin Input stream.
	 */
	StdInputBinStream(std::istream & rin) : in(rin), pos(0) {}
	//! Destructor.
	virtual ~StdInputBinStream() {}

	//! Get current position
	/**
	 * @return Current posistion.
	 */
	#if XERCES_VERSION_MAJOR >= 3
	virtual XMLFilePos curPos() const { return pos; }
	#else
	virtual unsigned int curPos() const { return pos; }	
	#endif
	//! Read bytes from stream.
	/**
	 * Fill the provided buffer with input bytes from stream.
	 *
	 * @param toFill Buffer to fill.
	 * @param maxToRead Maximum number of bytes to be read.
	 * @return Number of bytes that have been read and
	 *         transfered into toFill.
	 */
	#if XERCES_VERSION_MAJOR >= 3
	virtual XMLSize_t readBytes(XMLByte * const toFill, const XMLSize_t maxToRead)
	#else
	virtual unsigned int readBytes(XMLByte * const toFill, const unsigned int maxToRead)
	#endif
	{
		in.read(reinterpret_cast<char *>(toFill),maxToRead);
		pos += in.gcount();
		return in.gcount();
	}
	
	#if XERCES_VERSION_MAJOR >= 3
	virtual const XMLCh* getContentType() const
	{
		static const XMLCh type[] = {
			xercesc::chLatin_t, xercesc::chLatin_e, xercesc::chLatin_x, xercesc::chLatin_t,
			xercesc::chForwardSlash,
			xercesc::chLatin_x, xercesc::chLatin_m, xercesc::chLatin_l, xercesc::chNull };
		return &type[0];
	}
	#endif
};

//! Input stream wrapper for Xerces inputsource.
class StdISOInputSource : public xercesc::InputSource
{
	private:
	//! Input stream.
	std::istream & in;

	public:
	//! Constructor.
	/**
	 * @param rin Input stream.
	 */
	StdISOInputSource(std::istream & rin) : in(rin) {}
	//! Destructor.
	virtual ~StdISOInputSource() {}

	//! Make binary input stream.
	/**
	 * Make a binary input stream. The current object does
	 * not own the stream, it has to be deallocated by
	 * the caller.
	 * 
	 * @return Binary input stream.
	 */
	virtual xercesc::BinInputStream * makeStream() const {
		return new StdInputBinStream(in);
	}
};

#include <libmaus/util/ToUpperTable.hpp>

struct BlastNDocumentHandler : public xercesc::DocumentHandler, public xercesc::ErrorHandler
{
	libmaus::util::ToUpperTable const toup;
	
	std::map<std::string,std::string> const & ref;
	std::map<std::string,std::string> const & queries;

	XercesUtf8Transcoder utf8transcoder;
	
	bool readNameGatheringActive;
	std::string readName;
	bool readNameObtained;

	bool hitDefObtained;
	bool hitDefGatheringActive;
	std::string hitDef;

	bool hitLenObtained;
	bool hitLenGatheringActive;
	std::string hitLen;

	bool hspBitScoreObtained;
	bool hspBitScoreGatheringActive;
	std::string hspBitScore;

	bool hspScoreObtained;
	bool hspScoreGatheringActive;
	std::string hspScore;

	bool hspEvalueObtained;
	bool hspEvalueGatheringActive;
	std::string hspEvalue;

	bool hspQueryFromObtained;
	bool hspQueryFromGatheringActive;
	std::string hspQueryFrom;

	bool hspQueryToObtained;
	bool hspQueryToGatheringActive;
	std::string hspQueryTo;

	bool hspHitFromObtained;
	bool hspHitFromGatheringActive;
	std::string hspHitFrom;

	bool hspHitToObtained;
	bool hspHitToGatheringActive;
	std::string hspHitTo;

	bool hspQueryFrameObtained;
	bool hspQueryFrameGatheringActive;
	std::string hspQueryFrame;

	bool hspHitFrameObtained;
	bool hspHitFrameGatheringActive;
	std::string hspHitFrame;

	bool hspIdentityObtained;
	bool hspIdentityGatheringActive;
	std::string hspIdentity;

	bool hspPositiveObtained;
	bool hspPositiveGatheringActive;
	std::string hspPositive;

	bool hspGapsObtained;
	bool hspGapsGatheringActive;
	std::string hspGaps;

	bool hspAlignLenObtained;
	bool hspAlignLenGatheringActive;
	std::string hspAlignLen;

	bool hspQSeqObtained;
	bool hspQSeqGatheringActive;
	std::string hspQSeq;

	bool hspHSeqObtained;
	bool hspHSeqGatheringActive;
	std::string hspHSeq;
	
	uint64_t hspId;

	std::map<std::string,uint64_t> refnametoid;
	std::map<std::string,uint64_t> queriesnametoid;
	libmaus::bambam::BamWriter & bamwriter;
	
	double hitFirstScore;
	double hitFrac;

	#if 0
	<Hsp_identity>21727</Hsp_identity>
        <Hsp_positive>21727</Hsp_positive>
        <Hsp_gaps>17</Hsp_gaps>
        <Hsp_align-len>21778</Hsp_align-len>
	#endif            
	
	std::vector<libmaus::bambam::CramRange> const * ranges;
	
	bool inRange(std::string const & refname, int64_t const hitstart, int64_t const hitend)
	{
		if ( ! ranges )
			return true;
		
		for ( uint64_t i = 0; i < ranges->size(); ++i )
			if ( 
				!(((*ranges)[i]).intersect(libmaus::bambam::CramRange(refname,hitstart,hitend)).empty())
			)
				return true;
		
		std::cerr << "[E] dropping " << refname << ":" << hitstart << "-" << hitend << std::endl;
		
		return false;
	}

	BlastNDocumentHandler(
		std::map<std::string,std::string> const & rref,
		std::map<std::string,std::string> const & rqueries,
		std::map<std::string,uint64_t> rrefnametoid,
		std::map<std::string,uint64_t> rqueriesnametoid,
		libmaus::bambam::BamWriter & rbamwriter,
		double const rhitFrac,
		std::vector<libmaus::bambam::CramRange> const * rranges
	) : ref(rref), queries(rqueries), utf8transcoder(), readNameGatheringActive(false), readName(), readNameObtained(false), 
		hitDefObtained(false), hitDefGatheringActive(false), hitDef(),
		hitLenObtained(false), hitLenGatheringActive(false), hitLen(),	
		hspBitScoreObtained(false), hspBitScoreGatheringActive(false), hspBitScore(),
		hspScoreObtained(false), hspScoreGatheringActive(false), hspScore(),
		hspEvalueObtained(false), hspEvalueGatheringActive(false), hspEvalue(),
		hspQueryFromObtained(false), hspQueryFromGatheringActive(false), hspQueryFrom(),
		hspQueryToObtained(false), hspQueryToGatheringActive(false), hspQueryTo(),
		hspHitFromObtained(false), hspHitFromGatheringActive(false), hspHitFrom(),
		hspHitToObtained(false), hspHitToGatheringActive(false), hspHitTo(),
		hspQueryFrameObtained(false), hspQueryFrameGatheringActive(false), hspQueryFrame(),
		hspHitFrameObtained(false), hspHitFrameGatheringActive(false), hspHitFrame(),
		hspIdentityObtained(false), hspIdentityGatheringActive(false), hspIdentity(),
		hspPositiveObtained(false), hspPositiveGatheringActive(false), hspPositive(),
		hspGapsObtained(false), hspGapsGatheringActive(false), hspGaps(),
		hspAlignLenObtained(false), hspAlignLenGatheringActive(false), hspAlignLen(),
		hspQSeqObtained(false), hspQSeqGatheringActive(false), hspQSeq(),
		hspHSeqObtained(false), hspHSeqGatheringActive(false), hspHSeq(),
		hspId(0),
		refnametoid(rrefnametoid), queriesnametoid(rqueriesnametoid),
		bamwriter(rbamwriter),
		hitFirstScore(-1),
		hitFrac(rhitFrac),
		ranges(rranges)
	{
	
	}
	virtual ~BlastNDocumentHandler()
	{
	
	}

#if XERCES_VERSION_MAJOR < 3
	virtual void characters
	(
		const   XMLCh* const    chars, 
		const unsigned int    /* length */
	)
#else
	virtual void characters
	(
		const   XMLCh* const    chars, 
		XMLSize_t    /* length */
        )
#endif
	{
		if ( readNameGatheringActive )
			readName += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hitDefGatheringActive )
			hitDef += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hitLenGatheringActive )
			hitLen += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspBitScoreGatheringActive )
			hspBitScore += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspScoreGatheringActive )
			hspScore += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspEvalueGatheringActive )
			hspEvalue += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspQueryFromGatheringActive )
			hspQueryFrom += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspQueryToGatheringActive )
			hspQueryTo += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspHitFromGatheringActive )
			hspHitFrom += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspHitToGatheringActive )
			hspHitTo += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspQueryFrameGatheringActive )
			hspQueryFrame += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspHitFrameGatheringActive )
			hspHitFrame += utf8transcoder.transcodeStringToUtf8(chars);

		if ( hspIdentityGatheringActive )
			hspIdentity += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspPositiveGatheringActive )
			hspPositive += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspGapsGatheringActive )
			hspGaps += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspAlignLenGatheringActive )
			hspAlignLen += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspQSeqGatheringActive )
			hspQSeq += utf8transcoder.transcodeStringToUtf8(chars);
		if ( hspHSeqGatheringActive )
			hspHSeq += utf8transcoder.transcodeStringToUtf8(chars);

		#if 0
		<Hsp_identity>21727</Hsp_identity>
		<Hsp_positiv>21727</Hsp_positive>
                <Hsp_gaps>17</Hsp_gaps>
                <Hsp_align-len>21778</Hsp_align-len>
                #endif
	}

#if XERCES_VERSION_MAJOR < 3
	virtual void ignorableWhitespace    (const   XMLCh* const    /*chars*/, const unsigned int    /*length*/)
#else
	virtual void ignorableWhitespace    (const   XMLCh* const    /*chars*/, XMLSize_t    /*length*/)
#endif
	{
	
	}

	virtual void processingInstruction
	(
		const   XMLCh* const    /* target */
		, const XMLCh* const    /* data */
	)
	{
	
	}
	virtual void resetDocument()
	{
	
	}
	virtual void setDocumentLocator(const xercesc::Locator* const /*locator*/)
	{
	
	}
	virtual void startDocument()
	{
		// std::cerr << "start of document" << std::endl;
	}
	virtual void endDocument ()
	{
		// std::cerr << "end of document" << std::endl;
	}
	virtual void startElement
	(
		const   XMLCh* const    xname,
		xercesc::AttributeList&  /* xattrs */
	)
	{
		std::string const tag = utf8transcoder.transcodeStringToUtf8(xname);
		
		if ( tag == "Iteration" )
		{
			// std::cerr << "new query" << std::endl;
			readNameObtained = false;
			hspId = 0;
			// std::cerr << "***" << std::endl;
		}
		if ( tag == "Iteration_query-def" )
		{
			readNameGatheringActive = true;	
			readName.clear();
		}
		if ( tag == "Hit" )
		{
			hitDefObtained = false;	
			hitDef.clear();
			hitLenObtained = false;
			hitLen.clear();
			hitFirstScore = -1;
		}
		if ( tag == "Hit_num" )
		{
		
		}
		if ( tag == "Hit_def" )
		{
			hitDefGatheringActive = true;
		}
		if ( tag == "Hit_len" )
		{
			hitLenGatheringActive = true;
		}
		if ( tag == "Hsp" )
		{
			hspBitScoreObtained = false;
			hspBitScore.clear();
			hspScoreObtained = false;
			hspScore.clear();
			hspEvalueObtained = false;
			hspEvalue.clear();
			hspQueryFromObtained = false;
			hspQueryFrom.clear();
			hspQueryToObtained = false;
			hspQueryTo.clear();
			hspHitFromObtained = false;
			hspHitFrom.clear();
			hspHitToObtained = false;
			hspHitTo.clear();
			hspQueryFrame.clear();
			hspQueryFrameObtained = false;
			hspHitFrame.clear();
			hspHitFrameObtained = false;

			hspIdentity.clear();
			hspIdentityObtained = false;

			hspPositive.clear();
			hspPositiveObtained = false;

			hspGaps.clear();
			hspGapsObtained = false;

			hspAlignLen.clear();
			hspAlignLenObtained = false;

			hspQSeq.clear();
			hspQSeqObtained = false;

			hspHSeq.clear();
			hspHSeqObtained = false;
		}
		if ( tag == "Hsp_bit-score" )
		{
			hspBitScoreGatheringActive = true;
		}
		if ( tag == "Hsp_score" )
		{
			hspScoreGatheringActive = true;
		}
		if ( tag == "Hsp_evalue" )
		{
			hspEvalueGatheringActive = true;
		}
		if ( tag == "Hsp_query-from" )
		{
			hspQueryFromGatheringActive = true;
		}
		if ( tag == "Hsp_query-to" )
		{
			hspQueryToGatheringActive = true;
		}
		if ( tag == "Hsp_hit-from" )
		{
			hspHitFromGatheringActive = true;
		}
		if ( tag == "Hsp_hit-to" )
		{
			hspHitToGatheringActive = true;
		}
		if ( tag == "Hsp_query-frame" )
		{
			hspQueryFrameGatheringActive = true;
		}
		if ( tag == "Hsp_hit-frame" )
		{
			hspHitFrameGatheringActive = true;
		}
		if ( tag == "Hsp_identity" )
		{
			hspIdentityGatheringActive = true;
		}
		if ( tag == "Hsp_positive" )
		{
			hspPositiveGatheringActive = true;
		}
		if ( tag == "Hsp_gaps" )
		{
			hspGapsGatheringActive = true;
		}
		if ( tag == "Hsp_align-len" )
		{
			hspAlignLenGatheringActive = true;
		}
		if ( tag == "Hsp_qseq" )
		{
			hspQSeqGatheringActive = true;
		}
		if ( tag == "Hsp_hseq" )
		{
			hspHSeqGatheringActive = true;
		}


		#if 0
		<Hsp_identity>21727</Hsp_identity>
		<Hsp_positive>21727</Hsp_positive>
                <Hsp_gaps>17</Hsp_gaps>
                <Hsp_align-len>21778</Hsp_align-len>
                #endif

#if 0
<Hsp>
<Hsp_num>1</Hsp_num>
<Hsp_bit-score>17429.8</Hsp_bit-score>
<Hsp_score>9438</Hsp_score>
<Hsp_evalue>0</Hsp_evalue>
<Hsp_query-from>1</Hsp_query-from>
<Hsp_query-to>9603</Hsp_query-to>
<Hsp_hit-from>161746186</Hsp_hit-from>
<Hsp_hit-to>161736547</Hsp_hit-to>
#endif
                                                                        

		// std::cerr << "start of element " << utf8transcoder.transcodeStringToUtf8(xname) << std::endl;
	}
	
	template<typename number_type>
	static number_type parseNumber(std::string const & s)
	{
		std::istringstream istr(s);
		number_type i;
		istr >>i;
		
		if ( ! istr )
		{
			libmaus::exception::LibMausException lme;
			lme.getStream() << "Failed to parse " << s << " as number." << std::endl;
			lme.finish();
			throw lme;
		}
		
		return i;
	}
	
	virtual void endElement(const XMLCh* const xname)
	{
		std::string const tag = utf8transcoder.transcodeStringToUtf8(xname);
	
		if ( tag == "Iteration_query-def" )
		{
			readNameGatheringActive = false;
			readNameObtained = true;
			// std::cerr << "read name " << readName << std::endl;
		}
		if ( tag == "Hit_def" )
		{
			hitDefGatheringActive = false;
			hitDefObtained = true;
			// std::cerr << "hit def " << hitDef << std::endl;
		}
		if ( tag == "Hit_len" )
		{
			hitLenGatheringActive = false;
			hitLenObtained = true;
			// std::cerr << "hit len " << hitLen << std::endl;
		}
		if ( tag == "Hsp" )
		{
			bool ok = 
				hspBitScoreObtained &&
				hspScoreObtained &&
				hspEvalueObtained &&
				hspQueryFromObtained &&
				hspQueryToObtained &&
				hspHitFromObtained &&
				hspHitToObtained &&
				hspQueryFrameObtained &&
				hspHitFrameObtained &&
				hspIdentityObtained &&
				hspPositiveObtained &&
				hspGapsObtained &&
				hspAlignLenObtained &&
				hspQSeqObtained &&
				hspHSeqObtained &&
				(queries.find(readName) != queries.end());

			int64_t const thisHitScore = hspScoreObtained ?  parseNumber<int64_t>(hspScore) : -1;
				
			if ( hspId == 0 )
				hitFirstScore = thisHitScore;

			// reference
			std::map<std::string,std::string>::const_iterator hita = 
				ok ? ref.find(hitDef) : std::map<std::string,std::string>::const_iterator();
			// hit coord
			int64_t hitFrom = ok ? parseNumber<int64_t>(hspHitFrom) : -1;
			int64_t hitTo = ok ? parseNumber<int64_t>(hspHitTo) : -1;

			// hit start and end
			int64_t hitStart = std::min(hitFrom,hitTo)-1;
			int64_t hitEnd = std::max(hitFrom,hitTo)-1;
			int64_t hitLen = hitEnd-hitStart+1;
			
			if ( ok && 
				(hspId == 0 || (thisHitScore >= hitFrac * hitFirstScore)) && 
				inRange(hita->first, hitStart, hitEnd)
			)
			{
				// hita->first refseq

				int64_t queryFrame = parseNumber<int64_t>(hspQueryFrame);
				int64_t hitFrame = parseNumber<int64_t>(hspHitFrame);
				bool const rc = (queryFrame * hitFrame) < 0;
				// query coord
				int64_t queryFrom = parseNumber<int64_t>(hspQueryFrom);
				int64_t queryTo = parseNumber<int64_t>(hspQueryTo);
				
				
				// query start and end
				int64_t queryStart = std::min(queryFrom,queryTo)-1;
				int64_t queryEnd = std::max(queryFrom,queryTo)-1;
				int64_t queryLen = queryEnd-queryStart+1;
				int64_t queryFrontClip = queryStart;
				std::map<std::string,std::string>::const_iterator qita = queries.find(readName);
				assert ( qita != queries.end() );
				int64_t queryBackClip = qita->second.size() - (queryFrontClip + queryLen);
				
				std::cerr 
					<< readName << "[" << hspId << "]" << " queryFrame " << queryFrame << " hitFrame " << hitFrame 
					<< " query coord [" << hspQueryFrom << "," << hspQueryTo << "]"
					<< " hit coord [" << hspHitFrom << "," << hspHitTo << "]"
					<< std::endl;
					
				uint64_t qlen = 0;
				for ( uint64_t i = 0; i < hspQSeq.size(); ++i )
					qlen += hspQSeq[i] != '-';

				uint64_t hlen = 0;
				for ( uint64_t i = 0; i < hspHSeq.size(); ++i )
					hlen += hspHSeq[i] != '-';
						
				
				if ( qita != queries.end() && hita != ref.end() )
				{
					std::string const & H = hita->second;
					std::string const & Q = qita->second;
					
					std::string hsub = H.substr(hitStart,hitLen);
					std::string qsub = Q.substr(queryStart,queryLen);
					
					if ( hitFrame < 0 )
					{
						std::reverse(hsub.begin(),hsub.end());
						for ( uint64_t i = 0; i < hsub.size(); ++i )
							switch ( hsub[i] )
							{
								case 'a': case 'A': hsub[i] = 'T'; break;
								case 'c': case 'C': hsub[i] = 'G'; break;
								case 'g': case 'G': hsub[i] = 'C'; break;
								case 't': case 'T': hsub[i] = 'A'; break;
								default: hsub[i] = 'N'; break;
								break;
							}
					}
					if ( queryFrame < 0 )
					{
						std::reverse(qsub.begin(),qsub.end());
						for ( uint64_t i = 0; i < qsub.size(); ++i )
							switch ( qsub[i] )
							{
								case 'a': case 'A': qsub[i] = 'T'; break;
								case 'c': case 'C': qsub[i] = 'G'; break;
								case 'g': case 'G': qsub[i] = 'C'; break;
								case 't': case 'T': qsub[i] = 'A'; break;
								default: qsub[i] = 'N'; break;
								break;
							}
					}
					
					std::vector<char> ops;
					
					for ( int64_t i = 0; i < queryFrontClip; ++i )
						ops.push_back('S');

					uint64_t ih = 0, iq = 0;
					for ( uint64_t i = 0; i < hspQSeq.size(); ++i )
					{
						assert ( hspQSeq[i] != '-' || hspHSeq[i] != '-' );
						
						// symbol not in query sequence -> deleted from reference
						if ( hspQSeq[i] == '-' )
						{
							ops.push_back('D');
						}
						// symbol in query but not in reference -> insertion into reference
						else if ( hspHSeq[i] == '-' )
						{
							ops.push_back('I');
						}
						else if ( hspQSeq[i] == hspHSeq[i] )
						{
							// match
							ops.push_back('=');
						}
						else
						{
							ops.push_back('X');
						}
					
						if ( hspQSeq[i] != '-' )
						{
							bool const ok = ( 
								toup(static_cast<uint8_t>(hspQSeq[i]))
								== 
								toup(static_cast<uint8_t>(qsub[iq++]))
							);
							
							if ( ! ok )
								std::cerr << "[W] " << toup(static_cast<uint8_t>(hspQSeq[i])) << " != "
									<< toup(static_cast<uint8_t>(qsub[iq-1])) << std::endl;
						}
						if ( hspHSeq[i] != '-' )
							assert ( 
								toup(static_cast<uint8_t>(hspHSeq[i]))
								== 
								toup(static_cast<uint8_t>(hsub[ih++]))
							);
														
						#if 0
						std::cerr << "(" << hspQSeq[i] << "," << hspHSeq[i] << ",";
						
						if ( hspQSeq[i] == '-' )
							std::cerr << "-";
						else
							std::cerr << qsub[iq++];
						
						std::cerr << ",";

						if ( hspHSeq[i] == '-' )
							std::cerr << "-";
						else
							std::cerr << hsub[ih++];

						std::cerr << ")";
						#endif
					}

					for ( int64_t i = 0; i < queryBackClip; ++i )
						ops.push_back('S');
						
					if ( rc )
						std::reverse(ops.begin(),ops.end());

					std::string bamquery = qita->second;
					if ( rc )
						bamquery = libmaus::fastx::reverseComplementUnmapped(bamquery);
					
					std::vector < std::pair<char,uint64_t> > opruns;
					
					uint64_t low = 0;
					while ( low != ops.size() )
					{
						uint64_t high = low;
						while ( high != ops.size() && ops[high] == ops[low] )
							++high;
							
						opruns.push_back(std::pair<char,uint64_t>(ops[low],high-low));
							
						low = high;
					}
					
					std::ostringstream cigarostr;
					for ( uint64_t i = 0; i < opruns.size(); ++i )
					{
						std::cerr << "(" << opruns[i].first << "," << opruns[i].second << ")";
						cigarostr << opruns[i].second << opruns[i].first;
					}
					std::cerr << std::endl;

					bamwriter.encodeAlignment(
						readName,
						refnametoid.find(hita->first)->second,
						hitStart,
						0, // mapq
						(rc ? libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FREVERSE : 0)
						|
						(hspId == 0 ? 0 : libmaus::bambam::BamFlagBase::LIBMAUS_BAMBAM_FSECONDARY)
						, // flags
						cigarostr.str(),
						-1,
						-1,
						0,
						bamquery,
						std::string(qita->second.size(),255),
						0
					);
					bamwriter.putAuxNumber("AS", 'i', thisHitScore);
					bamwriter.putAuxNumber("ZA", 'f', parseNumber<double>(hspBitScore));
					bamwriter.putAuxNumber("ZB", 'i', parseNumber<int64_t>(hspScore));
					bamwriter.putAuxNumber("ZC", 'f', parseNumber<double>(hspEvalue));
					bamwriter.putAuxNumber("ZD", 'i', parseNumber<int64_t>(hspQueryFrom));
					bamwriter.putAuxNumber("ZE", 'i', parseNumber<int64_t>(hspQueryTo));
					bamwriter.putAuxNumber("ZF", 'i', parseNumber<int64_t>(hspHitFrom));
					bamwriter.putAuxNumber("ZG", 'i', parseNumber<int64_t>(hspHitTo));
					bamwriter.putAuxNumber("ZH", 'i', parseNumber<int64_t>(hspQueryFrame));
					bamwriter.putAuxNumber("ZI", 'i', parseNumber<int64_t>(hspHitFrame));
					bamwriter.putAuxNumber("ZJ", 'i', parseNumber<int64_t>(hspIdentity));
					bamwriter.putAuxNumber("ZK", 'i', parseNumber<int64_t>(hspPositive));
					bamwriter.putAuxNumber("ZL", 'i', parseNumber<int64_t>(hspGaps));
					bamwriter.putAuxNumber("ZM", 'i', parseNumber<int64_t>(hspAlignLen));
					bamwriter.commit();
					
					#if 0
					std::cerr << std::endl;
					#endif
				}
				
				#if 0
				if ( hitFrame < 0 )
				{
					hitFrame = -hitFrame;
					queryFrame = -queryFrame;
					
					std::reverse(hspQSeq.begin(),hspQSeq.end());
					std::reverse(hspHSeq.begin(),hspHSeq.end());
				}
				#endif
			}
			
			hspId += 1;
		}
		if ( tag == "Hsp_bit-score" )
		{
			hspBitScoreGatheringActive = false;
			hspBitScoreObtained = true;
			// std::cerr << "hspBitScore " << hspBitScore << std::endl;
		}
		if ( tag == "Hsp_score" )
		{
			hspScoreGatheringActive = false;
			hspScoreObtained = true;
			// std::cerr << "hspScore " << hspScore << std::endl;
		}
		if ( tag == "Hsp_evalue" )
		{
			hspEvalueGatheringActive = false;
			hspEvalueObtained = true;
			// std::cerr << "hspEvalue " << hspEvalue << std::endl;
		}
		if ( tag == "Hsp_query-from" )
		{
			hspQueryFromGatheringActive = false;
			hspQueryFromObtained = true;
			// std::cerr << "hspQueryFrom " << hspQueryFrom << std::endl;
		}
		if ( tag == "Hsp_query-to" )
		{
			hspQueryToGatheringActive = false;
			hspQueryToObtained = true;
			// std::cerr << "hspQueryTo " << hspQueryTo << std::endl;
		}
		if ( tag == "Hsp_hit-from" )
		{
			hspHitFromGatheringActive = false;
			hspHitFromObtained = true;
			// std::cerr << "hspHitFrom " << hspHitFrom << std::endl;
		}
		if ( tag == "Hsp_hit-to" )
		{
			hspHitToGatheringActive = false;
			hspHitToObtained = true;
			// std::cerr << "hspHitTo " << hspHitTo << std::endl;
		}
		if ( tag == "Hsp_query-frame" )
		{
			hspQueryFrameGatheringActive = false;
			hspQueryFrameObtained = true;
			// std::cerr << "hspQueryFrame " << hspQueryFrame << std::endl;
		}
		if ( tag == "Hsp_hit-frame" )
		{
			hspHitFrameGatheringActive = false;
			hspHitFrameObtained = true;
			// std::cerr << "hspHitFrame " << hspHitFrame << std::endl;
		}
		if ( tag == "Hsp_identity" )
		{
			hspIdentityGatheringActive = false;
			hspIdentityObtained = true;
			// std::cerr << "hspIdentity " << hspIdentity << std::endl;
		}
		if ( tag == "Hsp_positive" )
		{
			hspPositiveGatheringActive = false;
			hspPositiveObtained = true;
			// std::cerr << "hspPositive " << hspPositive << std::endl;
		}
		if ( tag == "Hsp_gaps" )
		{
			hspGapsGatheringActive = false;
			hspGapsObtained = true;
			// std::cerr << "hspGaps " << hspGaps << std::endl;
		}
		if ( tag == "Hsp_align-len" )
		{
			hspAlignLenGatheringActive = false;
			hspAlignLenObtained = true;
			// std::cerr << "hspAlignLen " << hspAlignLen << std::endl;
		}
		if ( tag == "Hsp_qseq" )
		{
			hspQSeqGatheringActive = false;
			hspQSeqObtained = true;
			
			#if 0
			std::cerr << "hspQSeq " << hspQSeq << " " << hspQSeq.size() << std::endl;

			std::map<std::string,std::string>::const_iterator ita = queries.find(readName);
			if ( ita != queries.end() )
			{
				std::cerr << "Found it." << std::endl;
				uint64_t const offset = atoi(hspQueryFrom.c_str());
				std::string const & query = ita->second;
				uint64_t j = 0;
				
				for ( uint64_t i = 0; i < hspQSeq.size(); ++i )
				{
					if ( hspQSeq[i] == '-' )
						std::cerr << "(-,-)";
					else
					{
						std::cerr << "(" << hspQSeq[i] << "," << query[offset-1+j++] << ")";
					}
				}
				std::cerr << std::endl;
			}
			#endif
		}
		if ( tag == "Hsp_hseq" )
		{
			hspHSeqGatheringActive = false;
			hspHSeqObtained = true;
			
			#if 0
			std::cerr << "hspHSeq " << hspHSeq << " " << hspHSeq.size() << std::endl;

			std::map<std::string,std::string>::const_iterator ita = ref.find(hitDef);
			if ( ita != ref.end() )
			{
				std::cerr << "Found it." << std::endl;
				uint64_t const offset = atoi(hspHitFrom.c_str());
				std::string const & refseq = ita->second;
				uint64_t j = 0;
				
				for ( uint64_t i = 0; i < hspHSeq.size(); ++i )
				{
					if ( hspHSeq[i] == '-' )
						std::cerr << "(-,-)";
					else
					{
						std::cerr << "(" << hspHSeq[i] << "," << refseq[offset-1+j++] << ")";
					}
				}
				std::cerr << std::endl;
			}
			#endif

		}
		// std::cerr << "end of element " << utf8transcoder.transcodeStringToUtf8(xname) << std::endl;	
	}

	virtual void warning(const xercesc::SAXParseException& toCatch)
	{
		std::string const msg = utf8transcoder.transcodeStringToUtf8(toCatch.getMessage());
		std::cerr << "[W] " << msg << std::endl;
	}

	virtual void error(const xercesc::SAXParseException& toCatch)
	{
		std::string const msg = utf8transcoder.transcodeStringToUtf8(toCatch.getMessage());
		libmaus::exception::LibMausException lme;
		lme.getStream() << "[E] XML parsing error: " << msg << std::endl;
		lme.finish();
		throw lme;
	}

	virtual void fatalError(const xercesc::SAXParseException& toCatch)
	{
		std::string const msg = utf8transcoder.transcodeStringToUtf8(toCatch.getMessage());
		libmaus::exception::LibMausException lme;
		lme.getStream() << "[E] XML parsing error: " << msg << std::endl;
		lme.finish();
		throw lme;	
	}

	virtual void resetErrors()
	{
	
	}
};

#include <xercesc/parsers/SAXParser.hpp>

#include <libmaus/fastx/FastAReader.hpp>

void loadFastAFile(std::string const & filename, std::map<std::string,std::string> & M, std::vector< std::pair<std::string,uint64_t> > & meta, std::map<std::string,uint64_t> & nametoid)
{
	libmaus::fastx::FastAReader fain(filename);
	libmaus::fastx::FastAReader::pattern_type pattern;
	
	while ( fain.getNextPatternUnlocked(pattern) )
	{
		M[pattern.sid] = pattern.spattern;	
		meta.push_back(std::pair<std::string,uint64_t>(pattern.sid,pattern.spattern.size()));
		uint64_t const id = nametoid.size();
		nametoid[pattern.sid] = id;
	}
}

std::string stripAfterSpace(std::string const & s)
{
	uint64_t firstspace = s.size();
	
	for ( uint64_t i = 0; i < s.size(); ++i )
		if ( isspace(s[i]) )
		{
			firstspace = i;
			break;
		}
		
	return s.substr(0,firstspace);
}

int main(int argc, char * argv[])
{
	int ret = EXIT_SUCCESS;
	bool xercesInitComplete = false;
	
	if ( ret == EXIT_SUCCESS )
		try
		{
			xercesc::XMLPlatformUtils::Initialize();
			xercesInitComplete = true;
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
			ret = EXIT_FAILURE;
		}

	if ( ret == EXIT_SUCCESS )
		try
		{
			libmaus::util::ArgInfo const arginfo(argc,argv);
			double const hitfrac = arginfo.getValue<double>("hitfrac",0.8);
			std::string const reffn = arginfo.restargs.at(0);
			std::string const queriesfn = arginfo.restargs.at(1);
			
			libmaus::util::unique_ptr< std::vector<libmaus::bambam::CramRange> >::type Pranges;
			std::vector<libmaus::bambam::CramRange> * ranges = 0;
			
			if ( arginfo.hasArg("range") )
			{
				libmaus::util::unique_ptr< std::vector<libmaus::bambam::CramRange> >::type Tranges(
					new std::vector<libmaus::bambam::CramRange>
				);
				Pranges = UNIQUE_PTR_MOVE(Tranges);
				
				Pranges->push_back(
					libmaus::bambam::CramRange(
						arginfo.getUnparsedValue("range",std::string())
					)
				);
				
				ranges = Pranges.get();
			}
			
			std::map<std::string,std::string> ref;
			std::vector< std::pair<std::string,uint64_t> > refmeta;
			std::map<std::string,uint64_t> refnametoid;
			loadFastAFile(reffn, ref, refmeta, refnametoid);

			std::map<std::string,std::string> queries;
			std::vector< std::pair<std::string,uint64_t> > queriesmeta;
			std::map<std::string,uint64_t> queriesnametoid;
			loadFastAFile(queriesfn, queries, queriesmeta, queriesnametoid);

			std::ostringstream headerostr;
			headerostr << "@HD\tVN:1.4\tSO:unknown\n";
			for ( uint64_t i = 0; i < refmeta.size(); ++i )
				headerostr << "@SQ\tSN:" << stripAfterSpace(refmeta[i].first) << "\tLN:" << refmeta[i].second << std::endl;
			headerostr 
				<< "@PG"<< "\t" 
				<< "ID:" << "blastnxmltobam" << "\t" 
				<< "PN:" << "blastnxmltobam" << "\t"
				<< "CL:" << arginfo.commandline << "\t"
				<< "VN:" << std::string(PACKAGE_VERSION)
				<< std::endl;

			::libmaus::bambam::BamHeader bamheader(headerostr.str());

			std::cerr << bamheader.text;
			libmaus::bambam::BamWriter writer(std::cout,bamheader);

			XercesUtf8Transcoder transc;
			StdISOInputSource in(std::cin);
			xercesc::SAXParser saxparser;
			saxparser.setValidationScheme(xercesc::SAXParser::Val_Never);
			saxparser.setLoadExternalDTD(false);
			BlastNDocumentHandler blasthandler(ref,queries,refnametoid,queriesnametoid,writer,hitfrac,ranges);
			saxparser.setDocumentHandler(&blasthandler);
			saxparser.setErrorHandler(&blasthandler);
			saxparser.parse(in);
			saxparser.setDocumentHandler(0);                      
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
			return EXIT_FAILURE;
		}

	if ( xercesInitComplete )
		try
		{
			xercesc::XMLPlatformUtils::Terminate();
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
			ret = EXIT_FAILURE;
		}
		
	return ret;
}
