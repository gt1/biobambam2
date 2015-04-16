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

#include <liftingwavelettransform/LiftingWaveletTransform.hpp>

#include <biobambam/Licensing.hpp>

static int getDefaultVerbose() { return 0; }

#include <algorithm>

/**
 * compute a single point for a convolution
 **/
template<typename accessor_type, typename filter_iterator> 
typename accessor_type::value_type interpolatingConvolveSingle(
	uint64_t const i, accessor_type const & M, 
	filter_iterator const F, int64_t const f, int64_t const sub,
	int64_t const bs, float const s
)
{
	typedef typename accessor_type::value_type value_type;
	value_type v = value_type();

	if ( f )
	{
		int64_t o = sub*bs;
		
		float vj = M(i,o);
		
		for ( int64_t j = 0; j < f; ++j )
		{
			o += bs;
			float const vj1 = M(i,o);
			float const vs = (vj+vj1)*.5f;
			v += vs * F[j];
			
			vj = vj1;
		}
	}

	return v*s;
}

/**
 * get support radius for filter of length f with block size bs
 **/
uint64_t getSupportRadius(uint64_t const f, uint64_t const bs)
{
	return f * bs;
}

/**
 * compute convolution and store result in a different vector
 **/
template<typename data_iterator, typename filter_iterator>
std::vector< typename std::iterator_traits<data_iterator>::value_type >
	interpolatingConvolve(data_iterator const A, uint64_t const n, filter_iterator const F, uint64_t const f, uint64_t const bs, float const s = 1.0f/::std::sqrt(2))
{
	typedef typename std::iterator_traits<data_iterator>::value_type value_type;
	std::vector< value_type > R(n);
	LiftingWaveletTransform::MirrorAccessor<data_iterator> const M(A,n);
	int64_t const sub = -(f/2);
	
	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for ( uint64_t i = 0; i < n; ++i )
		R[i] = interpolatingConvolveSingle(i,M,F,f,sub,bs,s);

	return R;
}

/**
 * compute convolution using as little additional space as possible
 **/
template<typename data_iterator, typename filter_iterator>
void
	interpolatingConvolveInPlace(data_iterator const A, uint64_t const n, filter_iterator const F, int64_t const f, uint64_t const bs, float const s = 1.0f/::std::sqrt(2))
{
	typedef typename std::iterator_traits<data_iterator>::value_type value_type;
	LiftingWaveletTransform::MirrorAccessor<data_iterator> const M(A,n);
	int64_t const sub = -(f/2);
	
	uint64_t const suprad = getSupportRadius(f,bs);

	if ( n < suprad+1 )
	{
		// use external array if input is too small
		std::vector< value_type > R = interpolatingConvolve(A,n,F,f,bs,s);
		std::copy(R.begin(),R.end(),A);
	}
	else
	{
		// precompute first suprad+1 elements and store them in additional space
		std::vector< value_type > V0(suprad+1);
		for ( uint64_t i = 0; i < suprad+1; ++i )
			V0[i] = interpolatingConvolveSingle(i,M,F,f,sub,bs,s);
		// compute rest of elements and store them at the front of the array
		for ( uint64_t i = 0; i < n-(suprad+1); ++i )
			A[i] = interpolatingConvolveSingle(i+(suprad+1),M,F,f,sub,bs,s);
		// copy rest of elements to the back of the array
		std::copy_backward(A,A+(n-(suprad+1)),A+n);
		// copy first elements into their place
		std::copy(V0.begin(),V0.end(),A);
	}
}

struct PeakInfo
{
	uint64_t bs;
	int64_t i;
	float hpre;
	float hpost;
	
	PeakInfo() : bs(0), i(0), hpre(0), hpost(0) {}
	PeakInfo(uint64_t const rbs,
		int64_t const ri,
		float const rhpre,
		float const rhpost
	) : bs(rbs), i(ri), hpre(rhpre), hpost(rhpost)
	{
	
	}
	
	bool operator<(PeakInfo const & o) const
	{
		return i < o.i;
	}
	
	int64_t getLeft() const
	{
		return static_cast<int64_t>(i) - static_cast<int64_t>(bs);
	}

	int64_t getRight() const
	{
		return static_cast<int64_t>(i) + static_cast<int64_t>(bs);
	}
};

std::ostream & operator<<(std::ostream & out, PeakInfo const & O)
{
	out << "PeakInfo(bs=" << O.bs << ",i=" << O.i << ",hpre=" << O.hpre << ",hpost=" << O.hpost << ")";
	return out;
}

struct PeakInfoBsComparator
{
	bool operator()(PeakInfo const & A, PeakInfo const & B) const
	{
		if ( A.bs != B.bs )
			return A.bs < B.bs;
		else
			return A.i < B.i;
	}
};

struct PeakInfoLeftComparator
{
	bool operator()(PeakInfo const & A, PeakInfo const & B) const
	{
		int64_t const aleft = A.getLeft();
		int64_t const bleft = B.getLeft();
		
		if ( aleft != bleft )
		{
			return aleft < bleft;
		}
		else
		{
			return A.bs > B.bs;
		}
	}
};

std::vector< PeakInfo  > analyse(std::deque<float> const & Q, std::string const & name, float const peakthres, uint64_t const minpeakwidth, uint64_t const maxpeakwidth)
{
	std::cerr << name << "\t" << Q.size() << std::endl;

	// biorthogonal 3.1 scaling function coefficients
	float const bior_3_1_reconst_low [] = { 0.1767766953, 0.5303300859, 0.5303300859, 0.1767766953 };
	// biorthogonal 3.1 wavelet coefficients
	float const bior_3_1_reconst_high[] = { -0.3535533906,  -1.0606601718, 1.0606601718, 0.3535533906 };

	// length of filters	
	int64_t const f = sizeof(bior_3_1_reconst_low)/sizeof(bior_3_1_reconst_low[0]);
	// scaling factor for filtering
	float const s = 1/sqrt(2);
	// length of input sequence
	uint64_t const n = Q.size();

	// low pass
	std::vector< float > L(Q.begin(),Q.end());
	
	std::vector < PeakInfo > peaks;
	
	#if 0
	std::map< uint64_t, std::vector<float> > fullHMap;
	#endif
			
	for ( uint64_t bs = 1; bs <= maxpeakwidth; bs *= 2 )
	{
		// float const scalefactor = s * bs;
		float const scalefactor = s;
		// std::cerr << "bs=" << bs << " Q.size()=" << Q.size() << std::endl;
	
		if ( bs >= minpeakwidth )
		{
			LiftingWaveletTransform::MirrorAccessor< std::vector<float>::const_iterator > const M(L.begin(),n);
			float hpre = n ? interpolatingConvolveSingle(0,M,&bior_3_1_reconst_high[0],f,-f/2,bs,scalefactor) : 0;
			std::deque<float> hvec;
			hvec.push_back(hpre);
	
			#if 0
			std::vector<float> & fullH = fullHMap[bs];
			fullH.push_back(hpre);
			#endif
		
			for ( int64_t i = 1; i < static_cast<int64_t>(Q.size()); ++i )
			{
				float const hnext = interpolatingConvolveSingle(i,M,&bior_3_1_reconst_high[0],f,-f/2,bs,scalefactor);
				hvec.push_back(hnext);
			
				#if 0
				fullH.push_back(hnext);
				#endif
			
				int64_t const backpos = static_cast<int64_t>(hvec.size())-1;
				int64_t const difpos = backpos - static_cast<int64_t>(bs);
				
				if ( 
					difpos >= 0 
					&&
					hvec[difpos] > peakthres && hvec[backpos] < -peakthres
				)
				{
					peaks.push_back(PeakInfo(bs,i-bs/2,hvec[difpos],hvec[backpos]));
					hvec.pop_front();
				}
				
				#if 0
				if ( hpre > peakthres && hnext < -peakthres )
					peaks.push_back(PeakInfo(bs,i,hpre,hnext));					
			
				// hpre = hnext;
				#endif
			}
		}

		// compute low pass on L
		interpolatingConvolveInPlace(L.begin(),L.size(),&bior_3_1_reconst_low[0],f,bs,s);
	}

	std::sort(peaks.begin(),peaks.end(),PeakInfoBsComparator());
	// std::vector<PeakInfo> npeaks;
	uint64_t optr = 0;
	
	uint64_t low = 0;
	while ( low != peaks.size() )
	{
		uint64_t high = low;
		while ( high != peaks.size() && peaks[high].bs == peaks[low].bs )
			++high;
			
		uint64_t llow = low;
		while ( llow != high )
		{
			uint64_t lhigh = llow+1;

			while ( 
				lhigh != high 
				&& 
				peaks[lhigh].i - peaks[llow].i == static_cast<int64_t>(lhigh-llow)
			)
			{
				++lhigh;
			}
						
			#if 0
			uint64_t const mid = peaks[(llow+lhigh)/2].i;
			
			std::cerr << "===" << mid << "===" << std::endl;
			for ( uint64_t i = llow; i != lhigh; ++i )
				std::cerr << peaks[i] << std::endl;
			std::cerr << "===" << std::endl;
			#endif
			
			// npeaks.push_back(peaks[(llow+lhigh)/2]);
			peaks[optr++] = peaks[(llow+lhigh)/2];
			
			llow = lhigh;
		}
			
		low = high;
	}

	#if 0
	peaks = npeaks;
	npeaks.resize(0);
	#endif

	peaks.resize(optr);

	/* check for containment of smaller peaks in larger peaks (drop small contained ones) */
	std::sort(peaks.begin(),peaks.end(),PeakInfoLeftComparator());
	uint64_t checkptr = 0;
	optr = 0;
	for ( uint64_t i = 0; i < peaks.size(); ++i )
	{
		// advance checkptr if intervals do not overlap
		while ( checkptr < i && peaks[i].getLeft() > peaks[checkptr].getRight() )
			++checkptr;
		
		// i is the new check interval, keep it
		if ( i == checkptr )
		{
			peaks[optr++] = peaks[i];
			#if 0
			std::cerr << "keeping[1] " << peaks[i] << std::endl;
			#endif
		}
		else
		{
			// no overlap, make i the new check interval
			if ( peaks[i].getRight() > peaks[checkptr].getRight() )
			{
				peaks[optr++] = peaks[i];
				checkptr = i;
				#if 0
				std::cerr << "keeping[2] " << peaks[i] << std::endl;
				#endif
			}
			else // if ( peaks[i].getRight() <= peaks[checkptr].getRight() )
			{
				#if 0
				std::cerr << "dropping " << peaks[i] << " contained in " << peaks[checkptr] << std::endl;
				#endif
			}
		}
	}
	peaks.resize(optr);

	// sort remaining peaks by centre
	std::sort(peaks.begin(),peaks.end());

	// int64_t const kpre = 10;
	
	for ( uint64_t z = 0; z < peaks.size(); ++z )
	{
		// int64_t const i = peaks[z].i;
		// uint64_t const bs = peaks[z].bs;
		// float const hpre = peaks[z].hpre;
		// float const hpost = peaks[z].hpost;
		// int64_t const k = kpre + bs;
		
		std::ostringstream ostr;
		ostr << " " << peaks[z] << " ";
		std::string const ostrstr = ostr.str();
		std::deque<char> D(ostrstr.begin(),ostrstr.end());
		while ( D.size() < 80 )
			if ( D.size() & 1 )
				D.push_back('*');
			else
				D.push_front('*');
				
		std::cerr << std::string(D.begin(),D.end()) << std::endl;

		#if 0
		for ( int64_t j = std::max(static_cast<int64_t>(0),i-k); j <= std::min(static_cast<int64_t>(L.size()),i+k); ++j )
			std::cerr << "j=" << j  << "\tQ=" << Q[j  ] << std::endl;				
		#endif
	}
	
	return peaks;
}

void generateGPL(std::deque<float> const & Q, std::vector< PeakInfo > const & peaks, std::string const & seqname)
{	
	#if 0
	std::set<int64_t> peakcoord;
	for ( uint64_t i = 0; i < peaks.size(); ++i )
	{
		int64_t centre = peaks[i].i;
		int64_t low = centre - peaks[i].bs;
		int64_t high = centre + peaks[i].bs;
		
		for ( int64_t j = low; j <= high; ++j )
			peakcoord.insert(j);
	}
	#endif
	
	// uint64_t const numblocks = 1024;
	// uint64_t blocksize = (Q.size() + numblocks-1)/numblocks;
	
	std::ostringstream depthstr;
	// depthstr << "plot_" << header.getRefIDName(previd) << ".gpl";
	depthstr << "plot_" << seqname << ".gpl";
	std::ostringstream depthpeakstr;
	depthpeakstr << "plot_" << seqname << "_peaks.gpl";
	
	libmaus::aio::CheckedOutputStream depthostr(depthstr.str());
	libmaus::aio::CheckedOutputStream depthpeakostr(depthpeakstr.str());
	
	for ( uint64_t i = 0; i < peaks.size(); ++i )
	{
		int64_t centre = peaks[i].i;
		int64_t low = centre - peaks[i].bs;
		int64_t high = centre + peaks[i].bs;
		float maxval = 0;
		
		for ( int64_t j = low; j <= high; ++j )
			maxval = std::max(maxval,Q[j]);

		for ( int64_t j = low; j <= high; ++j )
			depthpeakostr << j << "\t" << maxval << std::endl;
	}

	for ( uint64_t i = 0; i < Q.size(); i ++ )
	{
		depthostr << (i) << "\t" << Q[i] << "\n";
	}
	
	#if 0
	for ( uint64_t i = 0; i < Q.size(); i += blocksize )
	{
		float avg = 0;
		#if 0
		bool mark = false;
		#endif
		for ( uint64_t j = i; j < std::min(i + blocksize,static_cast<uint64_t>(Q.size())); ++j )
		{
			avg += Q[j];
			#if 0
			if ( peakcoord.find(j) != peakcoord.end() )
				mark = true;
			#endif
		}
		
		if ( avg > 0 )
		{
			depthostr << (i+blocksize/2) << "\t" << avg/blocksize << "\n";
			#if 0
			if ( mark )
				depthpeakostr << (i+blocksize/2) << "\t" << avg/blocksize << "\n"; 
			#endif
		}
	}
	#endif
	
	depthostr.flush(); depthostr.close();
	depthpeakostr.flush(); depthpeakostr.close();
}

template<typename container_type>
void ztrim(container_type & Q, uint64_t const keep)
{
	uint64_t frontz = 0;
	while ( frontz < Q.size() && Q[frontz] == 0 )
		++frontz;
	while ( frontz > keep )
	{
		frontz--;
		Q.pop_front();
	}
	
	uint64_t backz = 0;
	while ( backz < Q.size() && Q[Q.size()-backz-1] == 0 )
		backz++;
		
	while ( backz > keep )
	{
		backz--;
		Q.pop_back();
	}
}

int bamrefdepth(libmaus::util::ArgInfo const & arginfo)
{
	int const verbose = arginfo.getValue<int>("verbose",getDefaultVerbose());

	libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type pdec(
	libmaus::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
	libmaus::bambam::BamAlignmentDecoder & bamdec = pdec->getDecoder();
	libmaus::bambam::BamAlignment & algn = bamdec.getAlignment();
	libmaus::bambam::BamHeader const & header = bamdec.getHeader();
	std::string const sortorder = libmaus::bambam::BamHeader::getSortOrderStatic(header.text);
	float const peakthres = arginfo.getValue<float>("peakthres",1.0f);
	uint64_t const minpeakwidth = arginfo.getValueUnsignedNumeric<uint64_t>("minpeakwidth",1);
	uint64_t const maxpeakwidth = arginfo.getValueUnsignedNumeric<uint64_t>("maxpeakwidth",1024);
	
	libmaus::bambam::BamAlignment prevalgn;
	bool hasprev = false;
	uint64_t c = 0;
	
	std::deque<float> Q;
	libmaus::autoarray::AutoArray<libmaus::bambam::cigar_operation> cigop;
	libmaus::autoarray::AutoArray<char> decread;
	int64_t previd = -1;
	
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
			uint64_t const len = header.getRefIDLength(previd);
			while ( Q.size() < len )
				Q.push_back(0);
				
			ztrim(Q,2);
				
			std::vector< PeakInfo > const peaks = analyse(Q,refnames[previd],peakthres,minpeakwidth,maxpeakwidth);
			generateGPL(Q,peaks,header.getRefIDName(previd));
			// std::cerr << refnames[prevalgn.getRefID()] << " size " << Q.size() << std::endl;
			Q.resize(0);
			
			std::cerr << "inter." << std::endl;
		}
		
		if ( algn.isMapped() )
		{
			uint32_t const numcigop = algn.getCigarOperations(cigop);
			uint64_t const readlen = algn.decodeRead(decread);
			int64_t pos = algn.getPos();
			uint64_t readpos = 0;
			
			// std::cerr << "Q.size()=" << Q.size() << std::endl;

			#if 0
			while ( Q.size() && pos > leftpos )
			{
				if ( Q[0] )
					std::cout << refnames[algn.getRefID()] << "\t" << leftpos << "\t" << Q[0] << std::endl;
				Q.pop_front();
				leftpos++;
			}
			#endif
			
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
						for ( uint64_t i = 0; i < static_cast<uint64_t>(cigop[cidx].second); ++i, ++pos, ++readpos )
						{
							while ( pos >= static_cast<int64_t>(Q.size()) )
								Q.push_back(0);
								
							Q[pos] += 1;
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
			
			previd = algn.getRefID();
		}
			
		prevalgn.swap(algn);
		hasprev = true;

		if ( verbose && ( ((++c) & ((1ull<<20)-1)) == 0 ) )
			std::cerr << "[V] " << c << std::endl;
	}

	// std::cerr << refnames[prevalgn.getRefID()] << " size " << Q.size() << std::endl;

	if ( Q.size() )
	{
		//std::cerr << refnames[prevalgn.getRefID()] << " size " << Q.size() << std::endl;
		uint64_t const len = header.getRefIDLength(previd);
		while ( Q.size() < len )
			Q.push_back(0);

		ztrim(Q,2);

		std::vector< PeakInfo > peaks = analyse(Q,refnames[previd],peakthres,minpeakwidth,maxpeakwidth);

		generateGPL(Q,peaks,header.getRefIDName(previd));

		Q.resize(0);
		
		std::cerr << "done." << std::endl;
	}

	if ( verbose )
		std::cerr << "[V] " << c << std::endl;

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	#if 0
	{
	// biorthogonal 3.1 scaling function coefficients
	float const bior_3_1_reconst_low [] = { 0.1767766953, 0.5303300859, 0.5303300859, 0.1767766953 };
	// biorthogonal 3.1 wavelet coefficients
	float const bior_3_1_reconst_high[] = { -0.3535533906,  -1.0606601718, 1.0606601718, 0.3535533906 };

	float A[1841];
	uint64_t const n = sizeof(A)/sizeof(A[0]);
	srand(5);
	
	for ( uint64_t j = 0; j < 1024; ++j )
	{
	for ( uint64_t i = 0; i < n; ++i )
		A[i] = rand() % 256;

	uint64_t const bs = 4;
	std::vector<float> H = interpolatingConvolve(
		&A[0],n,
		&bior_3_1_reconst_high[0],sizeof(bior_3_1_reconst_low)/sizeof(bior_3_1_reconst_high[0]),bs
	);
	interpolatingConvolveInPlace(
		&A[0],n,
		&bior_3_1_reconst_high[0],sizeof(bior_3_1_reconst_low)/sizeof(bior_3_1_reconst_high[0]),bs
	);
	
	for ( uint64_t i = 0; i < sizeof(A)/sizeof(A[0]); ++i )
		if ( A[i] != H[i] )
			std::cerr << i << "\t" << A[i] << "\t" << H[i] << std::endl;
	}

	return 0;
	}
	#endif

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

