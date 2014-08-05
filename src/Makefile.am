ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS=biobambam/BamBamConfig.hpp biobambam/Licensing.hpp \
	biobambam/ClipAdapters.hpp biobambam/AttachRank.hpp biobambam/ResetAlignment.hpp \
	biobambam/Split12.hpp biobambam/Strip12.hpp \
	biobambam/ClipReinsert.hpp biobambam/zzToName.hpp \
	biobambam/KmerPoisson.hpp

MANPAGES = programs/bamtofastq.1 programs/bamsort.1 programs/bammarkduplicates.1 programs/bamcollate.1 \
	programs/bammaskflags.1 programs/bamrecompress.1 programs/bamadapterfind.1 \
	programs/bamfilteraux.1 programs/bamauxsort.1 programs/bamadapterclip.1 \
	programs/bamcollate2.1 programs/bam12auxmerge.1 programs/bam12split.1 \
	programs/bam12strip.1 programs/bamreset.1 programs/bamrank.1 \
	programs/bamclipreinsert.1 programs/bamzztoname.1 \
	programs/bamindex.1 programs/bamfilterflags.1 \
	programs/bamfilterrg.1 programs/bamdownsamplerandom.1 \
	programs/bammarkduplicates2.1 \
	programs/bamcat.1 programs/bammerge.1 programs/bamsplit.1 \
	programs/bamsplitdiv.1 programs/bamchecksort.1 programs/fastqtobam.1 \
	programs/bamfixmateinformation.1 programs/normalisefasta.1 \
	programs/bamfilterheader.1 \
	programs/bamfilterheader2.1 \
	programs/bammdnm.1 \
	programs/bamseqchksum.1 programs/bamvalidate.1 \
	programs/bamflagsplit.1 programs/bamintervalcomment.1 \
	programs/bamstreamingmarkduplicates.1

man_MANS = ${MANPAGES}

EXTRA_DIST = ${MANPAGES}

bin_PROGRAMS = bamtofastq bammarkduplicates bamsort bamcollate bammaskflags bamrecompress \
	bamadapterfind \
	bamfilteraux \
	bamauxsort \
	bamadapterclip \
	bamcollate2 \
	bam12auxmerge \
	bam12split \
	bam12strip \
	bamreset \
	bamrank \
	bamclipreinsert \
	bamzztoname \
	kmerprob \
	bamindex \
	bamfilterflags \
	bamfilterrg \
	bamdownsamplerandom \
	bammarkduplicates2 \
	bamcat \
	bammerge \
	bamsplit \
	bamsplitdiv \
	bamchecksort \
	fastqtobam \
	bamfixmateinformation \
	bamseqchksum \
	normalisefasta \
	bamfilterheader \
	bamfilterheader2 \
	bammdnm @BAMPARSORTINSTEXP@ \
	bammapdist \
	bamvalidate \
	bamflagsplit \
	bamintervalcomment \
	bamintervalcommenthist \
	bamstreamingmarkduplicates \
	@BAMALIGNMENTOFFSETSINSTEXP@ @BLASTXMLTOBAMINSTEXP@

noinst_PROGRAMS = bamfilter bamfixmatecoordinates bamfixmatecoordinatesnamesorted bamcheckalignments bamtoname \
	bamdisthist bamrefdepth fastabgzfextract bamfilternames @BAMREFDEPTHPEAKS@ \
	@BAMPARSORTNOINSTEXP@ bamheap bamfrontback bamclipextract \
	@BAMALIGNMENTOFFSETSNOINSTEXP@ @BLASTXMLTOBAMNOINSTEXP@ bamrandomtag
	
EXTRA_PROGRAMS = bamrefdepthpeaks bamparsort blastnxmltobam bamalignmentoffsets

bamcollate_SOURCES = programs/bamcollate.cpp biobambam/Licensing.cpp
bamcollate_LDADD = ${LIBMAUSLIBS}
bamcollate_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcollate_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilter_SOURCES = programs/bamfilter.cpp biobambam/Licensing.cpp
bamfilter_LDADD = ${LIBMAUSLIBS}
bamfilter_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilter_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfixmatecoordinates_SOURCES = programs/bamfixmatecoordinates.cpp biobambam/Licensing.cpp
bamfixmatecoordinates_LDADD = ${LIBMAUSLIBS}
bamfixmatecoordinates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfixmatecoordinates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfixmatecoordinatesnamesorted_SOURCES = programs/bamfixmatecoordinatesnamesorted.cpp biobambam/Licensing.cpp
bamfixmatecoordinatesnamesorted_LDADD = ${LIBMAUSLIBS}
bamfixmatecoordinatesnamesorted_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfixmatecoordinatesnamesorted_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bammarkduplicates_SOURCES = programs/bammarkduplicates.cpp biobambam/Licensing.cpp
bammarkduplicates_LDADD = ${LIBMAUSLIBS}
bammarkduplicates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammarkduplicates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamstreamingmarkduplicates_SOURCES = programs/bamstreamingmarkduplicates.cpp biobambam/Licensing.cpp
bamstreamingmarkduplicates_LDADD = ${LIBMAUSLIBS}
bamstreamingmarkduplicates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamstreamingmarkduplicates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bammarkduplicates2_SOURCES = programs/bammarkduplicates2.cpp biobambam/Licensing.cpp
bammarkduplicates2_LDADD = ${LIBMAUSLIBS}
bammarkduplicates2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammarkduplicates2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bammaskflags_SOURCES = programs/bammaskflags.cpp biobambam/Licensing.cpp
bammaskflags_LDADD = ${LIBMAUSLIBS}
bammaskflags_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammaskflags_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamsort_SOURCES = programs/bamsort.cpp biobambam/Licensing.cpp
bamsort_LDADD = ${LIBMAUSLIBS}
bamsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamtofastq_SOURCES = programs/bamtofastq.cpp biobambam/Licensing.cpp
bamtofastq_LDADD = ${LIBMAUSLIBS}
bamtofastq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS}
bamtofastq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamcheckalignments_SOURCES = programs/bamcheckalignments.cpp biobambam/Licensing.cpp
bamcheckalignments_LDADD = ${LIBMAUSLIBS}
bamcheckalignments_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcheckalignments_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamtoname_SOURCES = programs/bamtoname.cpp biobambam/Licensing.cpp
bamtoname_LDADD = ${LIBMAUSLIBS}
bamtoname_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS}
bamtoname_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamdisthist_SOURCES = programs/bamdisthist.cpp biobambam/Licensing.cpp
bamdisthist_LDADD = ${LIBMAUSLIBS}
bamdisthist_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS}
bamdisthist_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamrecompress_SOURCES = programs/bamrecompress.cpp biobambam/Licensing.cpp
bamrecompress_LDADD = ${LIBMAUSLIBS}
bamrecompress_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamrecompress_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamadapterfind_SOURCES = programs/bamadapterfind.cpp biobambam/Licensing.cpp \
	biobambam/ClipAdapters.cpp biobambam/KmerPoisson.cpp
bamadapterfind_LDADD = ${LIBMAUSLIBS}
bamadapterfind_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamadapterfind_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilteraux_SOURCES = programs/bamfilteraux.cpp biobambam/Licensing.cpp
bamfilteraux_LDADD = ${LIBMAUSLIBS}
bamfilteraux_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilteraux_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamauxsort_SOURCES = programs/bamauxsort.cpp biobambam/Licensing.cpp
bamauxsort_LDADD = ${LIBMAUSLIBS}
bamauxsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamauxsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamadapterclip_SOURCES = programs/bamadapterclip.cpp biobambam/Licensing.cpp \
	biobambam/ClipAdapters.cpp
bamadapterclip_LDADD = ${LIBMAUSLIBS}
bamadapterclip_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamadapterclip_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamcollate2_SOURCES = programs/bamcollate2.cpp biobambam/Licensing.cpp biobambam/AttachRank.cpp \
	biobambam/ResetAlignment.cpp
bamcollate2_LDADD = ${LIBMAUSLIBS}
bamcollate2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcollate2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bam12auxmerge_SOURCES = programs/bam12auxmerge.cpp biobambam/Licensing.cpp biobambam/Split12.cpp \
	biobambam/Strip12.cpp biobambam/ClipReinsert.cpp biobambam/zzToName.cpp
bam12auxmerge_LDADD = ${LIBMAUSLIBS}
bam12auxmerge_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bam12auxmerge_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bam12split_SOURCES = programs/bam12split.cpp biobambam/Licensing.cpp biobambam/Split12.cpp
bam12split_LDADD = ${LIBMAUSLIBS}
bam12split_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bam12split_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bam12strip_SOURCES = programs/bam12strip.cpp biobambam/Licensing.cpp biobambam/Strip12.cpp
bam12strip_LDADD = ${LIBMAUSLIBS}
bam12strip_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bam12strip_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamreset_SOURCES = programs/bamreset.cpp biobambam/Licensing.cpp biobambam/ResetAlignment.cpp
bamreset_LDADD = ${LIBMAUSLIBS}
bamreset_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamreset_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamrank_SOURCES = programs/bamrank.cpp biobambam/Licensing.cpp biobambam/AttachRank.cpp
bamrank_LDADD = ${LIBMAUSLIBS}
bamrank_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamrank_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamclipreinsert_SOURCES = programs/bamclipreinsert.cpp biobambam/Licensing.cpp \
	biobambam/ClipReinsert.cpp
bamclipreinsert_LDADD = ${LIBMAUSLIBS}
bamclipreinsert_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamclipreinsert_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamzztoname_SOURCES = programs/bamzztoname.cpp biobambam/Licensing.cpp biobambam/zzToName.cpp
bamzztoname_LDADD = ${LIBMAUSLIBS}
bamzztoname_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamzztoname_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

kmerprob_SOURCES = programs/kmerprob.cpp biobambam/Licensing.cpp biobambam/KmerPoisson.cpp
kmerprob_LDADD = ${LIBMAUSLIBS}
kmerprob_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
kmerprob_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamindex_SOURCES = programs/bamindex.cpp biobambam/Licensing.cpp
bamindex_LDADD = ${LIBMAUSLIBS}
bamindex_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamindex_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilterflags_SOURCES = programs/bamfilterflags.cpp biobambam/Licensing.cpp
bamfilterflags_LDADD = ${LIBMAUSLIBS}
bamfilterflags_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilterflags_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilterrg_SOURCES = programs/bamfilterrg.cpp biobambam/Licensing.cpp
bamfilterrg_LDADD = ${LIBMAUSLIBS}
bamfilterrg_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilterrg_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamdownsamplerandom_SOURCES = programs/bamdownsamplerandom.cpp biobambam/Licensing.cpp
bamdownsamplerandom_LDADD = ${LIBMAUSLIBS}
bamdownsamplerandom_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamdownsamplerandom_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamcat_SOURCES = programs/bamcat.cpp biobambam/Licensing.cpp
bamcat_LDADD = ${LIBMAUSLIBS}
bamcat_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcat_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bammerge_SOURCES = programs/bammerge.cpp biobambam/Licensing.cpp
bammerge_LDADD = ${LIBMAUSLIBS}
bammerge_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammerge_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamsplit_SOURCES = programs/bamsplit.cpp biobambam/Licensing.cpp
bamsplit_LDADD = ${LIBMAUSLIBS}
bamsplit_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamsplit_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamsplitdiv_SOURCES = programs/bamsplitdiv.cpp biobambam/Licensing.cpp
bamsplitdiv_LDADD = ${LIBMAUSLIBS}
bamsplitdiv_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamsplitdiv_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamchecksort_SOURCES = programs/bamchecksort.cpp biobambam/Licensing.cpp
bamchecksort_LDADD = ${LIBMAUSLIBS}
bamchecksort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamchecksort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamrefdepth_SOURCES = programs/bamrefdepth.cpp biobambam/Licensing.cpp
bamrefdepth_LDADD = ${LIBMAUSLIBS}
bamrefdepth_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamrefdepth_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamrefdepthpeaks_SOURCES = programs/bamrefdepthpeaks.cpp biobambam/Licensing.cpp
bamrefdepthpeaks_LDADD = ${LIBMAUSLIBS}
bamrefdepthpeaks_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamrefdepthpeaks_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

fastqtobam_SOURCES = programs/fastqtobam.cpp biobambam/Licensing.cpp
fastqtobam_LDADD = ${LIBMAUSLIBS}
fastqtobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
fastqtobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfixmateinformation_SOURCES = programs/bamfixmateinformation.cpp biobambam/Licensing.cpp
bamfixmateinformation_LDADD = ${LIBMAUSLIBS}
bamfixmateinformation_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfixmateinformation_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamseqchksum_SOURCES = programs/bamseqchksum.cpp biobambam/Licensing.cpp
bamseqchksum_LDADD = ${LIBMAUSLIBS}
bamseqchksum_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamseqchksum_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

normalisefasta_SOURCES = programs/normalisefasta.cpp biobambam/Licensing.cpp
normalisefasta_LDADD = ${LIBMAUSLIBS}
normalisefasta_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
normalisefasta_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilterheader_SOURCES = programs/bamfilterheader.cpp biobambam/Licensing.cpp
bamfilterheader_LDADD = ${LIBMAUSLIBS}
bamfilterheader_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilterheader_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilterheader2_SOURCES = programs/bamfilterheader2.cpp biobambam/Licensing.cpp
bamfilterheader2_LDADD = ${LIBMAUSLIBS}
bamfilterheader2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilterheader2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bammdnm_SOURCES = programs/bammdnm.cpp biobambam/Licensing.cpp
bammdnm_LDADD = ${LIBMAUSLIBS}
bammdnm_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammdnm_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

fastabgzfextract_SOURCES = programs/fastabgzfextract.cpp biobambam/Licensing.cpp
fastabgzfextract_LDADD = ${LIBMAUSLIBS}
fastabgzfextract_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
fastabgzfextract_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfilternames_SOURCES = programs/bamfilternames.cpp biobambam/Licensing.cpp
bamfilternames_LDADD = ${LIBMAUSLIBS}
bamfilternames_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilternames_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamparsort_SOURCES = programs/bamparsort.cpp biobambam/Licensing.cpp
bamparsort_LDADD = ${LIBMAUSLIBS}
bamparsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamparsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bammapdist_SOURCES = programs/bammapdist.cpp biobambam/Licensing.cpp
bammapdist_LDADD = ${LIBMAUSLIBS}
bammapdist_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammapdist_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamheap_SOURCES = programs/bamheap.cpp biobambam/Licensing.cpp
bamheap_LDADD = ${LIBMAUSLIBS}
bamheap_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamheap_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamfrontback_SOURCES = programs/bamfrontback.cpp biobambam/Licensing.cpp biobambam/ResetAlignment.cpp
bamfrontback_LDADD = ${LIBMAUSLIBS}
bamfrontback_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfrontback_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamclipextract_SOURCES = programs/bamclipextract.cpp biobambam/Licensing.cpp biobambam/ResetAlignment.cpp
bamclipextract_LDADD = ${LIBMAUSLIBS}
bamclipextract_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamclipextract_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

blastnxmltobam_SOURCES = programs/blastnxmltobam.cpp biobambam/Licensing.cpp
blastnxmltobam_LDADD = ${LIBMAUSLIBS} @xerces_c_LIBS@
blastnxmltobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
blastnxmltobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} @xerces_c_CFLAGS@

bamalignmentoffsets_SOURCES = programs/bamalignmentoffsets.cpp biobambam/Licensing.cpp
bamalignmentoffsets_LDADD = ${LIBMAUSLIBS}
bamalignmentoffsets_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamalignmentoffsets_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamvalidate_SOURCES = programs/bamvalidate.cpp biobambam/Licensing.cpp
bamvalidate_LDADD = ${LIBMAUSLIBS}
bamvalidate_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamvalidate_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamflagsplit_SOURCES = programs/bamflagsplit.cpp biobambam/Licensing.cpp
bamflagsplit_LDADD = ${LIBMAUSLIBS}
bamflagsplit_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamflagsplit_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamintervalcomment_SOURCES = programs/bamintervalcomment.cpp biobambam/Licensing.cpp
bamintervalcomment_LDADD = ${LIBMAUSLIBS}
bamintervalcomment_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamintervalcomment_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamintervalcommenthist_SOURCES = programs/bamintervalcommenthist.cpp biobambam/Licensing.cpp
bamintervalcommenthist_LDADD = ${LIBMAUSLIBS}
bamintervalcommenthist_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamintervalcommenthist_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamrandomtag_SOURCES = programs/bamrandomtag.cpp biobambam/Licensing.cpp
bamrandomtag_LDADD = ${LIBMAUSLIBS}
bamrandomtag_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamrandomtag_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
