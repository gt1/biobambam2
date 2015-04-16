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
	programs/bamstreamingmarkduplicates.1 \
	programs/bamsormadup.1

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
	bamheap2 \
	bamalignfrac \
	bamlastfilter \
	bamfilternames \
	bamsormadup \
	bamexplode \
	fastqtobampar \
	@BAMALIGNMENTOFFSETSINSTEXP@ @BLASTXMLTOBAMINSTEXP@

noinst_PROGRAMS = bamfilter bamfixmatecoordinates bamfixmatecoordinatesnamesorted bamcheckalignments bamtoname \
	bamdisthist bamrefdepth fastabgzfextract @BAMREFDEPTHPEAKS@ \
	@BAMPARSORTNOINSTEXP@ bamheap bamfrontback bamclipextract \
	@BAMALIGNMENTOFFSETSNOINSTEXP@ @BLASTXMLTOBAMNOINSTEXP@ bamrandomtag
	
EXTRA_PROGRAMS = bamrefdepthpeaks bamparsort blastnxmltobam bamalignmentoffsets

bamcollate_SOURCES = programs/bamcollate.cpp biobambam/Licensing.cpp
bamcollate_LDADD = ${LIBMAUS2LIBS}
bamcollate_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamcollate_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilter_SOURCES = programs/bamfilter.cpp biobambam/Licensing.cpp
bamfilter_LDADD = ${LIBMAUS2LIBS}
bamfilter_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilter_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfixmatecoordinates_SOURCES = programs/bamfixmatecoordinates.cpp biobambam/Licensing.cpp
bamfixmatecoordinates_LDADD = ${LIBMAUS2LIBS}
bamfixmatecoordinates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfixmatecoordinates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfixmatecoordinatesnamesorted_SOURCES = programs/bamfixmatecoordinatesnamesorted.cpp biobambam/Licensing.cpp
bamfixmatecoordinatesnamesorted_LDADD = ${LIBMAUS2LIBS}
bamfixmatecoordinatesnamesorted_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfixmatecoordinatesnamesorted_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bammarkduplicates_SOURCES = programs/bammarkduplicates.cpp biobambam/Licensing.cpp
bammarkduplicates_LDADD = ${LIBMAUS2LIBS}
bammarkduplicates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bammarkduplicates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamstreamingmarkduplicates_SOURCES = programs/bamstreamingmarkduplicates.cpp biobambam/Licensing.cpp
bamstreamingmarkduplicates_LDADD = ${LIBMAUS2LIBS}
bamstreamingmarkduplicates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamstreamingmarkduplicates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bammarkduplicates2_SOURCES = programs/bammarkduplicates2.cpp biobambam/Licensing.cpp
bammarkduplicates2_LDADD = ${LIBMAUS2LIBS}
bammarkduplicates2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bammarkduplicates2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bammaskflags_SOURCES = programs/bammaskflags.cpp biobambam/Licensing.cpp
bammaskflags_LDADD = ${LIBMAUS2LIBS}
bammaskflags_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bammaskflags_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamsort_SOURCES = programs/bamsort.cpp biobambam/Licensing.cpp
bamsort_LDADD = ${LIBMAUS2LIBS}
bamsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamtofastq_SOURCES = programs/bamtofastq.cpp biobambam/Licensing.cpp
bamtofastq_LDADD = ${LIBMAUS2LIBS}
bamtofastq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
bamtofastq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamcheckalignments_SOURCES = programs/bamcheckalignments.cpp biobambam/Licensing.cpp
bamcheckalignments_LDADD = ${LIBMAUS2LIBS}
bamcheckalignments_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamcheckalignments_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamtoname_SOURCES = programs/bamtoname.cpp biobambam/Licensing.cpp
bamtoname_LDADD = ${LIBMAUS2LIBS}
bamtoname_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
bamtoname_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamdisthist_SOURCES = programs/bamdisthist.cpp biobambam/Licensing.cpp
bamdisthist_LDADD = ${LIBMAUS2LIBS}
bamdisthist_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
bamdisthist_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamrecompress_SOURCES = programs/bamrecompress.cpp biobambam/Licensing.cpp
bamrecompress_LDADD = ${LIBMAUS2LIBS}
bamrecompress_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamrecompress_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamadapterfind_SOURCES = programs/bamadapterfind.cpp biobambam/Licensing.cpp \
	biobambam/ClipAdapters.cpp biobambam/KmerPoisson.cpp
bamadapterfind_LDADD = ${LIBMAUS2LIBS}
bamadapterfind_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamadapterfind_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilteraux_SOURCES = programs/bamfilteraux.cpp biobambam/Licensing.cpp
bamfilteraux_LDADD = ${LIBMAUS2LIBS}
bamfilteraux_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilteraux_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamauxsort_SOURCES = programs/bamauxsort.cpp biobambam/Licensing.cpp
bamauxsort_LDADD = ${LIBMAUS2LIBS}
bamauxsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamauxsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamadapterclip_SOURCES = programs/bamadapterclip.cpp biobambam/Licensing.cpp \
	biobambam/ClipAdapters.cpp
bamadapterclip_LDADD = ${LIBMAUS2LIBS}
bamadapterclip_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamadapterclip_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamcollate2_SOURCES = programs/bamcollate2.cpp biobambam/Licensing.cpp biobambam/AttachRank.cpp \
	biobambam/ResetAlignment.cpp
bamcollate2_LDADD = ${LIBMAUS2LIBS}
bamcollate2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamcollate2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bam12auxmerge_SOURCES = programs/bam12auxmerge.cpp biobambam/Licensing.cpp biobambam/Split12.cpp \
	biobambam/Strip12.cpp biobambam/ClipReinsert.cpp biobambam/zzToName.cpp
bam12auxmerge_LDADD = ${LIBMAUS2LIBS}
bam12auxmerge_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bam12auxmerge_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bam12split_SOURCES = programs/bam12split.cpp biobambam/Licensing.cpp biobambam/Split12.cpp
bam12split_LDADD = ${LIBMAUS2LIBS}
bam12split_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bam12split_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bam12strip_SOURCES = programs/bam12strip.cpp biobambam/Licensing.cpp biobambam/Strip12.cpp
bam12strip_LDADD = ${LIBMAUS2LIBS}
bam12strip_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bam12strip_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamreset_SOURCES = programs/bamreset.cpp biobambam/Licensing.cpp biobambam/ResetAlignment.cpp
bamreset_LDADD = ${LIBMAUS2LIBS}
bamreset_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamreset_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamrank_SOURCES = programs/bamrank.cpp biobambam/Licensing.cpp biobambam/AttachRank.cpp
bamrank_LDADD = ${LIBMAUS2LIBS}
bamrank_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamrank_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamclipreinsert_SOURCES = programs/bamclipreinsert.cpp biobambam/Licensing.cpp \
	biobambam/ClipReinsert.cpp
bamclipreinsert_LDADD = ${LIBMAUS2LIBS}
bamclipreinsert_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamclipreinsert_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamzztoname_SOURCES = programs/bamzztoname.cpp biobambam/Licensing.cpp biobambam/zzToName.cpp
bamzztoname_LDADD = ${LIBMAUS2LIBS}
bamzztoname_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamzztoname_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

kmerprob_SOURCES = programs/kmerprob.cpp biobambam/Licensing.cpp biobambam/KmerPoisson.cpp
kmerprob_LDADD = ${LIBMAUS2LIBS}
kmerprob_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
kmerprob_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamindex_SOURCES = programs/bamindex.cpp biobambam/Licensing.cpp
bamindex_LDADD = ${LIBMAUS2LIBS}
bamindex_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamindex_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilterflags_SOURCES = programs/bamfilterflags.cpp biobambam/Licensing.cpp
bamfilterflags_LDADD = ${LIBMAUS2LIBS}
bamfilterflags_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilterflags_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilterrg_SOURCES = programs/bamfilterrg.cpp biobambam/Licensing.cpp
bamfilterrg_LDADD = ${LIBMAUS2LIBS}
bamfilterrg_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilterrg_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamdownsamplerandom_SOURCES = programs/bamdownsamplerandom.cpp biobambam/Licensing.cpp
bamdownsamplerandom_LDADD = ${LIBMAUS2LIBS}
bamdownsamplerandom_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamdownsamplerandom_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamcat_SOURCES = programs/bamcat.cpp biobambam/Licensing.cpp
bamcat_LDADD = ${LIBMAUS2LIBS}
bamcat_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamcat_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bammerge_SOURCES = programs/bammerge.cpp biobambam/Licensing.cpp
bammerge_LDADD = ${LIBMAUS2LIBS}
bammerge_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bammerge_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamsplit_SOURCES = programs/bamsplit.cpp biobambam/Licensing.cpp
bamsplit_LDADD = ${LIBMAUS2LIBS}
bamsplit_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamsplit_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamsplitdiv_SOURCES = programs/bamsplitdiv.cpp biobambam/Licensing.cpp
bamsplitdiv_LDADD = ${LIBMAUS2LIBS}
bamsplitdiv_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamsplitdiv_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamchecksort_SOURCES = programs/bamchecksort.cpp biobambam/Licensing.cpp
bamchecksort_LDADD = ${LIBMAUS2LIBS}
bamchecksort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamchecksort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamrefdepth_SOURCES = programs/bamrefdepth.cpp biobambam/Licensing.cpp
bamrefdepth_LDADD = ${LIBMAUS2LIBS}
bamrefdepth_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamrefdepth_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamrefdepthpeaks_SOURCES = programs/bamrefdepthpeaks.cpp biobambam/Licensing.cpp
bamrefdepthpeaks_LDADD = ${LIBMAUS2LIBS}
bamrefdepthpeaks_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamrefdepthpeaks_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fastqtobam_SOURCES = programs/fastqtobam.cpp biobambam/Licensing.cpp
fastqtobam_LDADD = ${LIBMAUS2LIBS}
fastqtobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
fastqtobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfixmateinformation_SOURCES = programs/bamfixmateinformation.cpp biobambam/Licensing.cpp
bamfixmateinformation_LDADD = ${LIBMAUS2LIBS}
bamfixmateinformation_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfixmateinformation_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamseqchksum_SOURCES = programs/bamseqchksum.cpp biobambam/Licensing.cpp
bamseqchksum_LDADD = ${LIBMAUS2LIBS} @GMPLIBS@
bamseqchksum_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} @GMPLDFLAGS@ ${AM_LDFLAGS}
bamseqchksum_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} @GMPCPPFLAGS@

normalisefasta_SOURCES = programs/normalisefasta.cpp biobambam/Licensing.cpp
normalisefasta_LDADD = ${LIBMAUS2LIBS}
normalisefasta_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
normalisefasta_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilterheader_SOURCES = programs/bamfilterheader.cpp biobambam/Licensing.cpp
bamfilterheader_LDADD = ${LIBMAUS2LIBS}
bamfilterheader_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilterheader_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilterheader2_SOURCES = programs/bamfilterheader2.cpp biobambam/Licensing.cpp
bamfilterheader2_LDADD = ${LIBMAUS2LIBS}
bamfilterheader2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilterheader2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bammdnm_SOURCES = programs/bammdnm.cpp biobambam/Licensing.cpp
bammdnm_LDADD = ${LIBMAUS2LIBS}
bammdnm_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bammdnm_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fastabgzfextract_SOURCES = programs/fastabgzfextract.cpp biobambam/Licensing.cpp
fastabgzfextract_LDADD = ${LIBMAUS2LIBS}
fastabgzfextract_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
fastabgzfextract_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilternames_SOURCES = programs/bamfilternames.cpp biobambam/Licensing.cpp
bamfilternames_LDADD = ${LIBMAUS2LIBS}
bamfilternames_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfilternames_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamparsort_SOURCES = programs/bamparsort.cpp biobambam/Licensing.cpp
bamparsort_LDADD = ${LIBMAUS2LIBS}
bamparsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamparsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bammapdist_SOURCES = programs/bammapdist.cpp biobambam/Licensing.cpp
bammapdist_LDADD = ${LIBMAUS2LIBS}
bammapdist_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bammapdist_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamheap_SOURCES = programs/bamheap.cpp biobambam/Licensing.cpp
bamheap_LDADD = ${LIBMAUS2LIBS}
bamheap_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamheap_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfrontback_SOURCES = programs/bamfrontback.cpp biobambam/Licensing.cpp biobambam/ResetAlignment.cpp
bamfrontback_LDADD = ${LIBMAUS2LIBS}
bamfrontback_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamfrontback_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamclipextract_SOURCES = programs/bamclipextract.cpp biobambam/Licensing.cpp biobambam/ResetAlignment.cpp
bamclipextract_LDADD = ${LIBMAUS2LIBS}
bamclipextract_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamclipextract_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

blastnxmltobam_SOURCES = programs/blastnxmltobam.cpp biobambam/Licensing.cpp
blastnxmltobam_LDADD = ${LIBMAUS2LIBS} @xerces_c_LIBS@
blastnxmltobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
blastnxmltobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} @xerces_c_CFLAGS@

bamalignmentoffsets_SOURCES = programs/bamalignmentoffsets.cpp biobambam/Licensing.cpp
bamalignmentoffsets_LDADD = ${LIBMAUS2LIBS}
bamalignmentoffsets_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamalignmentoffsets_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamvalidate_SOURCES = programs/bamvalidate.cpp biobambam/Licensing.cpp
bamvalidate_LDADD = ${LIBMAUS2LIBS}
bamvalidate_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamvalidate_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamflagsplit_SOURCES = programs/bamflagsplit.cpp biobambam/Licensing.cpp
bamflagsplit_LDADD = ${LIBMAUS2LIBS}
bamflagsplit_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamflagsplit_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamintervalcomment_SOURCES = programs/bamintervalcomment.cpp biobambam/Licensing.cpp
bamintervalcomment_LDADD = ${LIBMAUS2LIBS}
bamintervalcomment_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamintervalcomment_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamintervalcommenthist_SOURCES = programs/bamintervalcommenthist.cpp biobambam/Licensing.cpp
bamintervalcommenthist_LDADD = ${LIBMAUS2LIBS}
bamintervalcommenthist_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamintervalcommenthist_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamrandomtag_SOURCES = programs/bamrandomtag.cpp biobambam/Licensing.cpp
bamrandomtag_LDADD = ${LIBMAUS2LIBS}
bamrandomtag_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamrandomtag_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamheap2_SOURCES = programs/bamheap2.cpp biobambam/Licensing.cpp
bamheap2_LDADD = ${LIBMAUS2LIBS}
bamheap2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamheap2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamalignfrac_SOURCES = programs/bamalignfrac.cpp biobambam/Licensing.cpp
bamalignfrac_LDADD = ${LIBMAUS2LIBS}
bamalignfrac_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamalignfrac_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamlastfilter_SOURCES = programs/bamlastfilter.cpp biobambam/Licensing.cpp
bamlastfilter_LDADD = ${LIBMAUS2LIBS}
bamlastfilter_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamlastfilter_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamsormadup_SOURCES = programs/bamsormadup.cpp biobambam/Licensing.cpp
bamsormadup_LDADD = ${LIBMAUS2LIBS}
bamsormadup_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamsormadup_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamexplode_SOURCES = programs/bamexplode.cpp biobambam/Licensing.cpp
bamexplode_LDADD = ${LIBMAUS2LIBS}
bamexplode_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamexplode_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fastqtobampar_SOURCES = programs/fastqtobampar.cpp biobambam/Licensing.cpp
fastqtobampar_LDADD = ${LIBMAUS2LIBS}
fastqtobampar_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
fastqtobampar_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}
