ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS=biobambam/BamBamConfig.hpp biobambam/Licensing.hpp \
	biobambam/ClipAdapters.hpp biobambam/AttachRank.hpp biobambam/ResetAlignment.hpp \
	biobambam/Split12.hpp biobambam/Strip12.hpp

man_MANS = programs/bamtofastq.1 programs/bamsort.1 programs/bammarkduplicates.1 programs/bamcollate.1 \
	programs/bammaskflags.1 programs/bamrecompress.1 programs/bamadapterfind.1 \
	programs/bamfilteraux.1 programs/bamauxsort.1 programs/bamadapterclip.1 \
	programs/bamcollate2.1 programs/bam12auxmerge.1 programs/bam12split.1 \
	programs/bam12strip.1 programs/bamreset.1 programs/bamrank.1 \
	programs/bamclipreinsert.1 programs/bamzztoname.1

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
	bamzztoname
	
noinst_PROGRAMS = bamfilter bamfixmatecoordinates bamfixmatecoordinatesnamesorted bamcheckalignments bamtoname \
	bamdisthist

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
	biobambam/ClipAdapters.cpp
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
	biobambam/Strip12.cpp
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

bamclipreinsert_SOURCES = programs/bamclipreinsert.cpp biobambam/Licensing.cpp
bamclipreinsert_LDADD = ${LIBMAUSLIBS}
bamclipreinsert_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamclipreinsert_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamzztoname_SOURCES = programs/bamzztoname.cpp biobambam/Licensing.cpp
bamzztoname_LDADD = ${LIBMAUSLIBS}
bamzztoname_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamzztoname_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
