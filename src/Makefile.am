ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS=\
	biobambam/bamToFastQ.hpp \
	biobambam/BamFile.hpp \
	biobambam/PutFastQBase.hpp \
	biobambam/OutputBuffer.hpp \
	biobambam/BamBamConfig.hpp \
	biobambam/Licensing.hpp

bin_PROGRAMS = bamtofastq bammarkduplicates bamsort bamcollate bammaskflags @SCRAMTOFASTQ@
noinst_PROGRAMS = bamfilter bamfixmatecoordinates bamfixmatecoordinatesnamesorted bamcheckalignments
EXTRA_PROGRAMS = scramtofastq

bamcollate_SOURCES = programs/bamcollate.cpp biobambam/Licensing.cpp
bamcollate_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcollate_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamcollate_LDADD = ${LIBMAUSLIBS}

bamfilter_SOURCES = programs/bamfilter.cpp biobambam/Licensing.cpp
bamfilter_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilter_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamfilter_LDADD = ${LIBMAUSLIBS}

bamfixmatecoordinates_SOURCES = programs/bamfixmatecoordinates.cpp biobambam/Licensing.cpp
bamfixmatecoordinates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfixmatecoordinates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamfixmatecoordinates_LDADD = ${LIBMAUSLIBS}

bamfixmatecoordinatesnamesorted_SOURCES = programs/bamfixmatecoordinatesnamesorted.cpp biobambam/Licensing.cpp
bamfixmatecoordinatesnamesorted_LDADD = ${LIBMAUSLIBS}
bamfixmatecoordinatesnamesorted_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${SAMTOOLSLDFLAGS}  ${AM_LDFLAGS}
bamfixmatecoordinatesnamesorted_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${SAMTOOLSCPPFLAGS}

bammarkduplicates_SOURCES = programs/bammarkduplicates.cpp biobambam/Licensing.cpp
bammarkduplicates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammarkduplicates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bammarkduplicates_LDADD = ${LIBMAUSLIBS}

bammaskflags_SOURCES = programs/bammaskflags.cpp biobambam/Licensing.cpp
bammaskflags_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammaskflags_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bammaskflags_LDADD = ${LIBMAUSLIBS}

bamsort_SOURCES = programs/bamsort.cpp biobambam/Licensing.cpp
bamsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamsort_LDADD = ${LIBMAUSLIBS}

bamtofastq_SOURCES = programs/bamtofastq.cpp biobambam/bamToFastQ.cpp biobambam/Licensing.cpp
bamtofastq_LDADD = ${SAMTOOLSLIBS} ${LIBMAUSLIBS}
bamtofastq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${SAMTOOLSLDFLAGS}
bamtofastq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${SAMTOOLSCPPFLAGS}

bamcheckalignments_SOURCES = programs/bamcheckalignments.cpp biobambam/Licensing.cpp
bamcheckalignments_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcheckalignments_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamcheckalignments_LDADD = ${LIBMAUSLIBS}

scramtofastq_SOURCES = programs/bamtofastq.cpp biobambam/bamToFastQ.cpp biobambam/Licensing.cpp
scramtofastq_LDADD = ${LIBMAUSLIBS}
scramtofastq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS}
scramtofastq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} -DBAMTOFASTQ_USE_LIBMAUS_IO_LIB
