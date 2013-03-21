ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS=\
	bambam/bamToFastQ.hpp \
	bambam/BamFile.hpp \
	bambam/PutFastQBase.hpp \
	bambam/OutputBuffer.hpp \
	bambam/BamBamConfig.hpp \
	bambam/Licensing.hpp

bin_PROGRAMS = bamtofastq bammarkduplicates bamsort bamcollate bammaskflags
noinst_PROGRAMS = bamfilter bamfixmatecoordinates bamfixmatecoordinatesnamesorted bamcheckalignments

bamcollate_SOURCES = programs/bamcollate.cpp bambam/Licensing.cpp
bamcollate_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcollate_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamcollate_LDADD = ${LIBMAUSLIBS}

bamfilter_SOURCES = programs/bamfilter.cpp bambam/Licensing.cpp
bamfilter_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfilter_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamfilter_LDADD = ${LIBMAUSLIBS}

bamfixmatecoordinates_SOURCES = programs/bamfixmatecoordinates.cpp bambam/Licensing.cpp
bamfixmatecoordinates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamfixmatecoordinates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamfixmatecoordinates_LDADD = ${LIBMAUSLIBS}

bamfixmatecoordinatesnamesorted_SOURCES = programs/bamfixmatecoordinatesnamesorted.cpp bambam/Licensing.cpp
bamfixmatecoordinatesnamesorted_LDADD = ${LIBMAUSLIBS}
bamfixmatecoordinatesnamesorted_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${SAMTOOLSLDFLAGS}  ${AM_LDFLAGS}
bamfixmatecoordinatesnamesorted_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${SAMTOOLSCPPFLAGS}

bammarkduplicates_SOURCES = programs/bammarkduplicates.cpp bambam/Licensing.cpp
bammarkduplicates_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammarkduplicates_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bammarkduplicates_LDADD = ${LIBMAUSLIBS}

bammaskflags_SOURCES = programs/bammaskflags.cpp bambam/Licensing.cpp
bammaskflags_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bammaskflags_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bammaskflags_LDADD = ${LIBMAUSLIBS}

bamsort_SOURCES = programs/bamsort.cpp bambam/Licensing.cpp
bamsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamsort_LDADD = ${LIBMAUSLIBS}

bamtofastq_SOURCES = programs/bamtofastq.cpp bambam/bamToFastQ.cpp bambam/Licensing.cpp
bamtofastq_LDADD = ${SAMTOOLSLIBS} ${LIBMAUSLIBS}
bamtofastq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${SAMTOOLSLDFLAGS}
bamtofastq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${SAMTOOLSCPPFLAGS}

bamcheckalignments_SOURCES = programs/bamcheckalignments.cpp bambam/Licensing.cpp
bamcheckalignments_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcheckalignments_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
bamcheckalignments_LDADD = ${LIBMAUSLIBS}