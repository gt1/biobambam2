TESTS=\
	testfastqbamloop.sh \
	testshortsortcoordinate.sh \
	testshortsortqueryname.sh \
	testshortsort.sh \
	testdupsingle.sh \
	testdupsinglemarkedsortedqreset.sh
TEST_ENVIRONMENT= 
LOG_COMPILER=/bin/bash
EXTRA_DIST= dupsingle.sh dupsinglemarked.sh sorttestshort.sh dupsinglemarkedsortedqreset.sh \
	testfastqbamloop.sh testshortsortcoordinate.sh testshortsortqueryname.sh testshortsort.sh testdupsingle.sh \
	testdupsinglemarkedsortedqreset.sh base64decode.sh #

check_PROGRAMS=bamcmp bamtosam

bamcmp_SOURCES = bamcmp.cpp
bamcmp_LDADD = ${LIBMAUSLIBS}
bamcmp_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcmp_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}

bamtosam_SOURCES = bamtosam.cpp
bamtosam_LDADD = ${LIBMAUSLIBS}
bamtosam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamtosam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
