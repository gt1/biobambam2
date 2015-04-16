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
bamcmp_LDADD = ${LIBMAUS2LIBS}
bamcmp_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamcmp_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamtosam_SOURCES = bamtosam.cpp
bamtosam_LDADD = ${LIBMAUS2LIBS}
bamtosam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamtosam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}
