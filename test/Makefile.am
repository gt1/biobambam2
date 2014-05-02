TESTS=\
	testfastqbamloop.sh \
	testshortsortcoordinate.sh \
	testshortsortqueryname.sh \
	testshortsort.sh \
	testdupsingle.sh
TEST_ENVIRONMENT= 
LOG_COMPILER=/bin/bash
EXTRA_DIST= sorttestshort.sh testfastqbamloop.sh testshortsortcoordinate.sh testshortsortqueryname.sh testshortsort.sh testdupsingle.sh # 

check_PROGRAMS=bamcmp

bamcmp_SOURCES = bamcmp.cpp
bamcmp_LDADD = ${LIBMAUSLIBS}
bamcmp_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS} ${LIBMAUSLDFLAGS} ${AM_LDFLAGS}
bamcmp_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUSCPPFLAGS}
