AC_INIT(biobambam,0.0.191,[gt1@sanger.ac.uk],[biobambam],[http://www.sanger.ac.uk])
AC_CANONICAL_SYSTEM
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE([std-options subdir-objects parallel-tests])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_SYSTEM
AC_MSG_CHECKING([whether we are cross compiling])
ISNATIVE=
if test \( "${ac_cv_build}" = "${ac_cv_host}" \) -a \( "${ac_cv_build}" = "${ac_cv_target}" \) ; then
	ISNATIVE=1
	AC_MSG_RESULT([no])
else
	ISNATIVE=0
	AC_MSG_RESULT([yes])
fi

AC_ARG_ENABLE(optimization,
        AS_HELP_STRING([--enable-optimization],[use compiler optimization (default yes)]),
        [optimization=${enableval}],[optimization=yes])
AC_ARG_ENABLE(debug,
        AS_HELP_STRING([--enable-debug],[use compiler debug flags (default no)]),
        [debug=${enableval}],[debug=no])
AC_ARG_ENABLE(profile,
        AS_HELP_STRING([--enable-profile],[use compiler profiling flags (default no)]),
        [profile=${enableval}],[profile=no])
AC_ARG_ENABLE(fast,
        AS_HELP_STRING([--enable-fast],[disable evaluation of assertions (default no)]),
        [fast=${enableval}],[fast=no])
AC_ARG_ENABLE(install-experimental,
        AS_HELP_STRING([--enable-install-experimental],[install experimental programs (default no)]),
        [install_experimental=${enableval}],[install_experimental=no])

if test "${install_experimental}" = "yes" ; then
	BAMPARSORTINSTEXP=bamparsort
	BAMALIGNMENTOFFSETSINSTEXP=bamalignmentoffsets
else
	BAMPARSORTNOINSTEXP=bamparsort
	BAMALIGNMENTOFFSETSNOINSTEXP=bamalignmentoffsets
fi


# check for C compiler
CFLAGS_SAVE="${CFLAGS}"
CFLAGS=""
AC_PROG_CC
CFLAGS="${CFLAGS_SAVE}"

CFLAGS_SAVE="${CFLAGS}"
CFLAGS=""
CXXFLAGS_SAVE="${CXXFLAGS}"
CXXFLAGS=""
AM_PROG_CC_C_O
CFLAGS="${CFLAGS_SAVE}"
CXXFLAGS="${CXXFLAGS_SAVE}"

# check for C++ compiler
CXXFLAGS_SAVE="${CXXFLAGS}"
CXXFLAGS=""
AC_PROG_CXX
CXXFLAGS="${CXXFLAGS_SAVE}"

if test "${debug}" = "yes" ; then
	CXXDEBUGFLAGS="-g -O0 -rdynamic"
	CDEBUGFLAGS="-g -O0 -rdynamic"
        CXXFLAGS="${CXXFLAGS} ${CXXDEBUGFLAGS}"
        CFLAGS="${CFLAGS} ${CDEBUGFLAGS}"
	openmp="no"

	AC_MSG_NOTICE([Adding debug options ${CDEBUGFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])
	AC_MSG_NOTICE([Adding debug options ${CXXDEBUGFLAGS} to C compiler ${CXX} flags, full list now ${CXXFLAGS}])
else
	if test "${profile}" = "yes" ; then
	        CXXPROFFLAGS="-g -pg -rdynamic"
	        CPROFFLAGS="-g -pg -rdynamic"

	        CXXFLAGS="${CXXPROFFLAGS} -g -pg -rdynamic"
	        CFLAGS="${CPROFFLAGS} -g -pg -rdynamic"

		AC_MSG_NOTICE([Adding profiling options ${CPROFFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])
		AC_MSG_NOTICE([Adding profiling options ${CXXPROFFLAGS} to C compiler ${CXX} flags, full list now ${CXXFLAGS}])
	else
		if test "${optimization}" = "yes" ; then
			COPTFLAGS=
			case ${CC} in
				gcc)
					COPTFLAGS="-O3 -rdynamic"
					;;
				*-mingw32msvc-gcc)
					COPTFLAGS="-O3 -rdynamic"
					;;
				cl.exe)
					COPTFLAGS="-O2 -Ob2 -Ot -Oy"
					;;
			esac
			CFLAGS="${CFLAGS} ${COPTFLAGS}"
			AC_MSG_NOTICE([Adding optimisation options ${COPTFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])

			CXXOPTFLAGS=
			case ${CXX} in
				g++)
					CXXOPTFLAGS="-O3 -rdynamic"
					;;
				*-mingw32msvc-g++)
					CXXOPTFLAGS="-O3 -rdynamic"
					;;
				cl.exe)
					CXXOPTFLAGS="-O2 -Ob2 -Ot -Oy"
					;;
			esac
			CXXFLAGS="${CXXFLAGS} ${CXXOPTFLAGS}"

			AC_MSG_NOTICE([Adding optimisation options ${CXXOPTFLAGS} to C++ compiler ${CXX} flags, full list now ${CXXFLAGS}])
		fi
	fi
fi


# check for pkg-config tool
PKG_PROG_PKG_CONFIG

if test "${fast}" = "yes" ; then
	CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi

# check for standard headers
AC_LANG_PUSH([C++])
AC_HEADER_STDC
AC_LANG_POP

if test "${CXX}" = "g++" ; then
	CPPFLAGS="${CPPFLAGS} -W -Wall"
fi

AC_ARG_WITH([libmaus],
            [AS_HELP_STRING([--with-libmaus@<:@=PATH@:>@], [path to installed libmaus library @<:@default=@:>@])],
            [with_libmaus=${withval}],
            [with_libmaus=])

if test ! -z "${with_libmaus}" ; then
	PKGCONFIGPATHSAVE="${PKG_CONFIG_PATH}"
	if test -z "${PKG_CONFIG_PATH}" ; then
		export PKG_CONFIG_PATH="${with_libmaus}/lib/pkgconfig"
	else
		export PKG_CONFIG_PATH="${with_libmaus}/lib/pkgconfig:${PKG_CONFIG_PATH}"
	fi
fi

PKG_CHECK_MODULES([libmaus],[libmaus >= 0.0.196])
PKG_CHECK_MODULES([libmausdigests],[libmausdigests >= 0.0.196])

if test ! -z "${with_libmaus}" ; then
	if test ! -z "${PKGCONFIGPATHSAVE}" ; then
		export PKG_CONFIG_PATH="${PKGCONFIGPATHSAVE}"
	fi
fi

LIBMAUSCPPFLAGS="${libmaus_CFLAGS} ${libmausdigests_CFLAGS}"
LIBMAUSLIBS="${libmaus_LIBS} ${libmausdigests_LIBS}"

CPPFLAGS_SAVE="${CPPFLAGS}"
LDFLAGS_SAVE="${LDFLAGS}"
LIBS_SAVE="${LIBS}"

CPPFLAGS="${CPPFLAGS} ${libmaus_CFLAGS}"
LIBS="${LIBS} ${libmaus_LIBS}"

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether we can compile a program using libmaus])
	AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <libmaus/util/ArgInfo.hpp>

				int main(int argc, char * argv[[]]) {
					::libmaus::util::ArgInfo const arginfo(argc,argv);
					return 0;
			}])],
			have_libmaus=yes,
			have_libmaus=no
		)
AC_LANG_POP
AC_MSG_RESULT($have_libmaus)

if test "${have_libmaus}" = "no" ; then
	AC_MSG_ERROR([Required libmaus is not available.])
fi

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether libmaus has support for io_lib])
	AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <libmaus/LibMausConfig.hpp>

				int main(int argc, char * argv[[]]) {
					#if ! defined(LIBMAUS_HAVE_IO_LIB)
					#error "installed libmaus does not support io_lib"
					#endif
					return 0;
			}])],
			have_libmaus_iolib=yes,
			have_libmaus_iolib=no
		)
AC_LANG_POP
AC_MSG_RESULT($have_libmaus_iolib)

LIBMAUSIOLIBDEFINE=
if test "${have_libmaus_iolib}" = "yes" ; then
	LIBMAUSIOLIBDEFINE="#define BIOBAMBAM_LIBMAUS_HAVE_IO_LIB"
fi

## zzzz

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether libmaus has support for iRODS])
	AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <libmaus/LibMausConfig.hpp>

				int main(int argc, char * argv[[]]) {
					#if ! defined(LIBMAUS_HAVE_IRODS)
					#error "installed libmaus does not support iRODS"
					#endif
					return 0;
			}])],
			have_libmaus_irods=yes,
			have_libmaus_irods=no
		)
AC_LANG_POP
AC_MSG_RESULT($have_libmaus_irods)

LIBMAUSIRODSDEFINE=
LIBMAUSIRODSCPPFLAGS=
LIBMAUSIRODSLIBS=
if test "${have_libmaus_irods}" = "yes" ; then
	if test ! -z "${with_libmaus}" ; then
		PKGCONFIGPATHSAVE="${PKG_CONFIG_PATH}"
		if test -z "${PKG_CONFIG_PATH}" ; then
			export PKG_CONFIG_PATH="${with_libmaus}/lib/pkgconfig"
		else
			export PKG_CONFIG_PATH="${with_libmaus}/lib/pkgconfig:${PKG_CONFIG_PATH}"
		fi
	fi

	PKG_CHECK_MODULES([libmausirods],[libmausirods >= 0.0.196])

	LIBMAUSIRODSCPPFLAGS="${libmausirods_CFLAGS}"
	LIBMAUSIRODSLIBS="${libmausirods_LIBS}"
	
	LIBMAUSCPPFLAGS="${LIBMAUSCPPFLAGS} ${LIBMAUSIRODSCPPFLAGS}"
	LIBMAUSLIBS="${LIBMAUSLIBS} ${LIBMAUSIRODSLIBS}"
	CPPFLAGS="${CPPFLAGS} ${LIBMAUSIRODSCPPFLAGS}"
	LIBS="${LIBS} ${LIBMAUSIRODSLIBS}"

	if test ! -z "${with_libmaus}" ; then
		if test ! -z "${PKGCONFIGPATHSAVE}" ; then
			export PKG_CONFIG_PATH="${PKGCONFIGPATHSAVE}"
		fi
	fi
	
	LIBMAUSIRODSDEFINE="#define BIOBAMBAM_LIBMAUS_HAVE_IRODS"
fi

CPPFLAGS="${CPPFLAGS_SAVE}"
LDFLAGS="${LDFLAGS_SAVE}"
LIBS="${LIBS_SAVE}"

AC_ARG_WITH([liftingwavelettransform],
            [AS_HELP_STRING([--with-liftingwavelettransform@<:@=PATH@:>@], [path to installed LiftingWaveletTransform library @<:@default=no@:>@])],
            [with_liftingwavelettransform=${withval}],
            [with_liftingwavelettransform=no])

if test ! -z "${with_liftingwavelettransform}" ; then
	PKGCONFIGPATHSAVE="${PKG_CONFIG_PATH}"
	if test -z "${PKG_CONFIG_PATH}" ; then
		export PKG_CONFIG_PATH="${with_liftingwavelettransform}/lib/pkgconfig"
	else
		export PKG_CONFIG_PATH="${with_liftingwavelettransform}/lib/pkgconfig:${PKG_CONFIG_PATH}"
	fi
fi

if test "${with_liftingwavelettransform}" != "no" ; then
	PKG_CHECK_MODULES([liftingwavelettransform],[liftingwavelettransform >= 0.0.0])
	LIFTINGWAVELETTRANSFORMDEFINE="#define BIOBAMBAM_HAVE_LIFTINGWAVELETTRANSFORM"
	BAMREFDEPTHPEAKS="bamrefdepthpeaks"

	LIFTINGWAVELETTRANSFORMCPPFLAGS="${liftingwavelettransform_CFLAGS}"
	LIFTINGWAVELETTRANSFORMLIBS="${liftingwavelettransform_LIBS}"

	CPPFLAGS="${CPPFLAGS} ${LIFTINGWAVELETTRANSFORMCPPFLAGS}"
	LIBS="${LIBS} ${LIFTINGWAVELETTRANSFORMLIBS}"
fi

AC_PATH_PROG(BASHPROG,bash,no)
AC_PATH_PROG(BASE64PROG,base64,no)

TESTDIR=
if test \( "${ISNATIVE}" -eq 1 \) -a \( "${BASHPROG}" != "no" \) -a \( "${BASE64PROG}" != "no" \) ; then
	TESTDIR=test
fi

AC_ARG_WITH([xerces-c],
            [AS_HELP_STRING([--with-xerces-c@<:@=PATH@:>@], [path to installed xerces-c library @<:@default=no@:>@])],
            [with_xerces_c=${withval}],
            [with_xerces_c=no])

BIOBAMBAM_HAVE_XERCES_C=
xerces_c_CFLAGS=
xerces_c_LIBS=
BLASTNXMLTOBAM=
if test "${with_xerces_c}" != "no" ; then
	if test ! -z "${with_xerces_c}" ; then
		PKGCONFIGPATHSAVE="${PKG_CONFIG_PATH}"
		if test -z "${PKG_CONFIG_PATH}" ; then
			export PKG_CONFIG_PATH="${with_xerces_c}/lib/pkgconfig"
		else
			export PKG_CONFIG_PATH="${with_xerces_c}/lib/pkgconfig:${PKG_CONFIG_PATH}"
		fi
	fi

	PKG_CHECK_MODULES([xerces_c],[xerces-c >= 3.1.1])

	if test ! -z "${with_xerces_c}" ; then
		if test ! -z "${PKGCONFIGPATHSAVE}" ; then
			export PKG_CONFIG_PATH="${PKGCONFIGPATHSAVE}"
		fi
	fi
	
	BIOBAMBAM_HAVE_XERCES_C="#define BIOBAMBAM_HAVE_XERCES_C"
	BLASTNXMLTOBAM="blastnxmltobam"
fi

AC_ARG_WITH([gmp],
            [AS_HELP_STRING([--with-gmp@<:@=PATH@:>@], [path to installed gmp library @<:@default=no@:>@])],
            [with_gmp=${withval}],
            [with_gmp=no])

BIOBAMBAM_HAVE_GMP=
gmp_CFLAGS=
gmp_LIBS=
if test "${with_gmp}" != "no" ; then

	
	AC_MSG_CHECKING([whether we can compile a program using the gmp library])
	LDFLAGS_SAVE=${LDFLAGS}
	CPPFLAGS_SAVE=${CPPFLAGS}
	LIBS_SAVE=${LIBS}
	
	if test \( ! -z "${with_gmp}" \) -a \( "${with_gmp}" != "yes" \) ; then
		LDFLAGS="-L${with_gmp}/lib"
		CPPFLAGS="-I${with_gmp}/include"
	fi

	LIBS="-lgmp ${LIBS}"

	AC_LANG_PUSH([C++])
	AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <gmp.h>

int main(int argc, char * argv[[]]) {
	return 0;
}])],
			have_gmp=yes,
			have_gmp=no
		)
	AC_LANG_POP
	AC_MSG_RESULT($have_gmp)

	LDFLAGS=${LDFLAGS_SAVE}
	CPPFLAGS=${CPPFLAGS_SAVE}
	LIBS=${LIBS_SAVE}

	if test "${have_gmp}" = "yes" ; then
		if test "${with_gmp}" != "yes" ; then
			GMPLDFLAGS="-L${with_gmp}/lib"
			GMPCPPFLAGS="-I${with_gmp}/include"
		else
			GMPLDFLAGS=
			GMPCPPFLAGS=
		fi
		BIOBAMBAM_HAVE_GMP="#define BIOBAMBAM_HAVE_GMP"
		GMPLIBS="-lgmp"
	else
		GMPLDFLAGS=
		GMPCPPFLAGS=
		GMPLIBS=
	fi
fi


if test "${install_experimental}" = "yes" ; then
	BLASTXMLTOBAMINSTEXP=${BLASTNXMLTOBAM}
else
	BLASTXMLTOBAMNOINSTEXP=${BLASTNXMLTOBAM}
fi

AC_MSG_NOTICE([Using flags ${CFLAGS} for C compiler ${CC}])
AC_MSG_NOTICE([Using flags ${CXXFLAGS} for C++ compiler ${CXX}])

PACKAGE_NAME=${PACKAGE}
PACKAGE_VERSION=${VERSION}

AC_SUBST([PACKAGE_NAME])
AC_SUBST([PACKAGE_VERSION])
AC_SUBST([BAMBAM_CFLAGS])
AC_SUBST([BAMBAM_LDFLAGS])
AC_SUBST([BAMBAM_LIBS])
# libmaus
AC_SUBST([LIBMAUSCPPFLAGS])
AC_SUBST([LIBMAUSLIBS])
# 
AC_SUBST([LIBMAUSIOLIBDEFINE])
AC_SUBST([LIBMAUSIRODSDEFINE])
#
AC_SUBST([LIFTINGWAVELETTRANSFORMDEFINE])
AC_SUBST([BAMREFDEPTHPEAKS])
#
AC_SUBST([BAMPARSORTINSTEXP])
AC_SUBST([BAMPARSORTNOINSTEXP])
AC_SUBST([BAMALIGNMENTOFFSETSINSTEXP])
AC_SUBST([BAMALIGNMENTOFFSETSNOINSTEXP])
AC_SUBST([BLASTXMLTOBAMINSTEXP])
AC_SUBST([BLASTXMLTOBAMNOINSTEXP])
AC_SUBST([TESTDIR])
#
AC_SUBST([xerces_c_CFLAGS])
AC_SUBST([xerces_c_LIBS])
AC_SUBST([BIOBAMBAM_HAVE_XERCES_C])
AC_SUBST([BLASTNXMLTOBAM])
AC_SUBST([BIOBAMBAM_HAVE_GMP])
AC_SUBST([GMPLDFLAGS])
AC_SUBST([GMPCPPFLAGS])
AC_SUBST([GMPLIBS])
# 
AC_OUTPUT(Makefile src/Makefile test/Makefile src/biobambam/BamBamConfig.hpp)
