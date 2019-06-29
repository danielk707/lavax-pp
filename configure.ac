#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([lavax], [0.1], [danielk5@kth.se])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
dnl m4_include([m4/ac_check_blitz.m4])
dnl m4_include([m4/boost.m4])
dnl m4_include([m4/ax_cxx_check_lib.m4])
dnl m4_include([m4/ax_boost_base.m4])
dnl m4_include([m4/ax_boost_system.m4])
dnl m4_include([m4/ax_boost_filesystem.m4])

dnl AC_LIB_LINKFLAGS([boost_system])
dnl AC_LIB_LINKFLAGS([boost_filesystem])
dnl AC_LIB_LINKFLAGS([blitzz])


dnl AC_CHECK_BLITZ([],[AC_MSG_ERROR_BLITZ])



AM_INIT_AUTOMAKE([foreign])
dnl AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
dnl AC_LANG_PUSH([C++])
#AC_CHECK_HEADERS([blitz/array.h], [], [AC_MSG_ERROR([Unable to find required blitz headers])])
# AC_CHECK_HEADERS([boost/units/base_units/metric/angstrom.hpp], [],
#                 [AC_MSG_ERROR([Unable to find required boost headers])])
dnl AX_CXX_CHECK_LIB(boost_filesystem, [boost::filesystem::path::begin() const])
dnl AX_BOOST_SYSTEMx
                 dnl [AC_MSG_ERROR([Unable to find required boost shared library])])
dnl AC_LANG_POP([C++])

AC_SEARCH_LIBS([pthread_create], [pthread])
dnl AC_SEARCH_LIBS([chmod], [boost_filesystem])

${CXXFLAGS=""}
AC_PROG_CXX

dnl AC_LANG_PUSH([C++])
dnl AX_BOOST_BASE([1.48], [], [])
dnl AX_BOOST_SYSTEM
dnl AX_BOOST_FILESYSTEM
dnl AC_LANG_POP([C++])

dnl AC_PROG_LIBTOOL
dnl BOOST_REQUIRE([1.48])
dnl BOOST_ARRAY
dnl BOOST_DATE_TIME
dnl BOOST_SYSTEM
dnl BOOST_THREAD
dnl BOOST_FILESYSTEM([p])

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])

AC_OUTPUT
