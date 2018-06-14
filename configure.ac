#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([lavax], [0.1], [danielk5@kth.se])
AM_INIT_AUTOMAKE([foreign])
dnl AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

AC_SEARCH_LIBS([pthread_create], [pthread])

${CXXFLAGS=""}
AC_PROG_CXX

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
