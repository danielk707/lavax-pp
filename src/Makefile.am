# AM_CPPFLAGS = -std=c++14 -I$(top_srcdir)/include -I/usr/local/include/eigen3 -O0 -lblitz -lm
# bin_PROGRAMS = qc mkbcc PKAvel mklammps_data
# bin_PROGRAMS = lavax unit_test #test1 test2 test3
# lavax_SOURCES = main_lavax.cpp lavax.cpp lavax.hpp
# lavax_CXXFLAGS = -lpthread

pkgdata_DATA = lavax.conf predictor.in

# AM_CPPFLAGS = -std=c++14 -lblitz -lm  -fpermissive -lpthread -DDATADIR=\"$(pkgdatadir)\"
# AM_CXXFLAGS = -std=c++14 -fpermissive -lpthread -I/cfs/klemming/scratch/d/danielk5/local/include -L/cfs/klemming/scratch/d/danielk5/local/lib -I$BOOST_HOME/include -L$BOOST_HOME/lib -DDATADIR=\"$(pkgdatadir)\"
AM_CXXFLAGS = -std=c++14 -fpermissive -lpthread -DDATADIR=\"$(pkgdatadir)\"
# AM_CPPFLAGS = -std=c++1y -lblitz -lm  -fpermissive -lpthread -nostdlib -nostdinc -DDATADIR=\"/home/danielk5/share/lavax\"
#AM_CPPFLAGS = -std=c++1y -fpermissive -DDATADIR=\"/home/danielk5/local/share/lavax\"
#AM_LDFLAGS = -lboost_system -lboost_filesystem
#AM_LDFLAGS=-L/home/mcg04/local/lib
#AM_LDFLAGS=-L/usr/local/lib
#AM_LDFLAGS=-L/home/danielk5/local/lib

bin_PROGRAMS = lavax PKAvel
# unit_test_SOURCES = unit_test.cpp lavax.cpp init_check.cpp
#unit_test_LDFLAGS = -lboost_system -lboost_filesystem
# unit_test_LDADD = -lboost_system -lboost_filesystem

lavax_SOURCES = main_lavax.cpp lavax.cpp init_check.cpp
lavax_LDADD = -lboost_system -lboost_filesystem -lblitz -lm #-lpthread

PKAvel_SOURCES = main_PKAvel.cpp
# lavax_LDFLAGS = -lboost_system -lboost_filesystem

# mklammps_data_SOURCES = mklammps_data.cpp qc.cpp qc.hpp

# # -------------------- JAVA --------------------
# javadir = $(datadir)/java
# java_JAVA = lavax.java Vec3.java Particle.java
