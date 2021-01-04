# ==========================================================================
# Copyright(c)'2007 by PRiSM Laboratory Opale Team
# see the copyright file for more information
# Author: B. Le Cun
# Makefile.am for example/Qap Bob++ library
#
# ==========================================================================

if USE_ATHA
  PROG_ATHA=pfsatha
endif

if USE_MPI
  PROG_MPI=pfsmpi
endif

EXTRA_DIST = 


noinst_PROGRAMS= pfs pfsthr ${PROG_ATHA} ${PROG_MPI}

CLEANFILES=${noinst_PROGRAMS} *.o

LDADD=-L${prefix}/lib -lbobpp-${GENERIC_API_VERSION}
AM_CPPFLAGS=-I${prefix}/include/${GENERIC_LIBRARY_NAME}-${GENERIC_API_VERSION} -g

pfs_SOURCES = pfs.cpp pfs.h lb.h
pfs_CPPFLAGS = $(AM_CPPFLAGS) ${INCGLOP}
pfs_LDFLAGS = 

pfsthr_SOURCES = pfs.cpp pfs.h lb.h
pfsthr_CPPFLAGS = $(AM_CPPFLAGS) ${INCGLOP} -DThreaded
pfsthr_LDFLAGS = 


if USE_ATHA
  pfsatha_SOURCES=pfs.cpp pfs.h lb.h
  pfsatha_LDFLAGS = ${LIBATHA} -lbobpp_atha
  pfsatha_CPPFLAGS = $(AM_CPPFLAGS) ${INCATHA} -DAtha
endif

if USE_MPI
MPICXXCOMPILE = $(MPICXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
MPICXXLINK = $(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
       --mode=link $(MPICXX) $(AM_CXXFLAGS) $(CXXFLAGS) \
       -lbobpp_mpi -o $@
  pfsmpi_SOURCES = pfs.cpp pfs.h lb.h
  pfsmpi_CPPFLAGS = $(AM_CPPFLAGS) ${INCGLOP} -DMPxMPI
  pfsmpi_LINK = $(MPICXXLINK) ${LIBGLOP}
pfsmpi-pfs.o: pfs.cpp
	$(MPICXXCOMPILE) ${INCGLOP} -DMPxMPI -c -o $@ $<
endif

