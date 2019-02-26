AC_INIT([mivdesy-tests], [0.3], [zdenek.matej@maxiv.lu.se])
PKG_PROG_PKG_CONFIG
AC_CONFIG_MACRO_DIR([m4])

AM_EXTRA_RECURSIVE_TARGETS([modules])

AC_LANG([Fortran])

AM_INIT_AUTOMAKE([subdir-objects])

AC_ARG_VAR([IORBSIZE], AS_HELP_STRING([Block size for ior -b, ], [use e.g. IORBSIZE=8M]))
AC_ARG_VAR([H5PW_GN],  AS_HELP_STRING([Image size in pixels for pwrite, ], [use e.g. H5PW_GN=8000]))
AC_ARG_VAR([H5PW_CHN], AS_HELP_STRING([Chunk size in pixels for pwrite, ], [use e.g. H5PW_CHN=1000]))
AC_ARG_VAR([MPICMD],   AS_HELP_STRING([Command to execute mpi progrmas, ], [use e.g. MPICMD=/usr/bin/mpirun]))
AC_ARG_VAR([TESTFSPATH], AS_HELP_STRING([Path to testing filesystem, ], [use e.g. TESTFSPATH=/gpfs/gpfs8m/]))
AC_ARG_VAR([HDF5_USE_SHLIB], AS_HELP_STRING([Use hdf5 shared libraries, ], [use e.g. HDF5_USE_SHLIB=yes]))

# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],[compile with MPI (parallelization) support. If none is found,MPI is not used. Default: auto])],,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])

AC_PROG_CC([h5pcc])
AC_PROG_FC([h5pfc])

# default values
if test x$IORBSIZE = x; then
	IORBSIZE=8M
fi

if test x$H5PW_GN = x; then
        H5PW_GN=8000
fi

if test x$H5PW_CHN = x; then
        H5PW_CHN=1000
fi

if test x$MPICMD = x; then
        MPICMD=$(which mpirun)
fi

if test x$TESTFSPATH = x; then
        TESTFSPATH=$(pwd)
fi

AC_SUBST(IORBSIZE)
AC_SUBST(H5PW_GN)
AC_SUBST(H5PW_CHN)
AC_SUBST(MPICMD)
AC_SUBST(TESTFSPATH)
FC=$(which h5pfc)
AC_SUBST(FC)
AC_SUBST(PATH)
AC_SUBST(HDF5_USE_SHLIB)


AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

echo
echo "Summary of configuration"

echo "  CC         : $CC"
echo "  FC         : $FC"
echo "  mpicc      : $(which mpicc)"
echo "  h5pfc      : $(which h5pfc)"
echo "  mpicmd     : $MPICMD"

echo "  IORBSIZE   : $IORBSIZE"
echo "  H5PW_GN    : $H5PW_GN"
echo "  H5PW_CHN   : $H5PW_CHN"
echo "  LIBS       : $LIBS"
echo "  TESTFSPATH : $TESTFSPATH"

mkdir -p conf
echo $IORBSIZE > conf/iorbsize
echo $H5PW_GN  > conf/h5pw_gn
echo $H5PW_CHN > conf/h5pw_chn
echo $TESTFSPATH > conf/testfspath
echo $MPICMD > conf/mpicmd

