AC_INIT(test, version-1)

m4_include([m4/check_openmp.m4])
m4_include([m4/check_autotype.m4])


AC_LANG(C++)
AC_PROG_CXX

LIBS="-lm "

dnl Small macro to print out a nice, pretty section title.
define(SECTION_TITLE,
[
	echo >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
	echo '   $1' >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
])


SECTION_TITLE([Checking for modern C++ features])

AC_CHECK_HEADERS(initializer_list)
AX_CXX_HAVE_AUTOTYPE

SECTION_TITLE([Checking for zlib])

z_avail=0
AC_CHECK_HEADERS(zlib.h,[], [z_avail=1])
AC_SEARCH_LIBS(deflate, z, [], [z_avail=1])

if test $z_avail == 0
then
	AC_DEFINE(HAVE_ZLIB)
	AC_MSG_NOTICE([zlib available!])
fi

SECTION_TITLE([Checking for OpenMP])

omp=0
AX_OPENMP([omp=1])

if test $omp == 1
then
	AC_DEFINE(HAVE_OPENMP)
	CXXFLAGS+=" $OPENMP_CXXFLAGS"
	LIBS+=" $OPENMP_CXXFLAGS"
	
else
	AC_MSG_NOTICE([NO OMP])

fi

SECTION_TITLE([Configuration Results])
echo "CXXFLAGS=$CXXFLAGS" >& AS_MESSAGE_FD
echo "LDFLAGS=$LDFLAGS" >& AS_MESSAGE_FD
echo "LIBS=$LIBS" >& AS_MESSAGE_FD

echo
echo

echo "Optional components:"

if test $z_avail == 0
then
	echo 'WITH zlib: adios can seamlessly use gzipped files' >& AS_MESSAGE_FD
else
	echo "WITHOUT zlib: adios can not use gzipped files" >& AS_MESSAGE_FD
fi

if test $omp == 1
then
	echo "WITH OpenMP: adios can use multiple threads" >& AS_MESSAGE_FD
else
	echo "WITHOUT OpenMP: adios can only run single threaded."
	echo
	echo "This can happen in recent versions of Apple OSX, where the default install of clang does not include OpenMP" >& AS_MESSAGE_FD
	echo "On OSX try \`brew install clang-omp\` to install clang with OpenMP" >& AS_MESSAGE_FD
fi
echo
echo


dnl process config.h
AC_CONFIG_HEADERS(config.h)

dnl process Makefile.in
AC_OUTPUT(Makefile)
