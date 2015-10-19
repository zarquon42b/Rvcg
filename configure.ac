AC_INIT([RVCG], 0.6) 

m4_include([src/m4/ax_cxx_compile_stdcxx_11.m4])
dnl AC_CONFIG_MACRO_DIR([m4])
#get R specific compiler and flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
echo "using C++ Compiler ...$CXX"
dnl CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_LANG_CPLUSPLUS
# check C++11 compliance and flags
AX_CXX_COMPILE_STDCXX_11([noext],)

# check availability of TR1 headers
AC_DEFUN([AC_HEADER_STDCXX_TR1], [
  AC_CACHE_CHECK(for ISO C++ TR1 include files,
  ac_cv_cxx_stdcxx_tr1,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([
  #include <tr1/unordered_set>
  #include <tr1/unordered_map>
  ],,
  ac_cv_cxx_stdcxx_tr1=yes, ac_cv_cxx_stdcxx_tr1=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_stdcxx_tr1" = yes; then
    AC_DEFINE(HAVE_TR1,,[Define if ISO C++ TR1 header files are present. ])
  fi
])
AC_HEADER_STDCXX_TR1

# check for inttypes.h
AC_DEFUN([AC_HEADER_STDCXX_INTTYPES], [
  AC_CACHE_CHECK(for inttypes.h,
  ac_cv_cxx_stdcxx_inttypes,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([
  #include <inttypes.h>
  ],,
  ac_cv_cxx_stdcxx_inttypes=yes, ac_cv_cxx_stdcxx_inttypes=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_stdcxx_inttypes" = yes; then
    AC_DEFINE(HAVE_INTTYPES,,[Define if inttypes.h header files are present. ])
  fi
])
AC_HEADER_STDCXX_INTTYPES

AC_SUBST(CXXFLAGS)
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
