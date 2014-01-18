AC_INIT([RVCG], 0.6) 
dnl m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
m4_include([m4/ax_cxx_compile_stdcxx_11.m4])

dnl AC_CONFIG_MACRO_DIR([m4])
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
echo "using C++ Compiler ...$CXX"
dnl CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_LANG_CPLUSPLUS
AX_CXX_COMPILE_STDCXX_11([noext],)

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

AC_DEFUN([AC_HEADER_STDCXX11_MAP], [
  AC_CACHE_CHECK(for ISO C++11 unordered_*.h include files,
  ac_cv_cxx_stdcxx11_map,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([
  #include <unordered_set>
  #include <unordered_map>
  
  ],,
  ac_cv_cxx_stdcxx11_map=yes, ac_cv_cxx_stdcxx11_map=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_stdcxx11_map" = yes; then
    AC_DEFINE(HAVE_CXX11_MAP,,[Define if ISO C++11 unordered_* header files are present. ])
  fi
])
AC_HEADER_STDCXX11_MAP

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

dnl AX_CHECK_COMPILE_FLAG([-std=gnu++11], [CXXFLAGS="$CXXFLAGS -std=gnu++11"])
AC_SUBST(CXXFLAGS)
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
