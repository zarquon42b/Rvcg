AC_INIT([RVCG], 0.6) 

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

dnl AX_CHECK_COMPILE_FLAG([-std=c++0x], [CXXFLAGS="$CXXFLAGS -std=c++0x"])
dnl AC_SUBST(CXXFLAGS)
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
