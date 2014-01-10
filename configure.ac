AC_INIT([RVCG], 0.6) 
#PKG_CXXFLAGS=`$(R_HOME)/bin/Rscript -e "Rcpp:::CxxFlags()"` -std=c++0x
#PKG_CXXFLAGS=-std=c++0x
#     AC_CONFIG_AUX_DIR([subdir])
AX_CHECK_COMPILE_FLAG([-std=c++0x], [CXXFLAGS="$CXXFLAGS -std=c++0x"])
AC_OUTPUT
AC_CANONICAL_HOST

have_windows="no"
case "$host_os" in
   mingw*)
      have_windows="yes"
   ;;
esac 
if test "${have_windows}" = no;then
   PKG_LIBS='`$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"`'
else 				    
   PKG_LIBS='$(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript.exe" -e "Rcpp:::LdFlags()")'

fi
#AX_CXX_COMPILE_STDCXX_0X
AC_SUBST(CXXFLAGS)
AC_SUBST(PKG_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
