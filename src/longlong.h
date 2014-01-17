#include <config.h>
#ifdef HAVE_INTTYPES
  #include <inttypes.h>
  #ifdef __int64
   #undef __int64
  #endif
  #define __int64 int64_t
  //#define __cdecl 
 #else
 #ifdef __int64
  #undef __int64
 #endif
  #define __int64 long long
  //#define __cdecl 
 #endif
