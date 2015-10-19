#include <config.h>
#ifdef WIN32
 #ifndef __MINGW32__
  #include <hash_map>
  #include <hash_set>
  #define STDEXT stdext
 #else
  #include <tr1/unordered_map>
  #include <tr1/unordered_set>
  #define  hash_map unordered_map
  #define  hash_set unordered_set
  #define  hash_multimap unordered_multimap
  #define STDEXT std::tr1
 #endif
#else
 #if defined HAVE_CXX11
  #include <unordered_map>
  #include <unordered_set>
  #define STDEXT std
  #define hash_map unordered_map
  #define hash_set unordered_set
  #define  hash_multimap unordered_multimap
 #elif defined HAVE_TR1
  #include <tr1/unordered_map>
  #include <tr1/unordered_set>
  #define STDEXT std::tr1
  #define  hash_map unordered_map
  #define  hash_set unordered_set
  #define  hash_multimap unordered_multimap
#else 
  #include <map>
  #include <set>
  #define STDEXT std
  #define hash_map map
  #define hash_set set
  #define hash_multimap multimap
 #endif
#endif
