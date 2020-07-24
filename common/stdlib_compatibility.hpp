#ifndef STDLIB_COMPATIBLITY_HPP
#define STDLIB_COMPATIBLITY_HPP

#ifndef __APPLE__
  #if __cplusplus > 199711L
    #include <unordered_map>
    #define U_MAP std::unordered_map
  #else
    #include <tr1/unordered_map>
    #define U_MAP std::tr1::unordered_map
  #endif
#else
  #include <unordered_map>
  #define U_MAP std::unordered_map
#endif

#ifdef __clang__
  #include <random>
#else
  #if GCC_VERSION < 40300
    #include <tr1/random>
  #else
    #include <random>
  #endif
#endif

#endif // STDLIB_COMPATIBLITY_HPP
