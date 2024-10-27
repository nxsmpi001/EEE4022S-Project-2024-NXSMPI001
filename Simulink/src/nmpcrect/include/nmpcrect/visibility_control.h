#ifndef NMPCRECT__VISIBILITY_CONTROL_H_
#define NMPCRECT__VISIBILITY_CONTROL_H_
#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define NMPCRECT_EXPORT __attribute__ ((dllexport))
    #define NMPCRECT_IMPORT __attribute__ ((dllimport))
  #else
    #define NMPCRECT_EXPORT __declspec(dllexport)
    #define NMPCRECT_IMPORT __declspec(dllimport)
  #endif
  #ifdef NMPCRECT_BUILDING_LIBRARY
    #define NMPCRECT_PUBLIC NMPCRECT_EXPORT
  #else
    #define NMPCRECT_PUBLIC NMPCRECT_IMPORT
  #endif
  #define NMPCRECT_PUBLIC_TYPE NMPCRECT_PUBLIC
  #define NMPCRECT_LOCAL
#else
  #define NMPCRECT_EXPORT __attribute__ ((visibility("default")))
  #define NMPCRECT_IMPORT
  #if __GNUC__ >= 4
    #define NMPCRECT_PUBLIC __attribute__ ((visibility("default")))
    #define NMPCRECT_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define NMPCRECT_PUBLIC
    #define NMPCRECT_LOCAL
  #endif
  #define NMPCRECT_PUBLIC_TYPE
#endif
#endif  // NMPCRECT__VISIBILITY_CONTROL_H_
// Generated 27-Oct-2024 09:47:09
 