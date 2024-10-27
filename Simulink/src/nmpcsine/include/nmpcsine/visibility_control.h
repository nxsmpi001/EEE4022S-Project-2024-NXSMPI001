#ifndef NMPCSINE__VISIBILITY_CONTROL_H_
#define NMPCSINE__VISIBILITY_CONTROL_H_
#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define NMPCSINE_EXPORT __attribute__ ((dllexport))
    #define NMPCSINE_IMPORT __attribute__ ((dllimport))
  #else
    #define NMPCSINE_EXPORT __declspec(dllexport)
    #define NMPCSINE_IMPORT __declspec(dllimport)
  #endif
  #ifdef NMPCSINE_BUILDING_LIBRARY
    #define NMPCSINE_PUBLIC NMPCSINE_EXPORT
  #else
    #define NMPCSINE_PUBLIC NMPCSINE_IMPORT
  #endif
  #define NMPCSINE_PUBLIC_TYPE NMPCSINE_PUBLIC
  #define NMPCSINE_LOCAL
#else
  #define NMPCSINE_EXPORT __attribute__ ((visibility("default")))
  #define NMPCSINE_IMPORT
  #if __GNUC__ >= 4
    #define NMPCSINE_PUBLIC __attribute__ ((visibility("default")))
    #define NMPCSINE_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define NMPCSINE_PUBLIC
    #define NMPCSINE_LOCAL
  #endif
  #define NMPCSINE_PUBLIC_TYPE
#endif
#endif  // NMPCSINE__VISIBILITY_CONTROL_H_
// Generated 27-Oct-2024 09:43:03
 