#ifndef CATAROB_INTERFACES__VISIBILITY_CONTROL_H_
#define CATAROB_INTERFACES__VISIBILITY_CONTROL_H_
#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define CATAROB_INTERFACES_EXPORT __attribute__ ((dllexport))
    #define CATAROB_INTERFACES_IMPORT __attribute__ ((dllimport))
  #else
    #define CATAROB_INTERFACES_EXPORT __declspec(dllexport)
    #define CATAROB_INTERFACES_IMPORT __declspec(dllimport)
  #endif
  #ifdef CATAROB_INTERFACES_BUILDING_LIBRARY
    #define CATAROB_INTERFACES_PUBLIC CATAROB_INTERFACES_EXPORT
  #else
    #define CATAROB_INTERFACES_PUBLIC CATAROB_INTERFACES_IMPORT
  #endif
  #define CATAROB_INTERFACES_PUBLIC_TYPE CATAROB_INTERFACES_PUBLIC
  #define CATAROB_INTERFACES_LOCAL
#else
  #define CATAROB_INTERFACES_EXPORT __attribute__ ((visibility("default")))
  #define CATAROB_INTERFACES_IMPORT
  #if __GNUC__ >= 4
    #define CATAROB_INTERFACES_PUBLIC __attribute__ ((visibility("default")))
    #define CATAROB_INTERFACES_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define CATAROB_INTERFACES_PUBLIC
    #define CATAROB_INTERFACES_LOCAL
  #endif
  #define CATAROB_INTERFACES_PUBLIC_TYPE
#endif
#endif  // CATAROB_INTERFACES__VISIBILITY_CONTROL_H_
// Generated 27-Oct-2024 09:47:09
 