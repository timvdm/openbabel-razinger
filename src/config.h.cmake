/* src/config.h.in.  Generated from configure.in by autoheader.  */

/* Where the data files are located */
#define BABEL_DATADIR "@BABEL_DATADIR@"

/* The version of Open Babel */
#define BABEL_VERSION "@BABEL_VERSION@"

/* The file extension used for shared modules */
#define MODULE_EXTENSION "@MODULE_EXTENSION@"

// If we are using a recent GCC version with visibility support use it
#ifdef HAVE_GCC_VISIBILITY
  #define OB_EXPORT __attribute__ ((visibility("default")))
  #define OB_IMPORT __attribute__ ((visibility("default")))
  #define OB_HIDDEN __attribute__ ((visibility("hidden")))
#elif defined(WIN32) && !defined(__MINGW32__)
 #define OB_EXPORT __declspec(dllexport)
 #define OB_IMPORT __declspec(dllimport)
 #define OB_HIDDEN
#else
 #define OB_EXPORT
 #define OB_IMPORT
 #define OB_HIDDEN
#endif

/* Used to export symbols for DLL / shared library builds */
#if defined(MAKE_OBDLL) // e.g. in src/main.cpp
 #ifndef EXTERN
  #define EXTERN   OB_EXPORT extern
 #endif
 #ifndef OBAPI
  #define OBAPI    OB_EXPORT
 #endif
 #ifndef OBCOMMON
  #define OBCOMMON OB_EXPORT
 #endif
 #ifndef OBCONV
  #define OBCONV   OB_EXPORT
 #endif
 #ifndef OBERROR
  #define OBERROR  OB_EXPORT
 #endif
 #ifndef OBFPRT
  #define OBFPRT   OB_EXPORT
 #endif
 #ifndef OBFPTR
  #define OBFPTR   OB_EXPORT
 #endif
 #ifndef OBMCDL
  #define OBMCDL   OB_EXPORT
 #endif
 #ifndef OBDEPICT
  #define OBDEPICT OB_EXPORT
 #endif

#else   // defined(MAKE_OBDLL)

 #ifndef EXTERN
  #define EXTERN   OB_IMPORT extern
 #endif
 #ifndef OBAPI
  #define OBAPI    OB_IMPORT
 #endif
 #ifndef OBCOMMON
  #define OBCOMMON OB_IMPORT
 #endif
 #ifndef OBCONV
  #define OBCONV   OB_IMPORT
 #endif
 #ifndef OBERROR
  #define OBERROR  OB_IMPORT
 #endif
 #ifndef OBFPRT
  #define OBFPRT   OB_IMPORT
 #endif
 #ifndef OBFPTR
  #define OBFPTR   OB_IMPORT
 #endif
 #ifndef OBMCDL
 #define OBMCDL    OB_IMPORT
  #ifndef OBDEPICT
 #define OBDEPICT  OB_IMPORT
 #endif
 
 #endif

#endif

#ifdef _MSC_VER
 // Supress warning on deprecated functions
 #pragma warning(disable : 4996)
 // Supress warning that compiler is ignoring C++ exception specification
 #pragma warning( disable : 4290 )
 // Supress warning on signed/unsigned comparison with < or > (harmless, but maybe should be fixed)
 #pragma warning( disable : 4018 )
 // Supress warning on unreferenced formal parameter
 #pragma warning( disable : 4100 )
 //
 #pragma warning( disable : 4251 )


 #include <crtdbg.h>

 #ifdef _DEBUG
 #define DEBUG_NEW new(_NORMAL_BLOCK, THIS_FILE, __LINE__)
 #else
  #define DEBUG_NEW new
 #endif
#endif  // _MSC_VER
/* have <conio.h> */
#cmakedefine HAVE_CONIO_H 1

/* have <sys/time.h> */
#cmakedefine HAVE_SYS_TIME_H 1

/* have <time.h> */
#cmakedefine HAVE_TIME_H 1

/* have <sstream> */
#cmakedefine HAVE_SSTREAM 1

/* have symbol clock_t */
#cmakedefine HAVE_CLOCK_T 1

/* have symbol rint */
#cmakedefine HAVE_RINT 1

/* have symbol snprintf */
#cmakedefine HAVE_SNPRINTF 1

/* have symbol sranddev */
#cmakedefine HAVE_SRANDDEV 1

/* have symbol strcasecmp */
#cmakedefine HAVE_STRCASECMP 1

/* have symbol strncasecmp */
#cmakedefine HAVE_STRNCASECMP 1

/* have struct clock_t */
#cmakedefine HAVE_CLOCK_T 1

#if defined(WIN32)
 #ifndef HAVE_ISFINITE
  #define isfinite _finite
  #define HAVE_ISFINITE 1
 #endif

 #ifndef HAVE_SNPRINTF
  #define snprintf _snprintf
  #define HAVE_SNPRINTF 1
 #endif

 #ifndef HAVE_STRCASECMP
  #define strcasecmp _stricmp
  #define HAVE_STRCASECMP 1
 #endif

 #ifndef HAVE_STRNCASECMP
  #define strncasecmp _strnicmp
  #define HAVE_STRNCASECMP 1
 #endif
#endif  // WIN32

#cmakedefine SCANDIR_NEEDS_CONST 1
#ifdef SCANDIR_NEEDS_CONST
 #define SCANDIR_CONST const
#else
 #define SCANDIR_CONST
#endif

#define OB_MODULE_PATH "@OB_MODULE_PATH@"

#ifndef TIME_WITH_SYS_TIME
  #ifdef HAVE_SYS_TIME
    #ifdef HAVE_TIME
      #define TIME_WITH_SYS_TIME 1
    #else
      #define TIME_WITH_SYS_TIME 0
    #endif
  #else
    #define TIME_WITH_SYS_TIME 0
  #endif
#endif

