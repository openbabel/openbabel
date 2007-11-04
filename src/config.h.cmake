/* src/config.h.in.  Generated from configure.in by autoheader.  */

/* Where the data files are located */
#define BABEL_DATADIR "@BABEL_DATADIR@"

/* The version of Open Babel */
#define BABEL_VERSION "@BABEL_VERSION@"

/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define EXTERN __declspec(dllimport) extern
 #else
  #define EXTERN __declspec(dllexport) extern
 #endif
#else //Everything else (behaviour as original)
 #define EXTERN extern 
#endif

/* The file extension used for shared modules */
#define MODULE_EXTENSION "@MODULE_EXTENSION@"

/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBAPI __declspec(dllimport)
 #else
  #define OBAPI __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBAPI 
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBCOMMON __declspec(dllimport)
 #else
  #define OBCOMMON __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBCOMMON 
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBCONV __declspec(dllimport)
 #else
  #define OBCONV __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBCONV
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBERROR __declspec(dllimport)
 #else
  #define OBERROR __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBERROR 
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(WIN32)
 #if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBFPTR __declspec(dllimport)
 #else
  #define OBFPTR __declspec(dllexport)
 #endif
#else //Everything else (behaviour as original)
 #define OBFPTR 
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
