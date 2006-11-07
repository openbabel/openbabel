#ifndef OB_BCONFIG_H
#define OB_BCONFIG_H

//For use with Visual C++ (version 8)
#define HAVE_CONIO_H 1
#define HAVE_CLOCK_T 1
#define HAVE_IOSTREAM	1
#define HAVE_FSTREAM 1
#define HAVE_SSTREAM 1
#define HAVE_SNPRINTF 1
#define HAVE_STRNCASECMP 1
//#define HAVE_LIBZ 1 Causes DLL builds to not work with any file input!
#define BABEL_VERSION  "2.1"

#define BABEL_DATADIR "."
#define snprintf _snprintf
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#define rindex(a,b) strrchr((a),(b))

// Supress warning on deprecated functions
#pragma warning(disable : 4996)
// Supress warning that compiler is ignoring C++ exception specification
#pragma warning( disable : 4290 )
// Supress warning on signed/unsigned comparison with < or > (harmless, but maybe should be fixed)
#pragma warning( disable : 4018 )
// Supress warning on unreferenced formal parameter
#pragma warning( disable : 4100 )

#ifndef __FUNCTION__
#define __FUNCTION__ __FILE__
#endif

//Makes unix input files compatible with VC++6
#define ALL_READS_BINARY

/* Export of functions and global variables from DLLs
In the header files for the OB core, mol.h etc, exported classes and variables 
have an OBAPI declaration specifiers.
When building any OBDLL.dll: define OBDLL_EXPORTS.
When building an application that links with OBDLL.lib and so uses OBDLL.dll:
define USING_DYNAMIC_LIBS.

In obconversion.h the exported classes have OBCONV declaration specifiers.
When building OBConv.dll: define OBCONV_EXPORTS.

When building Format dlls: no need for any defines since no exported functions.
  
In non-Windows systems OBAPI and OBCONV need to be defined as empty.
*/
#if defined(USING_DYNAMIC_LIBS)
 #if defined(OBDLL_EXPORTS) //OBDLL being built
  #define OBAPI __declspec(dllexport)
 #else
  #define OBAPI __declspec(dllimport)
 #endif
#else //Everything else (behaviour as original)
 #define OBAPI
#endif

#if defined(USING_DYNAMIC_LIBS)
 #if defined(OBFPRT_EXPORTS) //OBFPRT being built
  #define OBFPRT __declspec(dllexport)
 #else
  #define OBFPRT __declspec(dllimport)
 #endif
#else //Everything else (behaviour as original)
 #define OBFPRT
#endif


#if defined(USING_DYNAMIC_LIBS)
 #pragma warning (disable : 4251) //no dll interface for some templated classes
 #ifdef OBCONV_EXPORTS
  #define OBCONV __declspec(dllexport)
 #else
  #define OBCONV __declspec(dllimport)
 #endif
#else
	#define OBCONV //as nothing in non-Windows system
#endif

#if defined(USING_DYNAMIC_LIBS)
 #if defined(OBCOMMON_EXPORTS) //OBCommonFormats being built
  #define OBCOMMON __declspec(dllexport)
 #else
  #define OBCOMMON __declspec(dllimport)
 #endif
#else //Everything else (behaviour as original)
 #define OBCOMMON
#endif

#if defined(USING_DYNAMIC_LIBS)
 #if defined(OBERROR_EXPORTS) //OBError being built
  #define OBERROR __declspec(dllexport)
 #else
  #define OBERROR __declspec(dllimport)
 #endif
#else //Everything else (behaviour as original)
 #define OBERROR
#endif

#if defined(OBDLL_EXPORTS) //OBDLL being built
#  define EXTERN __declspec(dllexport) extern
#elif defined(USING_OBDLL) //program using OBDLL.dll being built
#  define EXTERN __declspec(dllimport) extern
#else //Everything else (behaviour as original)
#  define EXTERN extern
#endif

/*
#ifdef _DEBUG
void* __cdecl operator new(size_t nSize, const char* lpszFileName, int nLine);
void __cdecl operator delete(void* p, const char* lpszFileName, int nLine);
#define DEBUG_NEW new(THIS_FILE, __LINE__)
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
*/

#include <crtdbg.h>

#ifdef _DEBUG
#define DEBUG_NEW new(_NORMAL_BLOCK, THIS_FILE, __LINE__)
#else
 #define DEBUG_NEW new
#endif


#endif //OB_BCONFIG_H