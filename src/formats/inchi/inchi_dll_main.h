/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * May 16, 2006
 * Developed at NIST
 */

#ifndef __INCHI_DLL_MAIN_H__
#define __INCHI_DLL_MAIN_H__

#if _MSC_VER > 1000
#pragma once
#endif /* _MSC_VER > 1000 */

#if defined(_WIN32) && defined(_MSC_VER) && defined(_USRDLL)

/*#define WIN32_LEAN_AND_MEAN */  /* Exclude rarely-used stuff from Windows headers */
#include <windows.h>

#define  INCHI_DLLMAIN_TYPE APIENTRY

#else  /* not a Win32 DLL under MS VC++ */

#define  INCHI_DLLMAIN_TYPE

#endif


#endif /* __INCHI_DLL_MAIN_H__ */