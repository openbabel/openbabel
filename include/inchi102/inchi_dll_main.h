/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * October 31, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
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