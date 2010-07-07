/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
 */


/* inchi_dll_main.c : Defines the entry point for the DLL application. */

#if defined(_WIN32) && defined(_USRDLL) && defined(_DEBUG) && !(defined(__STDC__) && __STDC__ == 1)
#include "inchi_dll_main.h"
int INCHI_DLLMAIN_TYPE DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
                     )
{
    return TRUE;
}
#else
int dummy_inchi_dll_main=0;  /* avoid empty module to keep C compiler happy */
#endif

