/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * March 22, 2005
 * Developed at NIST
 */

/* inchi_dll_main.c : Defines the entry point for the DLL application. */
#if defined(_WIN32) && defined(_USRDLL) && defined(_DEBUG)
#include "inchi_dll_main.h"
//#include <windows.h>
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

