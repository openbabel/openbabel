//Optional initialisation routine for DLL
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#include <windows.h>

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
		char mes[] = "OBFormats.dll loaded";
    switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
			#ifdef _DEBUG
			::OutputDebugString(mes);
			#endif
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
    }
    return TRUE; 
};


