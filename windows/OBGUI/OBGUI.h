// OBGUI.h : main header file for the OBGUI application
//

#if !defined(AFX_OBGUI_H__F55EC011_D348_43F6_9EAE_A90C18AA7852__INCLUDED_)
#define AFX_OBGUI_H__F55EC011_D348_43F6_9EAE_A90C18AA7852__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// COBGUIApp:
// See OBGUI.cpp for the implementation of this class
//

class COBGUIApp : public CWinApp
{
public:
	COBGUIApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COBGUIApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(COBGUIApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OBGUI_H__F55EC011_D348_43F6_9EAE_A90C18AA7852__INCLUDED_)
