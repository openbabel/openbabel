// winbabeldialog.h : main header file for the WINBABELDIALOG application
//

#if !defined(AFX_WINBABELDIALOG_H__20564789_2421_11D6_9BEA_000000000000__INCLUDED_)
#define AFX_WINBABELDIALOG_H__20564789_2421_11D6_9BEA_000000000000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
  #error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CWinbabeldialogApp:
// See winbabeldialog.cpp for the implementation of this class
//

class CWinbabeldialogApp : public CWinApp
{
public:
	CWinbabeldialogApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWinbabeldialogApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CWinbabeldialogApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WINBABELDIALOG_H__20564789_2421_11D6_9BEA_000000000000__INCLUDED_)
