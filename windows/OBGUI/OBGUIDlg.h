// OBGUIDlg.h : header file
//

#if !defined(AFX_OBGUIDLG_H__6C28096C_05FC_4EA6_A521_BC59869475CB__INCLUDED_)
#define AFX_OBGUIDLG_H__6C28096C_05FC_4EA6_A521_BC59869475CB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <fstream>
#include <sstream>
#include "obconversion.h"
#include "DynamicOptions.h"

/////////////////////////////////////////////////////////////////////////////
// COBGUIDlg dialog

class COBGUIDlg : public CDialog
{
// Construction
public:
	COBGUIDlg(CWnd* pParent = NULL);	// standard constructor
	~COBGUIDlg();
// Dialog Data
	//{{AFX_DATA(COBGUIDlg)
	enum { IDD = IDD_OBGUI_DIALOG };
	CStatic	m_InPath;
	CStatic	m_OutCount;
	CButton	m_OutFixExt;
	CButton	m_InFixExt;
	CStatic	m_ConsoleInfo;
	CButton	m_NoInFile;
	CButton	m_InBox;
	CButton	m_OutBox;
	CButton	m_NoOutFile;
	CEdit	m_OutConsole;
	CEdit	m_OutputFile;
	CComboBox	m_OutputCombo;
	CEdit	m_InConsole;
	CEdit	m_InputFile;
	CStatic	m_Placeholder;
	CComboBox	m_InputCombo;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COBGUIDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;
	int OrigHeightDiff;
	int OrigConsoleWidth;
	CString InputFilterString, OutputFilterString;
	stringstream newerr;
	streambuf* Origclog;

	CDynamicOptions InputClassOptionCheckBoxes, InputOptionCheckBoxes, OutputOptionCheckBoxes,
		GeneralOptionCheckBoxes;
	//OpenBabel::OBConversion Conv;
  ifstream inFileStream;
	ofstream outFileStream;
	CRect NextOptionRect; //Placeholder for next dynamic option; relative to parent window.
	CRect OutputOptionRect, InputOptionRect; //for when in or out format chages

	void GetFilter(CString& Filter, CWnd& Wnd);
	void StreamtoOutConsole(stringstream& oss, bool bAppend=false);
	void DisplayPath(char* path);

	// Generated message map functions
	//{{AFX_MSG(COBGUIDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnInputFiles();
	afx_msg void OnChangeInputfile();
	afx_msg void OnInFormatInfo();
	afx_msg void OnOutFormatInfo();
	afx_msg void OnChangeOuputfile();
	afx_msg void OnChangeOutputformat();
	afx_msg void OnConvert();
	afx_msg void OnOutputFiles();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnManualoutput();
	afx_msg void OnManualinput();
	afx_msg void OnClose();
	afx_msg void OnChangeInputformat();
	afx_msg void OnDropFiles( HDROP hDropInfo );
//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OBGUIDLG_H__6C28096C_05FC_4EA6_A521_BC59869475CB__INCLUDED_)
