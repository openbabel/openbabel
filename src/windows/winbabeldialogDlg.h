// winbabeldialogDlg.h : header file
//

#if !defined(AFX_WINBABELDIALOGDLG_H__2056478B_2421_11D6_9BEA_000000000000__INCLUDED_)
#define AFX_WINBABELDIALOGDLG_H__2056478B_2421_11D6_9BEA_000000000000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CWinbabeldialogDlg dialog

class CWinbabeldialogDlg : public CDialog
{
// Construction
public:
	CWinbabeldialogDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CWinbabeldialogDlg)
	enum { IDD = IDD_WINBABELDIALOG_DIALOG };
	CComboBox	m_cbOutputType;
	CComboBox	m_cbInputType;
	CString	m_sInputFile;
	CString	m_sInputType;
	CString	m_sOutputFile;
	CString	m_sOutputType;
	BOOL	m_bSplit;
	int		m_iHydrogens;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWinbabeldialogDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CWinbabeldialogDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnPickInput();
	afx_msg void OnPickOutput();
	afx_msg void OnConvert();
	afx_msg void OnHelp();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WINBABELDIALOGDLG_H__2056478B_2421_11D6_9BEA_000000000000__INCLUDED_)
