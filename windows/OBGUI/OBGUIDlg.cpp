// OBGUIDlg.cpp : implementation file
//
#pragma warning (disable : 4786)

#include "stdafx.h"
#include <fstream>
#include <sstream>
#include <strstream>
#include <algorithm>
#include "OBGUI.h"
#include "OBGUIDlg.h"
#include "dlhandler.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace OpenBabel;
/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	CStatic	m_AboutCaption;
	CStatic	m_ExampleText;
	CStatic	m_AboutPlace;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	CDynamicOptions Opts;
	//{{AFX_MSG(CAboutDlg)
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	DDX_Control(pDX, IDC_ABOUTCAPTION, m_AboutCaption);
	DDX_Control(pDX, IDC_EXAMPLETEXT, m_ExampleText);
	DDX_Control(pDX, IDC_ABOUT_PLACEHOLDER, m_AboutPlace);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COBGUIDlg dialog

COBGUIDlg::COBGUIDlg(CWnd* pParent /*=NULL*/)
	: CDialog(COBGUIDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(COBGUIDlg)
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);

	//Window size info
	OrigConsoleWidth=340;
	OrigHeightDiff=170;
}

COBGUIDlg::~COBGUIDlg()
{
}


void COBGUIDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(COBGUIDlg)
	DDX_Control(pDX, IDC_OUTCOUNT, m_OutCount);
	DDX_Control(pDX, IDC_OUTFIXFORMAT, m_OutFixExt);
	DDX_Control(pDX, IDC_INFIXFORMAT, m_InFixExt);
	DDX_Control(pDX, IDC_CONSOLEINFO, m_ConsoleInfo);
	DDX_Control(pDX, IDC_MANUALINPUT, m_NoInFile);
	DDX_Control(pDX, IDC_INPUTBOX, m_InBox);
	DDX_Control(pDX, IDC_OUTPUTBOX, m_OutBox);
	DDX_Control(pDX, IDC_MANUALOUTPUT, m_NoOutFile);
	DDX_Control(pDX, IDC_OUPUTEDIT, m_OutConsole);
	DDX_Control(pDX, IDC_OUPUTFILE, m_OutputFile);
	DDX_Control(pDX, IDC_OUTPUTFORMAT, m_OutputCombo);
	DDX_Control(pDX, IDC_INPUTEDIT, m_InConsole);
	DDX_Control(pDX, IDC_INPUTFILE, m_InputFile);
	DDX_Control(pDX, IDC_PLACEHOLDER, m_Placeholder);
	DDX_Control(pDX, IDC_INPUTFORMAT, m_InputCombo);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(COBGUIDlg, CDialog)
	//{{AFX_MSG_MAP(COBGUIDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, OnInputFiles)
	ON_EN_CHANGE(IDC_INPUTFILE, OnChangeInputfile)
	ON_BN_CLICKED(IDC_IN_FORMAT_INFO, OnInFormatInfo)
	ON_BN_CLICKED(IDC_OUT_FORMAT_INFO, OnOutFormatInfo)
	ON_EN_CHANGE(IDC_OUPUTFILE, OnChangeOuputfile)
	ON_CBN_SELCHANGE(IDC_OUTPUTFORMAT, OnChangeOutputformat)
	ON_BN_CLICKED(IDC_CONVERT, OnConvert)
	ON_BN_CLICKED(IDC_BUTTON2, OnOutputFiles)
	ON_WM_SIZE()
	ON_BN_CLICKED(IDC_MANUALOUTPUT, OnManualoutput)
	ON_BN_CLICKED(IDC_MANUALINPUT, OnManualinput)
	ON_WM_CLOSE()
	ON_CBN_SELCHANGE(IDC_INPUTFORMAT, OnChangeInputformat)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COBGUIDlg message handlers

BOOL COBGUIDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
	m_NoInFile.SetCheck(AfxGetApp()->GetProfileInt("GUI","NoInFile",0));	
	OnManualinput();

	//Get data on available formats and add to comboboxes and to filter string
	char* str=NULL;
	OpenBabel::OBFormat* pFormat;
	int nInSel=0,nOutSel=0;
	InputFilterString="All Chemical Formats|*.";
	OutputFilterString = InputFilterString;
	Formatpos pos;
	while(Conv.GetNextFormat(pos,str,pFormat))
	{
		if(!str || !pFormat) break; //no formats available
		char* p = strstr(str," [");
		if(p) *p='\0'; //remove {Readonly] or [Writeonly]

		int n;
		CString txt=str;
		if(!(pFormat->Flags() & NOTREADABLE))
		{
			InputFilterString+=txt.Left(txt.Find(" "));
			InputFilterString+=";*.";
			n = m_InputCombo.AddString(str);
			m_InputCombo.SetItemDataPtr(n,pFormat);
			if(AfxGetApp()->GetProfileString("GUI","InFormat")==str)
				nInSel=n;
		}
		if(!(pFormat->Flags() & NOTWRITABLE))
		{
			OutputFilterString+=txt.Left(txt.Find(" "));
			OutputFilterString+=";*.";
			n = m_OutputCombo.AddString(str);
			m_OutputCombo.SetItemDataPtr(n,pFormat);
			if(AfxGetApp()->GetProfileString("GUI","OutFormat")==str)
				nOutSel=n;
		}		
	}
	m_InputCombo.SetCurSel(nInSel);
	m_OutputCombo.SetCurSel(nOutSel);


	InputFilterString = InputFilterString.Left(InputFilterString.GetLength()-3); //remove unneeded ;*.
	OutputFilterString = OutputFilterString.Left(OutputFilterString.GetLength()-3); //remove unneeded ;*.
	InputFilterString+="|AllFiles(*.*)|*.*||";
	OutputFilterString+="|AllFiles(*.*)|*.*||";
	
	CRect CtlRect, rect;
	m_Placeholder.GetWindowRect(CtlRect);
	GetWindowRect(rect);
	NextOptionRect = CtlRect-rect.TopLeft();

	//Construct checkboxes for General Options
	GeneralOptionCheckBoxes.Construct(Conv.Description(), this, NextOptionRect);	
	InputOptionRect=NextOptionRect;

	OnChangeInputformat(); //write input and then output options

	return TRUE;  // return TRUE  unless you set the focus to a control
}

//////////////////////////////////////////////////////////
void COBGUIDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}


// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void COBGUIDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR COBGUIDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

///////////////////////////////////////////
void COBGUIDlg::OnChangeOutputformat() 
{
	//Makes DynamicOptions Checkboxes appropriate for output format
	int nSel = m_OutputCombo.GetCurSel();
	if(nSel==CB_ERR) return;
	OBFormat* pFormat = (OBFormat*)m_OutputCombo.GetItemDataPtr(nSel);
	if(pFormat)
	{
		NextOptionRect=OutputOptionRect;
		InputOptionCheckBoxes.InsertText("", this, NextOptionRect); //Blank line
		OutputOptionCheckBoxes.Construct(pFormat->Description(), this, NextOptionRect);
	}
}

///////////////////////////////////////////
void COBGUIDlg::OnChangeInputfile() 
{
	const int BUFLEN=8000;
	char buf[BUFLEN];
	memset(buf,0,BUFLEN);
	CString FileName;
	m_InputFile.GetWindowText(FileName);
	FileName.TrimLeft();
	if(FileName.IsEmpty()) return;
	string stdfilename(FileName);
	
	vector<string> FileList;
	int nFiles = DLHandler :: findFiles (FileList, stdfilename);
	if(nFiles!=-1)
	{
		CString txt;
		txt.Format("** %d files found **\r\n",nFiles);
		vector<string>::iterator itr;
		for(itr=FileList.begin();itr!=FileList.end();itr++)
		{
			txt+= strrchr((*itr).c_str(),'\\')+1;
			txt+= "\r\n";
		}
		m_InConsole.SetWindowText(txt);
	}
	else
	{
		CFile iof;
		if(!iof.Open(FileName,CFile::modeRead |CFile::shareDenyNone ))
		{
			m_InConsole.SetWindowText("");
			return;
		}
		iof.Read(buf,BUFLEN);
		m_InConsole.SetWindowText(buf);
		iof.Close();
	}

	if(!m_InFixExt.GetCheck())
	{
		//Get extension
		OBFormat* pFormat = OBConversion::FormatFromExt(FileName);
		for(int iSel=0;iSel<m_InputCombo.GetCount();iSel++)
		{
			if(pFormat==m_InputCombo.GetItemDataPtr(iSel))
			{
				m_InputCombo.SetCurSel(iSel);
				break;
			}
		}
	}
}

////////////////////////////////////////
void COBGUIDlg::OnChangeOuputfile() 
{
	CString FileName;
	m_OutputFile.GetWindowText(FileName);
	FileName.TrimLeft();
	if(FileName.IsEmpty()) return;

	if(!m_OutFixExt.GetCheck())
	{
		//Get extension and adjust Format box
		OBFormat* pFormat = OBConversion::FormatFromExt(FileName);
		for(int iSel=0;iSel<m_OutputCombo.GetCount();iSel++)
		{
			if(pFormat==m_OutputCombo.GetItemDataPtr(iSel))
			{
				m_OutputCombo.SetCurSel(iSel);
				break;
			}
		}
		OnChangeOutputformat();
	}
}

///////////////////////////////////////////
void COBGUIDlg::OnInFormatInfo() 
{
	int nSel=m_InputCombo.GetCurSel();
	if(nSel<0) return;
	OBFormat* pFormat = (OBFormat*)m_InputCombo.GetItemDataPtr(nSel);
	CString mes =  pFormat->Description();
	CString url = pFormat->SpecificationURL();
	if(!url.IsEmpty())
		mes += "\nURL for specification: " + url;
	MessageBox(mes);		
}

///////////////////////////////////////////
void COBGUIDlg::OnOutFormatInfo() 
{
	int nSel=m_OutputCombo.GetCurSel();
	if(nSel<0) return;
	OBFormat* pFormat = (OBFormat*)m_OutputCombo.GetItemDataPtr(nSel);
	CString mes =  pFormat->Description();
	CString url = pFormat->SpecificationURL();
	if(!url.IsEmpty())
		mes += "\nURL for specification: " + url;
	MessageBox(mes);		
}

///////////////////////////////////////
void COBGUIDlg::OnConvert() 
{
	
	CString txt;
	m_InConsole.GetWindowText(txt);
	txt.GetLength();
	char* buf=txt.GetBuffer(txt.GetLength());//Default input streams from InConsole window
	strstream iss(buf,txt.GetLength());

	stringstream oss;//Default output stream. )Will be sent to OutConsole window)

	OBConversion Conv(&iss,&oss);

	int iSel = m_InputCombo.GetCurSel();
	if((iSel)<0) return;
	int oSel = m_OutputCombo.GetCurSel();
	if((oSel)<0) return;
	Conv.SetInAndOutFormats(
		(OBFormat*)m_InputCombo.GetItemDataPtr(iSel),
		(OBFormat*)m_OutputCombo.GetItemDataPtr(oSel));

	//GeneralOptions have data from GeneralOptionCheckboxes and InputOptionCheckboxes
	CString s = GeneralOptionCheckBoxes.GetOptions();
	s += InputOptionCheckBoxes.GetOptions();
	Conv.SetGeneralOptions(s);

	Conv.SetOptions(OutputOptionCheckBoxes.GetOptions());


  //redirect cerr
	stringstream serr;
  streambuf* cerr_sbuf = cerr.rdbuf();
  cerr.rdbuf(serr.rdbuf());
  // now everything written to 'std::cerr' is captured by serr

	CWaitCursor cw;

	CString FileName;
	m_OutputFile.GetWindowText(FileName);
	FileName.TrimLeft();
	string stdOutputFileName(FileName);

	//If you are trying to output with no filename what you really wanted
	//was to output to the OutConsole
	if(m_NoOutFile.GetCheck() || FileName.IsEmpty())
	{
		m_NoOutFile.SetCheck(1);
		OnManualoutput(); //disable file name box
		stdOutputFileName.erase();
	}

	vector<string> FileList, OutputFileList;
	if(m_NoInFile.GetCheck()==0)
	{
		m_InputFile.GetWindowText(FileName);
		FileName.TrimLeft();
		if(FileName.IsEmpty()) return;
		string stdfilename(FileName);

		DLHandler::findFiles(FileList,stdfilename);
	}

	int Count = Conv.FullConvert(FileList, stdOutputFileName, OutputFileList);

	CString mes;
	mes.Format("%d objects converted",Count);
	m_OutCount.SetWindowText(mes);

	StreamtoOutConsole(serr); //Clear window and display any errors 
	
	if(OutputFileList.size()>1)
	{
		CString oldtxt, txt;
		m_OutConsole.GetWindowText(oldtxt);
		if(!oldtxt.IsEmpty()) oldtxt+="\r\n";

		txt.Format("** %d files written **\r\n",OutputFileList.size());
		vector<string>::iterator itr;
		for(itr=OutputFileList.begin();itr!=OutputFileList.end();itr++)
		{
			txt+=(*itr).c_str();
			txt+= "\r\n\r\n";
		}
		txt += "The following is ";
		txt += OutputFileList[0].c_str();
		txt+= "\r\n\r\n";
		m_OutConsole.SetWindowText(oldtxt + txt);
	}

	if(Count>0)
	{
		if(m_NoOutFile.GetCheck()==0)
		{
			//Read back file and add to output console
			const int BUFLEN=8000;
			char buf[BUFLEN];
			memset(buf,0,BUFLEN);
			m_OutConsole.GetWindowText(buf,BUFLEN);

			if(!OutputFileList[0].empty())
			{
				CFile iof(OutputFileList[0].c_str(),CFile::modeRead |CFile::shareDenyNone );
				iof.Read(buf+strlen(buf),BUFLEN-strlen(buf)); //after current contents
				m_OutConsole.SetWindowText(buf);
				iof.Close();
			}
		}
		else
			StreamtoOutConsole(oss,true);
	}
  cerr.rdbuf(cerr_sbuf); // restore original stream buffer to avoid
                         // references to deleted objects
}

///////////////////////////////////////////////
void COBGUIDlg::StreamtoOutConsole(stringstream& oss, bool bAppend)
{
	string s=oss.str();

	//Change from lf to crlf needed for edit control apparently
	int pos=-2;
	while((pos=s.find(0x0a,pos+2))>=0)
		s.replace(pos,1,"\r\n");
	CString txt;
	if(bAppend)
		m_OutConsole.GetWindowText(txt);
	txt += s.c_str();
	m_OutConsole.SetWindowText(txt);
}
/////////////////////////////////////////////////
void COBGUIDlg::OnInputFiles() 
{
	//Entry for currently selected format

	CString Filter;
	GetFilter(Filter,m_InputCombo);
	CFileDialog dlg(TRUE,NULL,NULL,OFN_HIDEREADONLY,Filter+InputFilterString);
 	if(dlg.DoModal()==IDOK)
		m_InputFile.SetWindowText(dlg.GetPathName());
}

//////////////////////////////////////////////////
void COBGUIDlg::OnOutputFiles() 
{
	CString Filter;
	GetFilter(Filter,m_OutputCombo);
	CFileDialog dlg(FALSE,NULL,NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,Filter+OutputFilterString); 
 	if(dlg.DoModal()==IDOK)
		m_OutputFile.SetWindowText(dlg.GetPathName());	
}

///////////////////////////////////////////////////
void COBGUIDlg::GetFilter(CString& Filter, CWnd& Wnd)
{
	//Uses text from the window (input or output combo) to construct filter string
	CString txt;
	Wnd.GetWindowText(txt);
	Filter = txt;
	int n = txt.Find(" ");
	Filter = Filter.Mid(n);
	Filter.TrimLeft(" :-\t");
	Filter += " (*." + txt.Left(n) + ")|*." + txt.Left(n) + "|";
}

///////////////////////////////////////////////////
void COBGUIDlg::OnSize(UINT nType, int cx, int cy) 
{
	CDialog::OnSize(nType, cx, cy);
	//Makes consoles expand with main window
	int NewConsoleHeight=  + cy - OrigHeightDiff;
	if(::IsWindow(m_InConsole.m_hWnd))
	{
		const EXTRA=148;
		m_InConsole.SetWindowPos(NULL,0,0,OrigConsoleWidth,NewConsoleHeight,
			SWP_NOMOVE | SWP_NOZORDER | SWP_SHOWWINDOW);
		m_OutConsole.SetWindowPos(NULL,0,0,OrigConsoleWidth,NewConsoleHeight,
			SWP_NOMOVE | SWP_NOZORDER | SWP_SHOWWINDOW);
		m_InBox.SetWindowPos(NULL,0,0,OrigConsoleWidth,NewConsoleHeight + EXTRA,
			SWP_NOMOVE | SWP_NOZORDER | SWP_SHOWWINDOW);
		m_OutBox.SetWindowPos(NULL,0,0,OrigConsoleWidth,NewConsoleHeight+ EXTRA,
			SWP_NOMOVE | SWP_NOZORDER | SWP_SHOWWINDOW);
	}
}

////////////////////////////////////////////////
void COBGUIDlg::OnManualinput() 
{
	m_InputFile.EnableWindow(m_NoInFile.GetCheck()==0);
	if(m_NoInFile.GetCheck())
	{
		m_ConsoleInfo.ShowWindow(SW_SHOW);
		m_InConsole.SetFocus();
	}
	else
		m_ConsoleInfo.ShowWindow(SW_HIDE);
}

////////////////////////////////////////////////
void COBGUIDlg::OnManualoutput() 
{
	m_OutputFile.EnableWindow(m_NoOutFile.GetCheck()==0);	
}

////////////////////////////////////////////////
BOOL CAboutDlg::OnInitDialog() 
{
	CDialog::OnInitDialog();
	
char Caption[600] = " \
The application uses CDynamicOptions, which constructs a set of check,\n \
radio and edit boxes from text which looks like a description for\n \
command line options. The parsing is fairly flexible and the type of\n \
control is deduced from the text. Each option has a single letter;\n \
a call to GetOptions() returns a string with the letters of the active\n \
options and the edit text. The example returns  \
";
  
const char txt[]=	" \
This is some example text which is parsed to produce the controls below.\n \
Options (They start on the next line)\n \
a - A simple checkbox \n \
	b Another simple checkbox \n \
     c    This is checked (default appears)\n \
d Not checked because of no, not or none (default not)\n \
e This is the first mutually exclusive option (default) or\n \
f The next (note <or> is the last word) or\n \
g The last\n \
h<number> An editbox (punctuation, except '-' after option letter\n \
i# Another editbox default <Thurs>\n \
j# This edit is big because last character is :\n \
k The last option because the next line is blank, or end of text\n \
\n \
";								
	m_ExampleText.SetWindowText(txt);

	CRect Ctlrect, rect;
	GetWindowRect(rect);
	int Width=rect.Width();
	m_AboutPlace.GetWindowRect(Ctlrect);
	Width=Ctlrect.Width();
	Opts.Construct(txt,this,Ctlrect-rect.TopLeft());

	m_AboutCaption.SetWindowText(strcat(Caption,Opts.GetOptions()));
	return TRUE;  
}


void COBGUIDlg::OnClose() 
{
	//Save current format types in Registry to be loaded next startup
	CString txt;
	m_InputCombo.GetWindowText(txt);	
	AfxGetApp()->WriteProfileString("GUI","InFormat",txt);	
	m_OutputCombo.GetWindowText(txt);	
	AfxGetApp()->WriteProfileString("GUI","OutFormat",txt);	
	AfxGetApp()->WriteProfileInt("GUI","NoInFile",m_NoInFile.GetCheck());	
	CDialog::OnClose();
}

void COBGUIDlg::OnChangeInputformat() 
{
	//Makes DynamicOptions Checkboxes appropriate for input format
	int nSel = m_InputCombo.GetCurSel();
	if(nSel==CB_ERR) return;
	OBFormat* pFormat = (OBFormat*)m_InputCombo.GetItemDataPtr(nSel);
	if(pFormat)
	{
		NextOptionRect=InputOptionRect;
//		GeneralOptionCheckBoxes.InsertText("", this, NextOptionRect); //Blank line
		InputOptionCheckBoxes.Construct(pFormat->TargetClassDescription(), this, NextOptionRect);
	}

	OutputOptionRect=NextOptionRect;
	OnChangeOutputformat();  //now write output options
}

/*
/// If filename contains wildcard characters, populates FileList with
/// all the matching file names and returns the number found.
/// If no wildcards, copies filename to FileList[0] and returns -1 
int COBGUIDlg::TestWildcards(const char* filename, vector<string>* pFileList)
{
	pFileList->clear();
	pFileList->push_back(filename); //input file name, whatever it is 
	if(!strpbrk(filename,"*?")) return -1;
	pFileList->clear();

	WIN32_FIND_DATA file_data;
	HANDLE hF = FindFirstFile(filename,&file_data);
	if(hF==INVALID_HANDLE_VALUE) return 0;
	do
	{
		string fname = file_data.cFileName;
		pFileList->push_back(fname);
	}while(FindNextFile(hF,&file_data));
	return pFileList->size();
}

/////////////////////////////////////////////////
*/
