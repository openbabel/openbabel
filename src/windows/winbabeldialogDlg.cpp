// winbabeldialogDlg.cpp : implementation file
//

#include "stdafx.h"
#include "winbabeldialog.h"
#include "winbabeldialogDlg.h"

extern int NewMain(char*, char*, char*, char*, int, bool);

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
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
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWinbabeldialogDlg dialog

CWinbabeldialogDlg::CWinbabeldialogDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CWinbabeldialogDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CWinbabeldialogDlg)
	m_sInputFile = _T("");
	m_sInputType = _T("");
	m_sOutputFile = _T("");
	m_sOutputType = _T("");
	m_bSplit = FALSE;
	m_iHydrogens = -1;
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CWinbabeldialogDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CWinbabeldialogDlg)
	DDX_Control(pDX, IDC_OUTPUT_TYPE, m_cbOutputType);
	DDX_Control(pDX, IDC_INPUT_TYPE, m_cbInputType);
	DDX_Text(pDX, IDC_INPUT_FILE, m_sInputFile);
	DDX_CBString(pDX, IDC_INPUT_TYPE, m_sInputType);
	DDX_Text(pDX, IDC_OUTPUT_FILE, m_sOutputFile);
	DDX_CBString(pDX, IDC_OUTPUT_TYPE, m_sOutputType);
	DDX_Check(pDX, IDC_SPLIT, m_bSplit);
	DDX_Radio(pDX, IDC_H_NOCHANGE, m_iHydrogens);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CWinbabeldialogDlg, CDialog)
	//{{AFX_MSG_MAP(CWinbabeldialogDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_PICK_INPUT, OnPickInput)
	ON_BN_CLICKED(IDC_PICK_OUTPUT, OnPickOutput)
	ON_BN_CLICKED(IDCONVERT, OnConvert)
	ON_BN_CLICKED(IDHELP, OnHelp)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWinbabeldialogDlg message handlers

BOOL CWinbabeldialogDlg::OnInitDialog()
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
	
	// TODO: Add extra initialization here
	// init variables
	m_sInputFile = "";
	m_sOutputFile = "";
	m_bSplit = FALSE;
	m_iHydrogens = 0;

	// populate the combo boxes
	m_cbInputType.AddString("alc -- Alchemy file");
	m_cbInputType.AddString("prep -- Amber PREP file");
	m_cbInputType.AddString("bs -- Ball & Stick file");
	m_cbInputType.AddString("caccrt -- Cacao Cartesian file");
	m_cbInputType.AddString("ccc -- CCC file");
	m_cbInputType.AddString("box -- Dock 3.5 Box file");
	m_cbInputType.AddString("dmol -- DMol3 coordinates file");
	m_cbInputType.AddString("feat -- Feature file");
	m_cbInputType.AddString("gam -- GAMESS Output file");
	m_cbInputType.AddString("gamout -- GAMESS Output file");
	m_cbInputType.AddString("mm1gp -- Ghemical MM file");
	m_cbInputType.AddString("qm1gp -- Ghemical QM file");
	m_cbInputType.AddString("hin -- HyperChem HIN file");
	m_cbInputType.AddString("jout -- Jaguar Output file");
	m_cbInputType.AddString("bin -- OpenEye binary file");
	m_cbInputType.AddString("mmd -- MacroModel file");
	m_cbInputType.AddString("mmod -- MacroModel file");
	m_cbInputType.AddString("out -- MacroModel file");
	m_cbInputType.AddString("dat -- MacroModel file");
	m_cbInputType.AddString("car -- MSI Biosym/Insight II file");
	m_cbInputType.AddString("sdf -- MDL Isis SDF file");
	m_cbInputType.AddString("sd -- MDL ISIS SDF file");
	m_cbInputType.AddString("mdl -- MDL Molfile file");
	m_cbInputType.AddString("mol -- MDL Molfile file");
	m_cbInputType.AddString("mopcrt -- MOPAC Cartesian file");
	m_cbInputType.AddString("mopout -- MOPAC Output file");
	m_cbInputType.AddString("mpqc -- MPQC file");
	m_cbInputType.AddString("bgf -- MSI BGF file");
	m_cbInputType.AddString("nwo -- NWChem Output file");
	m_cbInputType.AddString("pdb -- PDB file");
	m_cbInputType.AddString("qcout -- QChem Output file");
	m_cbInputType.AddString("smi -- SMILES file");
	m_cbInputType.AddString("mol2 -- Sybyl Mol2 file");
	m_cbInputType.AddString("unixyz -- UniChem XYZ file");
	m_cbInputType.AddString("xyz -- XYZ file");
	m_cbOutputType.AddString("alc -- Alchemy file");
	//m_cbOutputType.AddString("prep -- Amber PREP file");
	m_cbOutputType.AddString("bs -- Ball & Stick file");
	m_cbOutputType.AddString("caccrt -- Cacao Cartesian file");
	m_cbOutputType.AddString("cacint -- Cacao Internal file");
	m_cbOutputType.AddString("cache -- CAChe MolStruct file");
	m_cbOutputType.AddString("ct -- ChemDraw Connection Table file");
	m_cbOutputType.AddString("cssr -- CSD CSSR file");
	//m_cbOutputType.AddString("ccc -- CCC file");
	m_cbOutputType.AddString("box -- Dock 3.5 Box file");
	m_cbOutputType.AddString("dmol -- DMol3 coordinates file");
	m_cbOutputType.AddString("feat -- Feature file");
	m_cbOutputType.AddString("fh -- Fenske-Hall Z-Matrix file");
	m_cbOutputType.AddString("gamin -- GAMESS Input file");
	m_cbOutputType.AddString("inp -- GAMESS Input file");
	m_cbOutputType.AddString("gcart -- Gaussian Cartesian file");
	m_cbOutputType.AddString("gau -- Gaussian Input file");
	m_cbOutputType.AddString("mm1gp -- Ghemical MM file");
	//m_cbOutputType.AddString("qm1gp -- Ghemical QM file");
	m_cbOutputType.AddString("gr96A -- Gromos (A) file");
	m_cbOutputType.AddString("gr96N -- Gromos (nm) file");
	m_cbOutputType.AddString("hin -- HyperChem HIN file");
	m_cbOutputType.AddString("jin -- Jaguar Input file");
	m_cbOutputType.AddString("bin -- OpenEye binary file");
	m_cbOutputType.AddString("mmd -- MacroModel file");
	m_cbOutputType.AddString("mmod -- MacroModel file");
	m_cbOutputType.AddString("out -- MacroModel file");
	m_cbOutputType.AddString("dat -- MacroModel file");
	//m_cbOutputType.AddString("car -- MSI Biosym/Insight II file");
	m_cbOutputType.AddString("sdf -- MDL Isis SDF file");
	m_cbOutputType.AddString("sd -- MDL ISIS SDF file");
	m_cbOutputType.AddString("mdl -- MDL Molfile file");
	m_cbOutputType.AddString("mol -- MDL Molfile file");
	m_cbOutputType.AddString("mopcrt -- MOPAC Cartesian file");
	//m_cbOutputType.AddString("mopout -- MOPAC Output file");
	//m_cbOutputType.AddString("mpqc -- MPQC file");
	m_cbOutputType.AddString("bgf -- MSI BGF file");
	m_cbOutputType.AddString("csr -- MSI Quanta CSR file");
	m_cbOutputType.AddString("nw -- NWChem Input file");
	m_cbOutputType.AddString("pdb -- PDB file");
	m_cbOutputType.AddString("report -- Report file");
	m_cbOutputType.AddString("qcin -- QChem Input file");
	m_cbOutputType.AddString("smi -- SMILES file");
	m_cbOutputType.AddString("fix -- SMILES Fix file");
	m_cbOutputType.AddString("mol2 -- Sybyl Mol2 file");
	m_cbOutputType.AddString("txyz -- Tinker XYZ file");
	m_cbOutputType.AddString("txt -- Titles file");
	m_cbOutputType.AddString("unixyz -- UniChem XYZ file");
	m_cbOutputType.AddString("xed -- XED file");
	m_cbOutputType.AddString("xyz -- XYZ file");


	UpdateData(FALSE);

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CWinbabeldialogDlg::OnSysCommand(UINT nID, LPARAM lParam)
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

void CWinbabeldialogDlg::OnPaint() 
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
HCURSOR CWinbabeldialogDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CWinbabeldialogDlg::OnPickInput() 
{
	CString ext;

	// TODO: Add your control notification handler code here
	CFileDialog m_ldFile(TRUE);
	if (m_ldFile.DoModal() == IDOK) {
		m_sInputFile = m_ldFile.GetPathName();
		UpdateData(FALSE);
		ext = m_ldFile.GetFileExt();
		if (ext == "alc")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "alc "));
		if (ext == "mol")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "mol "));
		if (ext == "mol2")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "mol2 "));
		if (ext == "smi")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "smi "));
		if (ext == "pdb")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "pdb "));
		if (ext == "xyz")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "xyz "));
		if (ext == "mmd")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "mmd "));
		if (ext == "mmod")
			m_cbInputType.SetCurSel(m_cbInputType.FindString(0, "mmod "));
		//m_sInputType = m_cbInputType.GetItemData();
	}
}

void CWinbabeldialogDlg::OnPickOutput() 
{
	// TODO: Add your control notification handler code here
	CString ext;

	CFileDialog m_ldFile(TRUE);
	if (m_ldFile.DoModal() == IDOK) {
		m_sOutputFile = m_ldFile.GetPathName();
		UpdateData(FALSE);
		ext = m_ldFile.GetFileExt();
		if (ext == "alc")
			m_cbOutputType.SelectString(0, "alc ");
		if (ext == "mol")
			m_cbOutputType.SelectString(0, "mol ");
		if (ext == "mol2")
			m_cbOutputType.SelectString(0, "mol2 ");
		if (ext == "smi")
			m_cbOutputType.SelectString(0, "smi ");
		if (ext == "pdb")
			m_cbOutputType.SelectString(0, "pdb ");
		if (ext == "xyz")
			m_cbOutputType.SelectString(0, "xyz ");
		if (ext == "mmd")
			m_cbOutputType.SelectString(0, "mmd ");
		if (ext == "mmod")
			m_cbOutputType.SelectString(0, "mmod ");
		//m_sOutputType = m_cbOutputType.GetItemData();
	}	
}

void CWinbabeldialogDlg::OnConvert() 
{
	// get data from dialog objects
	UpdateData(TRUE);	

	if (m_sInputFile == "") {
		MessageBox("Invalid input file name.", "Error", MB_ICONWARNING | MB_OK);
		return;
	}
	if (m_sInputType == "") {
		MessageBox("Invalid input file type.  Select from the combo box.", "Error", MB_ICONWARNING | MB_OK);
		return;
	}
	if (m_sOutputFile == "") {
		MessageBox("Invalid output file name.", "Error", MB_ICONWARNING | MB_OK);
		return;
	}
	if (m_sOutputType == "") {
		MessageBox("Invalid output file type.  Select from the combo box", "Error", MB_ICONWARNING | MB_OK);
		return;
	}

	// fix type selection
	m_sInputType = m_sInputType.Left(m_sInputType.Find(" --"));
	m_sOutputType = m_sOutputType.Left(m_sOutputType.Find(" --"));

	int x;
	char *cInFile = (char*)malloc(sizeof(char) * (m_sInputFile.GetLength() + 1));
	for (x = 0; x < m_sInputFile.GetLength(); x++) {
		cInFile[x] = m_sInputFile.GetAt(x);
	}
	cInFile[x] = 0;
	char *cInType = (char*)malloc(sizeof(char) * (m_sInputType.GetLength() + 1));
	for (x = 0; x < m_sInputType.GetLength(); x++) {
		cInType[x] = m_sInputType.GetAt(x);
	}
	cInType[x] = 0;
	char *cOutFile = (char*)malloc(sizeof(char) * (m_sOutputFile.GetLength() + 1));
	for (x = 0; x < m_sOutputFile.GetLength(); x++) {
		cOutFile[x] = m_sOutputFile.GetAt(x);
	}
	cOutFile[x] = 0;
	char *cOutType = (char*)malloc(sizeof(char) * (m_sOutputType.GetLength() + 1));
	for (x = 0; x < m_sOutputType.GetLength(); x++) {
		cOutType[x] = m_sOutputType.GetAt(x);
	}
	cOutType[x] = 0;
	int suc = NewMain(cInFile, cInType, cOutFile, cOutType, m_iHydrogens, m_bSplit);
	if (suc > 0)
		MessageBox("Conversion failed! Reason = " + sprintf(0,"%d",suc), "Error", MB_ICONWARNING | MB_OK);
	else
		MessageBox("Successful conversion.", "Success!", MB_ICONINFORMATION | MB_OK);
}

void CWinbabeldialogDlg::OnHelp() 
{
	// TODO: Add your control notification handler code here
	CAboutDlg about;
	
	about.DoModal();
}
