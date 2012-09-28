/**********************************************************************
OBGUI.cpp -  Implementation of a cross-platform Graphical User Interface

Copyright (C) 2006 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#define HAVE_SNPRINTF_DECL //needed to avoid "inconsistent linkage"
#include <openbabel/babelconfig.h>
#include <stdwx.h>
#include <wx/file.h>
//#include <sstream>
#include <openbabel/plugin.h>
#include <openbabel/obconversion.h>
//#include <openbabel/dlhandler.h>
#include <selformats.h>
#include <OBGUI.h>
#include <wx/splash.h>
#include <wx/imagpng.h>
#include <openbabel/tokenst.h>


#ifdef __WXMAC__
// As described at http://wiki.wxwidgets.org/WxMac_Issues
#include <ApplicationServices/ApplicationServices.h>
#endif

/*
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
*/
using namespace OpenBabel;

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them.
BEGIN_EVENT_TABLE(OBGUIFrame, wxFrame)
  EVT_MENU(ID_CONVERT,     OBGUIFrame::OnConvert)
  EVT_MENU(wxID_EXIT,      OBGUIFrame::OnQuit)
  EVT_MENU(wxID_SAVE,      OBGUIFrame::OnSaveInputText)
  EVT_MENU(ID_COPYTOINPUT, OBGUIFrame::OnCopyToInput)
  EVT_MENU(ID_SAVECONFIG,  OBGUIFrame::OnSaveConfig)
  EVT_MENU(ID_SELFORMATS,  OBGUIFrame::OnSelectFormats)
  EVT_MENU(ID_RESTRICTFORMATS,  OBGUIFrame::OnRestrictFormats)
  EVT_MENU(ID_SETDISPLAYFILE,  OBGUIFrame::OnSetDisplayFile)
  EVT_MENU(ID_MINSIZE,  OBGUIFrame::OnSetMinSize)
  EVT_MENU_RANGE(ID_PLUGINS,ID_PLUGINS+1000, OBGUIFrame::OnClickPlugin)
  EVT_MENU_RANGE(ID_SHOWCONVOPTIONS,ID_SHOWOUTOPTIONS, OBGUIFrame::OnChangeFormat)
  EVT_MENU(wxID_ABOUT, OBGUIFrame::OnAbout)
  EVT_MENU(wxID_HELP, OBGUIFrame::OnHelp)
//  EVT_UPDATE_UI(ID_PLUGINS, OBGUIFrame::OnPlugins)
  EVT_BUTTON(ID_INGETFILES, OBGUIFrame::OnGetInputFile)
  EVT_BUTTON(ID_OUTGETFILES, OBGUIFrame::OnGetOutputFile)
  EVT_BUTTON(ID_ININFO, OBGUIFrame::OnInFormatInfo)
  EVT_BUTTON(ID_OUTINFO, OBGUIFrame::OnOutFormatInfo)
  EVT_UPDATE_UI_RANGE(ID_OUTFILENAME,ID_OUTGETFILES, OBGUIFrame::OnOutFileNameUpdate)
  EVT_UPDATE_UI_RANGE(ID_INFILENAME, ID_INGETFILES, OBGUIFrame::OnInFileNameUpdate)
  EVT_BUTTON(ID_CONVERT, OBGUIFrame::OnConvert)
  EVT_CHECKBOX(ID_INPUTHERE,OBGUIFrame::OnChangeInputHere)
  EVT_CHOICE(ID_INFORMAT,OBGUIFrame::OnChangeFormat)
  EVT_CHOICE(ID_OUTFORMAT,OBGUIFrame::OnChangeFormat)
  EVT_MOUSEWHEEL(OBGUIFrame::OnMouseWheel)
  EVT_CLOSE(OBGUIFrame::OnClose)
 END_EVENT_TABLE()

IMPLEMENT_APP(OBGUIApp)

// 'Main program' equivalent: the program execution "starts" here
bool OBGUIApp::OnInit()
{
#ifdef __WXMAC__
  // As described at http://wiki.wxwidgets.org/WxMac_Issues
  ProcessSerialNumber PSN;
  GetCurrentProcess(&PSN);
  TransformProcessType(&PSN,kProcessTransformToForegroundApplication);
#endif

  OBConversion dummy; //needed for OBConversion to load plugin classes (including formats)
  //Read in the stored window sizes and extensions previously used (or use the defaults)
  wxConfig config(_T("OpenBabelGUI"));
  wxSize size;
  size.SetWidth(config.Read(_T("Width"),1020));
  size.SetHeight(config.Read(_T("Height"),918));
  wxPoint position;
  position.x = config.Read(_T("Left"),2);
  position.y = config.Read(_T("Top"),2);

  wxFileName help, spl;
  wxString pHelp;
  wxGetEnv(_T("BABEL_DATADIR"), &pHelp);
  help.MakeAbsolute(pHelp);
  help.RemoveLastDir();
  help.AppendDir(_T("doc"));
  spl = help;
  help.SetFullName(_T("OpenBabelGUI.html"));
  HelpFile = help.GetFullPath();

  // create the main application window
  OBGUIFrame *frame = new OBGUIFrame(_T("OpenBabelGUI"), position, size);

  // and show it (the frames, unlike simple controls, are not shown when
  // created initially)
  frame->Show(true);
  frame->SetInitialFocus(); //doesn't work!

  //Show the splash screen
  wxImage::AddHandler(new wxPNGHandler);
  wxBitmap bitmap;
  std::ifstream fs; //not actually used
  //Get splash screen from normal OB data directory
  std::string splfile = OpenDatafile(fs, "splash.png");
  if(!splfile.empty() && bitmap.LoadFile(wxString(splfile.c_str(), wxConvUTF8), wxBITMAP_TYPE_PNG))
  {
    wxSplashScreen* splash = new wxSplashScreen(bitmap,
      wxSPLASH_CENTRE_ON_PARENT|wxSPLASH_TIMEOUT,
      4000, NULL, -1, wxDefaultPosition, wxDefaultSize,
      wxBORDER_SIMPLE|wxSTAY_ON_TOP);
  }
  wxImage::RemoveHandler(wxT("PNG file"));
  wxYield();
  return true;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

OBGUIFrame::OBGUIFrame(const wxString& title, wxPoint position, wxSize size)
       : wxFrame(NULL,wxID_ANY,title, position,size), m_pGenOptsPanel(NULL)
{
  // set the frame icon
  SetIcon(wxICON(sample));
  SetDropTarget(new DnD(this)); //newer method for file drap and drop
  wxConfig config(_T("OpenBabelGUI"));

  // create a menu bar
  fileMenu = new wxMenu;
  viewMenu = new wxMenu;
  listMenu = new wxMenu;
  helpMenu = new wxMenu;
  viewMenu->AppendCheckItem(ID_RESTRICTFORMATS, _T("Use &restricted set of formats"));
  viewMenu->Append(ID_SELFORMATS, _T("&Select set of formats"));
  viewMenu->AppendSeparator();
  viewMenu->Append(ID_SETDISPLAYFILE, _T("Configure stucture display"));
  viewMenu->AppendSeparator();
  viewMenu->AppendCheckItem(ID_SHOWAPIOPTIONS, _T("&API Options"),_T("e.g. errorlevel"));
  viewMenu->AppendCheckItem(ID_SHOWCONVOPTIONS, _T("&Conversion Options"));
  viewMenu->AppendCheckItem(ID_SHOWOBJOPTIONS1, _T("Chemical Object Options(&single char)"),
    _T("relating to molecules or reactions\ne.g. -h Add hydrogens"));
  viewMenu->AppendCheckItem(ID_SHOWOBJOPTIONS2, _T("Chemical Object Options(&multichar)"),
    _T("e.g. --addtotitle"));
  viewMenu->AppendCheckItem(ID_SHOWINOPTIONS, _T("&Input Format Options"));
  viewMenu->AppendCheckItem(ID_SHOWOUTOPTIONS, _T("&Output Format Options"));
  viewMenu->AppendSeparator();
  viewMenu->AppendCheckItem(ID_INWRAPPED, _T("Wrap input text (on restart)"));
  viewMenu->AppendCheckItem(ID_OUTWRAPPED, _T("Wrap output text (on restart)"));
  viewMenu->AppendSeparator();
  viewMenu->Append(ID_MINSIZE, _T("Reduce window to minimum useful size"));

  fileMenu->Append(ID_CONVERT, _T("&Convert\tAlt-C"));
  fileMenu->Append(wxID_SAVE, _T("&Save Input Text As...\tCtrl+S"),
    _T("Brings up Save dialog"));
  fileMenu->Append(ID_COPYTOINPUT, _T("Copy Output To Input"),
    _T("Copies output text and format to input"));
  fileMenu->Append(ID_SAVECONFIG, _T("Save screen configuration\tAlt-S"));
  fileMenu->Append(wxID_EXIT, _T("E&xit\tAlt-X"), _T("Quit OpenBabelGUI"));

  helpMenu->Append(wxID_HELP, _T("&Instructions"),
    _T("Opens default browser to view help file"));
  helpMenu->Append(wxID_ABOUT, _T("&About..."), _T("Show about dialog"));

  MakePluginsMenu();

  bool chk;
  config.Read(_T("ShowConvOptions"),&chk,true);
  viewMenu->Check(ID_SHOWCONVOPTIONS,chk);
  config.Read(_T("ShowAPIOptions"),&chk,false);
  viewMenu->Check(ID_SHOWAPIOPTIONS,chk);
  config.Read(_T("ShowObjOptions1"),&chk,true);
  viewMenu->Check(ID_SHOWOBJOPTIONS1,chk);
  config.Read(_T("ShowObjOptions2"),&chk,true);
  viewMenu->Check(ID_SHOWOBJOPTIONS2,chk);
  config.Read(_T("ShowInOptions"),&chk,true);
  viewMenu->Check(ID_SHOWINOPTIONS,chk);
  config.Read(_T("ShowOutOptions"),&chk,true);
  viewMenu->Check(ID_SHOWOUTOPTIONS,chk);
  config.Read(_T("InWrapped"),&chk,true);
  viewMenu->Check(ID_INWRAPPED,chk);
  config.Read(_T("OutWrapped"),&chk,true);
  viewMenu->Check(ID_OUTWRAPPED,chk);
  m_DisplayFile = config.Read(_T("DisplayFile"), wxFileName::GetTempDir()+_T("/gui.svg"));
  m_DisplayCmd  = config.Read(_T("DisplayCmd"), _T("firefox file:///"));

  chk = m_ActiveFormats.ReadConfig(config);
  viewMenu->Check(ID_RESTRICTFORMATS,chk);

  // now append the freshly created menu to the menu bar...
  wxMenuBar *menuBar = new wxMenuBar();
  menuBar->Append(fileMenu, _T("&File"));
  menuBar->Append(viewMenu, _T("&View"));
  menuBar->Append(listMenu, _T("&Plugins"));
  menuBar->Append(helpMenu, _T("&Help"));

  // ... and attach this menu bar to the frame
  SetMenuBar(menuBar);

//******************************************************
//**************** Controls (in tab order)**************

  wxPanel* panel = new wxPanel(this, wxID_ANY);
  m_pOptsWindow =
    new wxScrolledWindow(panel, wxID_ANY,wxDefaultPosition,wxDefaultSize,wxVSCROLL | wxNO_BORDER);
  m_pOptsWindow->SetScrollRate(0,10);

  m_pInFormat    = new wxChoice(panel,ID_INFORMAT,	wxDefaultPosition,wxDefaultSize,
    0, static_cast<wxString*> (NULL));
  m_pInInfo      = new wxButton  (panel, ID_ININFO, wxT("?"),
        wxDefaultPosition,wxDefaultSize,wxBU_EXACTFIT);
  m_pForceInFormat  = new wxCheckBox(panel,ID_INFORCEFORMAT,
        wxT("Use this format for all &input files (ignore file extensions)"));
  m_pInPath      = new wxStaticText(panel,wxID_STATIC,wxT(""),
        wxDefaultPosition,wxDefaultSize,wxST_NO_AUTORESIZE );
  m_pInFilename  = new CFilenames(panel, ID_INFILENAME, wxEmptyString,
        wxDefaultPosition,wxSize(150,20));
  m_pInFiles     = new wxButton  (panel, ID_INGETFILES, wxT("..."),
        wxDefaultPosition,wxSize(35,20));
  m_pInputHere   = new wxCheckBox(panel,ID_INPUTHERE,
        wxT("Input below (ignore input file)"));
  long notwrapped = viewMenu->IsChecked(ID_INWRAPPED) ? 0 : wxTE_DONTWRAP;
  m_pInText      = new wxTextCtrl( panel, ID_INTEXT, _T(""),
        wxDefaultPosition, wxSize(195,200), wxTE_MULTILINE|wxTE_READONLY|notwrapped);

  m_pConvert     = new wxButton  (panel, ID_CONVERT, wxT("&CONVERT"),
        wxDefaultPosition,wxSize(80,40));
  MakeBold(m_pConvert);
  m_pConvert->SetToolTip(_T("Do conversion (Alt C)"));

  m_pOutFormat   = new wxChoice(panel,ID_OUTFORMAT,wxDefaultPosition,wxDefaultSize,
    0, static_cast<wxString*> (NULL));
  m_pOutInfo     = new wxButton  (panel, ID_OUTINFO, wxT("?"),
        wxDefaultPosition,wxDefaultSize,wxBU_EXACTFIT);
  m_pOutFilename = new wxTextCtrl(panel, ID_OUTFILENAME,wxEmptyString,
        wxDefaultPosition,wxSize(150,20));
  m_pOutFiles    = new wxButton  (panel, ID_OUTGETFILES, wxT("..."),
        wxDefaultPosition,wxSize(35,20));
  m_pNoOutFile   = new wxCheckBox(panel,ID_NOOUTFILE,
        wxT("Output below only (no output file)"));
  m_pDisplay     = new wxCheckBox(panel,ID_DISPLAY, wxT("Display in " + m_DisplayCmd.BeforeFirst(' ')));
  notwrapped = viewMenu->IsChecked(ID_OUTWRAPPED) ? 0 : wxTE_DONTWRAP;

  //Output windows: splitter with messages(clog) and output text(cout)
  m_pSplitter = new wxSplitterWindow(panel,wxID_ANY,wxDefaultPosition,wxSize(1955,200));
  m_pMessages = new wxTextCtrl(m_pSplitter,ID_MESSAGES,wxEmptyString,
    wxDefaultPosition,wxDefaultSize,wxTE_MULTILINE|wxTE_READONLY);
  m_pMessages->SetToolTip(_T("Message window. Drag the divider down to make bigger"));
  m_pfixedFont = new wxFont(
    8, //int pointSize
    wxFONTFAMILY_MODERN, //wxFontFamily family
    wxFONTSTYLE_NORMAL, //style
    wxFONTWEIGHT_NORMAL); //wxFontWeight weight
  m_pMessages->SetFont(*m_pfixedFont);

  m_pOutText = new wxTextCtrl(m_pSplitter, ID_OUTTEXT, _T(""),
        wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY|notwrapped);

  m_pSplitter->SplitHorizontally(m_pMessages, m_pOutText, 20);
  m_pSplitter->SetMinimumPaneSize(20);
  int messize;
  config.Read(_T("MessageWindowSize"),&messize,40);
  m_pSplitter->SetSashPosition(messize);

//******************************************************
//************** Layout with Sizers ********************

  topSizer		 = new wxBoxSizer(wxHORIZONTAL);
  InSizer      = new wxBoxSizer(wxVERTICAL);
  OutSizer     = new wxBoxSizer(wxVERTICAL);
  CenterSizer  = new wxBoxSizer(wxVERTICAL);
  OptionsSizer = new wxBoxSizer(wxVERTICAL);
  wxBoxSizer *InFilesSizer = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *OutFilesSizer = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *OutAuxSizer    = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *InFormatSizer = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *OutFormatSizer = new wxBoxSizer(wxHORIZONTAL);

  InFormatSizer->Add(m_pInFormat,1,wxEXPAND);
  InFormatSizer->Add(m_pInInfo,0,wxLEFT,5);

  InFilesSizer->Add(m_pInFilename,1,wxEXPAND);
  InFilesSizer->Add(m_pInFiles,0,wxLEFT,5);

  OutFilesSizer->Add(m_pOutFilename,1,wxEXPAND);
  OutFilesSizer->Add(m_pOutFiles,0,wxLEFT,5);

  OutAuxSizer->Add(m_pNoOutFile,0, wxLEFT|wxBOTTOM,5);
#ifndef __WXMAC__
  if(OBPlugin::GetPlugin(NULL, "0xout")) //display checkbox only if extra output capability is present
    OutAuxSizer->Add(m_pDisplay,0,wxLEFT|wxBOTTOM,5);
#endif

  OutFormatSizer->Add(m_pOutFormat,1,wxEXPAND);
  OutFormatSizer->Add(m_pOutInfo,0,wxLEFT,5);
  wxStaticText* pStatic = new wxStaticText(panel,wxID_STATIC,wxT("        ---- INPUT FORMAT ----"));
  MakeBold(pStatic);
  InSizer->Add(pStatic);
  InSizer->Add(InFormatSizer,0, wxLEFT|wxTOP,5);
  InSizer->Add(m_pForceInFormat,0,wxLEFT|wxBOTTOM,5);
  InSizer->Add(m_pInPath,0,wxEXPAND|wxLEFT|wxTOP,5);
  InSizer->Add(InFilesSizer,0,wxEXPAND|wxALL,5);
  InSizer->Add(m_pInputHere,0,wxLEFT|wxBOTTOM,5);
  InSizer->Add(m_pInText,
      1,        // make vertically stretchable
      wxEXPAND| // make horizontally stretchable
      wxALL,    // and make border all around
      5 );     // set border width

  m_pGenOptsPanel = new DynOptionswx(m_pOptsWindow, OptionsSizer);
  m_pAPIOptsPanel = new DynOptionswx(m_pOptsWindow, OptionsSizer);
  m_pConvOptsPanel = new DynOptionswx(m_pOptsWindow, OptionsSizer);
  m_pInOptsPanel = new DynOptionswx(m_pOptsWindow, OptionsSizer);
  m_pOutOptsPanel = new DynOptionswx(m_pOptsWindow, OptionsSizer);

  m_pOptsWindow->SetSizer(OptionsSizer);

  pStatic = new wxStaticText(panel,wxID_STATIC,wxT("        ---- OUTPUT FORMAT ----"));
  MakeBold(pStatic);
  OutSizer->Add(pStatic);
  OutSizer->Add(OutFormatSizer,0, wxLEFT|wxTOP,5);

  //gap where checkbox used to be
  int chkwidth,chkheight;
  m_pForceInFormat->GetSize(&chkwidth, &chkheight);
  OutSizer->AddSpacer(chkheight+3);

  OutSizer->Add(new wxStaticText(panel,wxID_STATIC,wxT("Output file")),0,wxLEFT|wxTOP,5);
  OutSizer->Add(OutFilesSizer,0,wxEXPAND|wxALL,5);
  OutSizer->Add(OutAuxSizer,0, wxLEFT|wxBOTTOM,5);
  OutSizer->Add(m_pSplitter, 1, wxEXPAND | wxALL, 5 );

  CenterSizer->Add(m_pConvert, 0, wxALL|wxALIGN_CENTER_HORIZONTAL,10);
  CenterSizer->Add(new wxStaticLine(panel),0,wxBOTTOM|wxEXPAND,5);
  CenterSizer->Add(m_pOptsWindow,1);

  topSizer->Add(InSizer,1,wxEXPAND);
  topSizer->Add(CenterSizer,0,wxEXPAND);
  topSizer->Add(OutSizer,1,wxEXPAND);

  panel->SetSizer( topSizer );     // use the sizer for layout
  topSizer->Fit( panel );          // fit the dialog to the contents
//  topSizer->SetSizeHints( panel ); // set hints to honor min size

  GetAvailableFormats();
  wxString inExt = config.Read(_T("InExt"),_T("smi"));
  wxString outExt = config.Read(_T("OutExt"),_T("smi"));
  SetChoice(m_pInFormat, inExt);
  SetChoice(m_pOutFormat,outExt);

  m_InFileBasePath = config.Read(_T("InputPath"), wxGetCwd());
  m_pInPath->SetLabel(ShortenedPath(m_InFileBasePath, *m_pInPath));
  ::wxSetWorkingDirectory(m_InFileBasePath);

  config.Read(_T("InputHere"),&chk,false);
  m_pInputHere->SetValue(chk);
  ChangeInputHere(chk);

  config.Read(_T("ForceInFormat"),&chk,false);
  m_pForceInFormat->SetValue(chk);

 //Display the options
  wxCommandEvent dum;
  OnChangeFormat(dum);
}

//******************************************************
//********** Event handlers ****************************

void OBGUIFrame::OnClose(wxCloseEvent& event)
{
  //Save the window size, in and out formats, and various options, for use next time.
  wxCommandEvent ev;
  OnSaveConfig(ev);

  delete m_pGenOptsPanel;
  delete m_pAPIOptsPanel;
  delete m_pConvOptsPanel;
  delete m_pInOptsPanel;
  delete m_pOutOptsPanel;
  delete m_pfixedFont;
  this->Destroy();
}

void OBGUIFrame::OnSaveConfig(wxCommandEvent& event)
{
    //Save the window size, in and out formats, and various options, for use next time.
  wxConfig config(_T("OpenBabelGUI"));

  int width, height, left, top;
  GetPosition(&left,&top);
  GetSize(&width, &height);
  config.Write(_T("Left"),left);
  config.Write(_T("Top"),top);
  config.Write(_T("Width"),width);
  config.Write(_T("Height"),height);
  config.Write(_T("MessageWindowSize"),m_pSplitter->GetSashPosition());

  config.Write(_T("InputPath"), m_InFileBasePath);

  wxString ext = m_pInFormat->GetStringSelection();
  int pos = ext.find_first_of(_T(" \t-"));
  config.Write(_T("InExt"), ext.substr(0,pos));

  ext = m_pOutFormat->GetStringSelection();
  pos = ext.find_first_of(_T(" \t-"));
  config.Write(_T("OutExt"), ext.substr(0,pos));
  config.Write(_T("DisplayFile"), m_DisplayFile);
  config.Write(_T("DisplayCmd"), m_DisplayCmd);

  config.Write(_T("ShowConvOptions"),viewMenu->IsChecked(ID_SHOWCONVOPTIONS));
  config.Write(_T("ShowAPIOptions"),viewMenu->IsChecked(ID_SHOWAPIOPTIONS));
  config.Write(_T("ShowObjOptions1"),viewMenu->IsChecked(ID_SHOWOBJOPTIONS1));
  config.Write(_T("ShowObjOptions2"),viewMenu->IsChecked(ID_SHOWOBJOPTIONS2));
  config.Write(_T("ShowInOptions"),viewMenu->IsChecked(ID_SHOWINOPTIONS));
  config.Write(_T("ShowOutOptions"),viewMenu->IsChecked(ID_SHOWOUTOPTIONS));
  config.Write(_T("InWrapped"),viewMenu->IsChecked(ID_INWRAPPED));
  config.Write(_T("OutWrapped"),viewMenu->IsChecked(ID_OUTWRAPPED));
  config.Write(_T("InputHere"), m_pInputHere->IsChecked());
  config.Write(_T("ForceInFormat"), m_pForceInFormat->IsChecked());

  m_ActiveFormats.WriteConfig(config);
  config.Write(_T("UseRestrictedFormats"), viewMenu->IsChecked(ID_RESTRICTFORMATS));
}

void OBGUIFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
  // true is to force the frame to close
   Close(true);
}

void OBGUIFrame::OnSaveInputText(wxCommandEvent& WXUNUSED(event))
{
  wxString filter = GetFilter(m_pInFormat) +_T("All files(*.*)|*.*||");
  wxFileDialog dialog(this, _T("Save input text"), m_InFileBasePath,
    _T("FromOB"), filter, wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
  if(dialog.ShowModal() == wxID_OK)
  {
    wxString filename = dialog.GetPath();
    if(!filename.empty())
      m_pInText->SaveFile(filename);
  }
}

void OBGUIFrame::OnCopyToInput(wxCommandEvent& WXUNUSED(event))
{
  //Copies contents of output textbox to input textbox,
  // and output format to input format.
  m_pInText->Clear();
  m_pInText->WriteText(m_pOutText->GetValue());
  SetChoice(m_pInFormat, m_pOutFormat->GetStringSelection());
  m_pInText->SetInsertionPoint(0);
  m_pInputHere->SetValue(true);
  ChangeInputHere(true);
}

void OBGUIFrame::OnHelp(wxCommandEvent& WXUNUSED(event))
{
//  // 1) Get File Type
  wxMimeTypesManager mimeType;
    wxFileType *c_type;
   c_type = mimeType.GetFileTypeFromExtension(_T("html"));
    if(!c_type) return ; //Couldn't find association

    // 2) Get Open Message
    wxString command;

      command = c_type->GetOpenCommand(((OBGUIApp*)wxTheApp)->HelpFile);
    if(!command) return; //No default program

    // 3) Execute message
    wxExecute(command);

    delete c_type;
//	::wxLaunchDefaultBrowser(_T("OpenBabelGUIHelp.html"));
  //wxExecute(_T("OpenBabelGUI.html"));
}
void OBGUIFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
  //(wxT not required after v2.90)  
  std::string cmsg =
    "OpenBabelGUI (C) 2006 by Chris Morley\n\n"
    "This program is part of the OpenBabel project,\n"
    "which is released under the GNU General Public License.\n"
    "See: http://openbabel.org/wiki/Main_Page\n\n"
    "For more detailed information see: "
    "http://openbabel.org/wiki/Windows_GUI\n\n"
    "Please cite the paper:\n"
    "\"Open Babel: An open chemical toolbox\"\n"
    "N O'Boyle, M Banck, C A James, C Morley,\n"
    "T Vandermeersch and G R Hutchison\n"
    "Journal of Cheminformatics 2011, 3:33\n"
    "doi:10.1186/1758-2946-3-33\n\n"
    "OpenBabel version ";

  wxString msg(cmsg.c_str(), wxConvUTF8);
  msg << _T(BABEL_VERSION);
  wxMessageBox(msg, _T("About OpenBabelGUI"), wxOK | wxICON_INFORMATION | wxCENTER, this);
}
///////////////////////////////////////////

void OBGUIFrame::OnSelectFormats(wxCommandEvent& event)
{
  if(m_ActiveFormats.SelectFormats())
  {
    //When formats have been selected always set the Use Restricted Set check
    viewMenu->Check(ID_RESTRICTFORMATS,true);
    GetAvailableFormats();
  }
}
void OBGUIFrame::OnRestrictFormats(wxCommandEvent& event)
{
  GetAvailableFormats();
}
void OBGUIFrame::OnSetDisplayFile(wxCommandEvent& event)
{
  wxTextEntryDialog dialog(this, _T("Enter display command and temporary display file on separate lines"),
    _T("Parameters for structure display"), m_DisplayCmd + _T('\n') + m_DisplayFile,
    wxTE_MULTILINE | wxOK | wxCANCEL);
  if (dialog.ShowModal() == wxID_OK)
  {
    m_DisplayCmd  = dialog.GetValue().BeforeFirst('\n');
    m_DisplayFile = dialog.GetValue().AfterFirst('\n');
    m_pDisplay->SetLabel(_T("Display in ") + m_DisplayCmd.BeforeFirst(' '));
  }
}

void OBGUIFrame::OnSetMinSize(wxCommandEvent& event)
{
  SetSize(872, 249);
}

void OBGUIFrame::OnConvert(wxCommandEvent& WXUNUSED(event))
{
  wxBusyCursor cw;
  m_pOutText->Clear();
  m_pMessages->Clear();

  //Default input is from input text box;
  std::stringstream ss(std::string(m_pInText->GetValue().mb_str()));
  //Default output is a string stream which is written to the Output text box at the end
  std::stringstream GUIostream;
  // LNK2005 error with VS10 (but not VS9) when stringstream used, possibly related to:
  // http://osdir.com/ml/OpenSceneGraph-Users/2010-08/msg00501.html
  // An unsatisfactory workaround is to set in obgui project
  // ConfigurationProperties/Linker/General/Force File Output to
  // Multiply Defined Symbol Only (/FORCE:MULTIPLE)

  OBConversion Conv(&ss, &GUIostream);

  int iSel = m_pInFormat->GetSelection();
  if((iSel)<0) return;
  int oSel = m_pOutFormat->GetSelection();
  if((oSel)<0) return;

  OBFormat* pInFormat = NULL;
  OBFormat* pOutFormat = (OBFormat*)m_pOutFormat->GetClientData(oSel);
  if(m_pForceInFormat->IsChecked() || m_pInputHere->IsChecked())
    pInFormat = (OBFormat*)m_pInFormat->GetClientData(iSel);
  Conv.SetInAndOutFormats( pInFormat,pOutFormat);

  DoOptions(Conv);

  //Setup output file
  std::string stdOutputFileName;
  wxString OutFileName = m_pOutFilename->GetValue();
  //If you are trying to output with no filename what you really wanted
  //was to output to the OutConsole
  if(m_pNoOutFile->IsChecked() || OutFileName.IsEmpty())
    m_pNoOutFile->SetValue(true);
  else
  {
    //If output filename has no path, use input path
    wxFileName fn(OutFileName);
    if(!fn.IsAbsolute())
    {
      fn.MakeAbsolute(m_InFileBasePath);
      OutFileName = fn.GetFullPath();
    }
    stdOutputFileName = OutFileName.mb_str();
  }
  OBFormat* pOutFileFormat = Conv.FormatFromExt(OutFileName.mb_str());
  if(!m_pNoOutFile->IsChecked() && pOutFileFormat && (pOutFileFormat!=pOutFormat))
    if(wxMessageBox(_T("The output file name extension does not correspond \
with the output format.\nDo you wish to continue the conversion?"),
      _T("Is the output filename correct?"), wxOK | wxCANCEL)!=wxOK)
      return;

  //Setup input file
  std::vector<std::string> FileList, OutputFileList;
  if(!m_pInputHere->IsChecked())
    m_pInFilename->Expand(FileList);

  //redirect cerr & clog & cout
#ifndef __WXMAC__
    wxStreamToTextRedirector cerrCapture(m_pMessages, &std::cerr);
    wxStreamToTextRedirector clogCapture(m_pMessages, &std::clog);
#endif
//		wxStreamToTextRedirector coutCapture(m_pOutText);

//	m_pOutText->Freeze();//Otherwise seems to be redrawn after each char from cout

  //2D depiction in svg (or other) format automatically sent to file
  if(m_pDisplay->IsChecked() && !m_DisplayFile.empty())
    Conv.AddOption("0xout", OBConversion::GENOPTIONS, m_DisplayFile.mb_str());

  int count = Conv.FullConvert(FileList, stdOutputFileName, OutputFileList);

//	m_pOutText->Thaw();

  Conv.ReportNumberConverted(count);

  if(OutputFileList.size()>1)
  {
    std::clog << '\n' << OutputFileList.size()
      << " files output. The first is " << OutputFileList[0];
  }

  if(count>0)
  {
    if(!m_pNoOutFile->IsChecked())
    {
      //Read back file and add to output console
      m_pOutText->Clear();
      m_pOutText->LoadFile(wxString(OutputFileList[0].c_str(), wxConvUTF8));
      //m_pOutText->LoadFile(_T(OutputFileList));
    }
    else
    {
      m_pOutText->WriteText(wxString(GUIostream.str().c_str(), wxConvUTF8));
      m_pOutText->SetInsertionPoint(0);
    }
  }
#ifndef __WXMAC__
  //Call Firefox to display the 2D structure
  if(m_pDisplay->IsChecked() && wxFile::Exists(m_DisplayFile.Trim()))
  {
    wxExecute(m_DisplayCmd + m_DisplayFile);
  }
#endif
}

///////////////////////////////////////////
void OBGUIFrame::OnInFormatInfo(wxCommandEvent& WXUNUSED(event))
{
  int nSel=m_pInFormat->GetSelection();
  if(nSel<0) return;
  OBFormat* pFormat = (OBFormat*)m_pInFormat->GetClientData(nSel);
  wxString mes(pFormat->Description(), wxConvUTF8);
  wxString url(pFormat->SpecificationURL(), wxConvUTF8);
  if(!url.IsEmpty())
    mes += _T("\nURL for specification: ") + url;
  wxMessageBox(mes, _T("Format info"));
}

///////////////////////////////////////////
void OBGUIFrame::OnOutFormatInfo(wxCommandEvent& WXUNUSED(event))
{
  int nSel=m_pOutFormat->GetSelection();
  if(nSel<0) return;
  OBFormat* pFormat = (OBFormat*)m_pOutFormat->GetClientData(nSel);
  wxString mes(pFormat->Description(), wxConvUTF8);
  wxString url(pFormat->SpecificationURL(), wxConvUTF8);
  if(!url.IsEmpty())
    mes += _T("\nURL for specification: ") + url;
  wxMessageBox(mes, _T("Format info"));
}

void OBGUIFrame::OnGetInputFile(wxCommandEvent& WXUNUSED(event))
{
  wxFileDialog dialog(this,_T("Choose Input File"),m_InFileBasePath,_T(""),
      GetFilter(m_pInFormat) + InputFilterString,
      wxFD_MULTIPLE | wxFD_FILE_MUST_EXIST );
  if(dialog.ShowModal() == wxID_OK)
  {
//		m_pInFilename->Clear();
    wxArrayString filepatharray;
    dialog.GetPaths(filepatharray);
    DisplayInputFiles(filepatharray);
  }
}
void OBGUIFrame::DisplayInputFiles(wxArrayString filepatharray)
{
  int i, endsel=0, startsel=0;
  if(!wxGetKeyState(WXK_CONTROL))
  {
    m_pInFilename->Clear();
    wxFileName filenamewx(filepatharray[0]);
    m_InFileBasePath = filenamewx.GetVolume() + filenamewx.GetVolumeSeparator()
      + filenamewx.GetPath(wxPATH_GET_SEPARATOR );
    m_pInPath->SetLabel(ShortenedPath(m_InFileBasePath, *m_pInPath));
    ::wxSetWorkingDirectory(m_InFileBasePath);
  }
  else
  {
    if(!m_pInFilename->GetValue().empty())
      m_pInFilename->AppendText(_T(';'));
    startsel  = m_pInFilename->GetLastPosition();
  }

  for(i=0;i<filepatharray.GetCount();++i)
  {
    if(i==1)
      endsel = m_pInFilename->GetLastPosition();
    if(i!=0)
      m_pInFilename->AppendText(_T(';'));
    wxFileName fnamewx(filepatharray[filepatharray.GetCount()-1-i]);
    fnamewx.MakeRelativeTo(m_InFileBasePath);
    m_pInFilename->AppendText(fnamewx.GetFullPath());

    if(wxGetKeyState(WXK_CONTROL) && i==0)
      endsel = m_pInFilename->GetLastPosition();
  }
  m_pInFilename->SetSelection(startsel,endsel);

  m_pInText->Clear();
  m_pInText->LoadFile(filepatharray[filepatharray.GetCount()-1]);
  m_pInFilename->SetFocus();
}

void OBGUIFrame::OnGetOutputFile(wxCommandEvent& WXUNUSED(event))
{
  wxFileDialog dialog(this,_T("Choose Output File"),_T(""),_T(""),
      GetFilter(m_pOutFormat) + OutputFilterString,
      wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
  if(dialog.ShowModal() == wxID_OK)
  {
    wxString filepath = dialog.GetPath();
    m_pOutFilename->Clear();
    m_pOutFilename->AppendText(filepath);
  }
}

void OBGUIFrame::OnOutFileNameUpdate(wxUpdateUIEvent& event)
{
  event.Enable(!m_pNoOutFile->IsChecked());
}
void OBGUIFrame::OnInFileNameUpdate(wxUpdateUIEvent& event)
{
  event.Enable(!m_pInputHere->IsChecked());
}

void OBGUIFrame::OnChangeInputHere(wxCommandEvent& event)
{
  ChangeInputHere(event.IsChecked());
}

void OBGUIFrame::ChangeInputHere(bool chk)
{
  wxColour bg = chk ? wxColour(250,255,210) : wxColour(250,255,255);
  m_pInText->SetEditable(chk);
  m_pInText->SetBackgroundColour(bg);
  m_pInText->Refresh();
  m_pInText->SetFocus();
  m_pInFilename->Enable(chk);
}

void OBGUIFrame::OnChangeFormat(wxCommandEvent& WXUNUSED(event))
{
  //Display the options
  m_pAPIOptsPanel->Clear();
  m_pConvOptsPanel->Clear();
  m_pGenOptsPanel->Clear();
  m_pInOptsPanel->Clear();
  m_pOutOptsPanel->Clear();

  if(viewMenu->IsChecked(ID_SHOWAPIOPTIONS))
  {
    OBFormat* pAPI= OBConversion::FindFormat("obapi");
    if(pAPI)
      m_pAPIOptsPanel->Construct(pAPI->Description());
  }

  if(viewMenu->IsChecked(ID_SHOWCONVOPTIONS))
    m_pConvOptsPanel->Construct(OBConversion::Description());

  OBFormat* pInFormat = (OBFormat*)m_pInFormat->GetClientData(m_pInFormat->GetSelection());
  OBFormat* pOutFormat = (OBFormat*)m_pOutFormat->GetClientData(m_pOutFormat->GetSelection());
  if(!pInFormat || !pOutFormat)
    return;
  if(viewMenu->IsChecked(ID_SHOWOBJOPTIONS1)) //Single char
    m_pGenOptsPanel->Construct(pOutFormat->TargetClassDescription(),NULL,1);
  if(viewMenu->IsChecked(ID_SHOWOBJOPTIONS2)) //Multi char
    m_pGenOptsPanel->Construct(pOutFormat->TargetClassDescription(),NULL,2);

  if(viewMenu->IsChecked(ID_SHOWINOPTIONS))
  {
    if(pInFormat && !m_pInOptsPanel->Construct(pInFormat->Description(),"input"))
      m_pInOptsPanel->Construct(pInFormat->Description(),"read");//try again
  }
  if(viewMenu->IsChecked(ID_SHOWOUTOPTIONS))
  {
    if(pOutFormat && !m_pOutOptsPanel->Construct(pOutFormat->Description(),"output"))
      m_pOutOptsPanel->Construct(pOutFormat->Description(), "write");//try again
  }

  CenterSizer->Fit(m_pOptsWindow);
  CenterSizer->Layout();
  topSizer->Layout();
}

void OBGUIFrame::SetInitialFocus()
{
  m_pInFilename->SetFocus();
}

//*******************************************************
//********** Local functions ****************************

/**	 Set the selection of the control to the string which starts with
     the file extension (case independent). The extension .gz is ignored.
     The string parameter can be the extension itself.
     Returns false if no match found.
 **/
bool OBGUIFrame::SetChoice(wxChoice* pChoice, const wxString& FileName)
{
  wxString ext(FileName);
  int pos = FileName.rfind('.');
  if(pos!=-1)
  {
    wxString ext = FileName.substr(pos+1);
    if(FileName.substr(pos)==_T(".gz"))
      pos = FileName.rfind('.',pos-1);
  }
  for(int iSel=0;iSel<pChoice->GetCount();iSel++)
  {
    wxString txt(pChoice->GetString(iSel).substr(0,ext.length()));
    if(txt.MakeUpper()==ext.MakeUpper())
    {
      pChoice->SetSelection(iSel);
      return true;
    }
  }
  return false;
}
///////////////////////////////////////////////////
wxString OBGUIFrame::GetFilter(wxChoice* pChoice)
{
  //Uses text from the window (input or output combo) to construct filter string
  wxString  txt = pChoice->GetStringSelection();
  int pos1 = txt.find(_T(" "));
  int pos2 = txt.find_first_not_of( _T(" :-\t"), pos1);
  return txt.substr(pos2) + _T(" (*.") + txt.substr(0, pos1) + _T(")|*.")
    + txt.substr(0, pos1) + _T(";*.") + txt.substr(0, pos1) + _T(".gz|");
}
/////////////////////////////////////////////////////
wxString OBGUIFrame::ShortenedPath(const wxString& path, const wxWindow& wnd, int wndwidth)
{
  /* Tries to leave first and last folder names, so that with at least
  four separators for the path to be shortened like:

  c:\My Documents\projects\wxWindows\test\readme.txt
  to
  c:\My Documents\...\wxWindows\test\readme.txt
  or
  c:\My Documents\...\test\readme.txt

  If the width is still too large or there are 3 or fewer separators,
  the path is returned unshortened.:
  */
  int txtwidth, txtheight, wndheight;
  wnd.GetTextExtent(path, &txtwidth, &txtheight);
  if(wndwidth<0)
    wnd.GetClientSize(&wndwidth,&wndheight);
  if(txtwidth <= wndwidth) return path;

  // Set pos1 at the second separator or return unchanged string if not possible
  size_t pos1 = path.find_first_of(_T("/\\"));
  if(pos1==wxNOT_FOUND || pos1+1 == path.length()) return path;
  pos1= path.find_first_of(_T("/\\"), pos1+1);
  if(pos1==wxNOT_FOUND || pos1+1 == path.length()) return path;

  wxString tpath(path);
  size_t pos2 = pos1;
  while(txtwidth > wndwidth)
  {
    pos2= tpath.find_first_of(_T("/\\"), pos2+1);
    if(pos2==wxNOT_FOUND || pos2+1 == tpath.length()) return path;
    //pos2 now has next separator
    //Ensure that there is at least one further directory
    if(tpath.find_first_of(_T("/\\"), pos2+1)==wxNOT_FOUND) return path;

    tpath.replace(pos1+1, pos2-pos1-1,_T("..."));
    pos2=pos1+4;
    wnd.GetTextExtent(tpath, &txtwidth, &txtheight);
  }
  return tpath;
}

void OBGUIFrame::DoOptions(OpenBabel::OBConversion& Conv)
{
  // Is a API directive, e.g.---errorlevel
  //Send to the pseudoformat "obapi" (without any leading -)
  OBFormat* pAPI= OpenBabel::OBConversion::FindFormat("obapi");
  if(pAPI)
  {
    OBConversion apiConv;
    if(m_pAPIOptsPanel->SetOptions(apiConv, OBConversion::GENOPTIONS))
    {
      apiConv.SetOutFormat(pAPI);
      apiConv.Write(NULL, &std::cout);
    }
  }
  m_pGenOptsPanel->SetOptions(Conv, OBConversion::GENOPTIONS);
  m_pConvOptsPanel->SetOptions(Conv, OBConversion::GENOPTIONS);
  m_pInOptsPanel->SetOptions(Conv, OBConversion::INOPTIONS);
  m_pOutOptsPanel->SetOptions(Conv, OBConversion::OUTOPTIONS);
}

void OBGUIFrame::GetAvailableFormats()
{
  //Get data on available formats and add to comboboxes and to filter string
  //OBConversion dummy; //needed for OBConversion to load format classes
  int nInSel=0,nOutSel=0;
  m_pInFormat->Clear();
  m_pOutFormat->Clear();
  InputFilterString=_T("All Chemical Formats|*.");
  OutputFilterString = InputFilterString;
  m_ActiveFormats.Clear();

  std::vector<std::string> vec;
  if(OBPlugin::ListAsVector("formats", NULL, vec))//check that there are some formats
  {
    OBFormat::PluginIterator itr;
    for(itr=OBFormat::Begin("formats");itr!=OBFormat::End("formats");++itr)
    {
      OBFormat* pFormat = static_cast<OBFormat*>(itr->second);
      if((pFormat->Flags() & NOTWRITABLE) && (pFormat->Flags() & NOTREADABLE))
        continue;

      std::string stxt;
      itr->second->Display(stxt, NULL, itr->first);
      wxString txt(stxt.c_str(), wxConvUTF8);
      int pos = txt.find('[');
      if(pos!=wxString::npos)
        txt.erase(pos);
      //Check whether is in restricted set of formats, if this is in use
      if(!m_ActiveFormats.Add(txt) && viewMenu->IsChecked(ID_RESTRICTFORMATS))
        continue;
      int n;
      if(!(pFormat->Flags() & NOTREADABLE))
      {
        n = m_pInFormat->Append(txt,pFormat);
        InputFilterString+=txt.Left(txt.Find(_T(" ")));
        InputFilterString+=_T(";*.");
      }
      if(!(pFormat->Flags() & NOTWRITABLE))
      {
        n = m_pOutFormat->Append(txt,pFormat);
        OutputFilterString+=txt.Left(txt.Find(_T(" ")));
        OutputFilterString+=_T(";*.");
      }
    }
  }
  if(m_pInFormat->GetCount()==0)
  {
    m_pInFormat->Append(_T("No input format in the selected set"));
    m_pInFormat->SetClientData(0,NULL);
  }
  if(m_pOutFormat->GetCount()==0)
  {
    m_pOutFormat->Append(_T("No output format in the selected set"));
    m_pOutFormat->SetClientData(0,NULL);
  }
  m_pInFormat->SetSelection(nInSel);
  m_pOutFormat->SetSelection(nOutSel);

  InputFilterString = InputFilterString.Left(InputFilterString.Length()-3); //remove unneeded ;*.
  OutputFilterString = OutputFilterString.Left(OutputFilterString.Length()-3); //remove unneeded ;*.
  InputFilterString+=_T("|AllFiles(*.*)|*.*||");
  OutputFilterString+=_T("|AllFiles(*.*)|*.*||");
}

void OBGUIFrame::DisplayInFile(wxString filename)
{
  m_pInText->Clear();
  wxFileName fn(filename);
  fn.MakeAbsolute(m_InFileBasePath);
  m_pInText->LoadFile(fn.GetFullPath());
}

void OBGUIFrame::MakeBold(wxWindow* pWnd)
{
  wxFont font = pWnd->GetFont();
  font.SetWeight(wxFONTWEIGHT_BOLD );
  pWnd->SetFont(font);
}

//**********************************************
BEGIN_EVENT_TABLE(CFilenames,wxTextCtrl)
  EVT_LEFT_DCLICK(CFilenames::OnDblClick)
  EVT_CHAR(CFilenames::OnKeyPress)
END_EVENT_TABLE()

void CFilenames::OnDblClick(wxMouseEvent& event)
{
  //extract double-clicked filename
  OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
  wxASSERT(frame);
  frame->DisplayInFile(SelectFilename());
}
/// Highlights and returns the filename at the current cursor position
wxString CFilenames::SelectFilename()
{
  wxString fname(GetValue());
  int endsel;
  long n = GetInsertionPoint();
  int pos = fname.find(';',n);
  endsel=pos;
  if(pos!=-1)
    fname.erase(pos);
  else
    endsel=GetLastPosition();
  pos = fname.rfind(';',n);
  if(pos!=-1)
    fname.erase(0,pos+1);
  SetSelection(pos+1,endsel);

  fname.Trim();
  fname.Trim(false);
  return fname;
}

void OBGUIFrame::OnMouseWheel(wxMouseEvent& event)
{
  int delta = event.GetWheelRotation() / event.GetWheelDelta();
  m_pInFilename->ToNextFile(-delta); //is direction correct?
}


void OBGUIFrame::MakePluginsMenu()
{
  int n=0;
  std::vector<std::string> topvec;
  OBPlugin::ListAsVector(NULL, NULL, topvec);//get 'formats', 'fingerprints', etc
  for(int itop=0; itop < topvec.size(); ++itop)
  {
    wxMenu* subMenu = new wxMenu();
    std::vector<std::string> subvec;
//    std::vector<std::string> verbosevec;
    OBPlugin::ListAsVector(topvec[itop].c_str(), NULL, subvec);//get each format, etc as single line
//    OBPlugin::ListAsVector(topvec[itop].c_str(), "verbose", verbosevec);//get full description of each format
    subMenu->Append(ID_HINT,_T("     (Click to see details and to copy ID to clipboard)"));
    for(int isub=0; isub < subvec.size(); ++isub)
    {
//      wxMenu* subsubMenu = new wxMenu();
//      subsubMenu->Append(ID_PLUGINS+n++,wxString(verbosevec[isub].c_str()));
//      subMenu->AppendSubMenu(subsubMenu, wxString(subvec[isub].c_str()));
      subMenu->Append(ID_PLUGINS+n++,wxString(subvec[isub].c_str(), wxConvUTF8));
    }
    listMenu->AppendSubMenu(subMenu, wxString(topvec[itop].c_str(), wxConvUTF8),
      _T("Plugin Classes"));
  }
}

void OBGUIFrame::OnClickPlugin(wxCommandEvent& event)
{
  //Display message box with info on plugin class
  int nID = event.GetId();
  wxString plugintype;
  wxMenu* parent;
  wxMenuItem* item = listMenu->FindItem(nID, &parent);
  if(item)
  {
    //Find the name of the plugin type. It seems difficult to go up a menu hierarchy.
    wxwxMenuItemListNode *node = listMenu->GetMenuItems().GetFirst();
    while (node)
    {
      wxMenuItem* itemm = node->GetData();
      if(itemm->GetSubMenu()==parent)
      {
        plugintype = itemm->GetItemLabel();
        break;
      }
      node = node->GetNext();
    }

    wxString id = item->GetItemLabel().BeforeFirst(' ');
    OBPlugin* plugin = OBPlugin::GetPlugin(plugintype.mb_str(), id.mb_str());
    if(plugin)
    {
      std::string txt;
      plugin->Display(txt, "verbose", id.mb_str());
      wxMessageBox(wxString(txt.c_str(), wxConvUTF8), _T("Plugin details"), wxOK | wxICON_INFORMATION | wxCENTER, this);
    }
    //And also put ID on clipboard
    if (wxTheClipboard->Open())
    {
      wxString itemText = item->GetItemLabel();
      wxTheClipboard->SetData( new wxTextDataObject(itemText.BeforeFirst(' ')) );
    }
    wxTheClipboard->Close();
  }
}

void CFilenames::OnKeyPress(wxKeyEvent& event)
{
  int delta=1;
  OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
  switch (event.GetKeyCode())
  {
  case 313:
  case WXK_PAGEDOWN: //why is code not correct?
    delta=-1;
  case 312:
  case WXK_PAGEUP:
    ToNextFile(delta);
    break;
  case WXK_TAB:
    if(event.ShiftDown())
      ToNextFile(-1);
    else
      ToNextFile(+1);
    break;

  case WXK_RETURN:
    if(event.ShiftDown())
      SetValue(nameWithWildcard);
    else
    {
      nameWithWildcard = GetValue();
      if(nameWithWildcard.find_first_of(_T("*?"))==wxNOT_FOUND)
      {
        wxCommandEvent dum;
        frame->OnConvert(dum);
      }
      else
      {
        //expand wildcards in each filename
        std::vector<std::string> filelist;
        int count = Expand(filelist);
        wxString mes;
        mes << count << _T(" files found");
        frame->DisplayMessage(mes);

        Clear();
        int i, endsel=0;
        for(i=0;i<filelist.size();++i)
        {
          if(i!=0)
          {
            if(i==1)
              endsel = GetLastPosition();
            AppendText(_T(';'));
          }
          wxString name(filelist[i].c_str(), wxConvUTF8);
          wxFileName fnamewx(name);
          fnamewx.MakeRelativeTo(frame->GetInFileBasePath());
          AppendText(fnamewx.GetFullPath());
        }
        SetSelection(0,endsel);
      }
    }
    break;
  default:
    event.Skip();
  }
}

int  CFilenames::Expand(std::vector<std::string>& filelist)
{
  //Adds full path names of all input files to filelist, expanding wildcards they are present.
  OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
  int namestart=0, count=0;
  wxString txt(GetValue());
  do //for each filename
  {
    int nameend = txt.find(';',namestart);
    wxString name = txt.substr(namestart, nameend-namestart);
    name.Trim().Trim(false);
    if(name.IsEmpty())
      break;
    wxFileName fn(name);
    namestart=nameend+1; // 0 at end
    fn.MakeAbsolute(frame->GetInFileBasePath());
    count += DLHandler::findFiles(filelist, std::string(fn.GetFullPath().mb_str()));
  }while(namestart);

  return count;
}

bool CFilenames::ToNextFile(int delta)
{
  wxString fname(GetValue());
  if (fname.IsEmpty())
    return false;
  long n = GetInsertionPoint();
  int pos;
  if(delta>0)
  {
    pos = fname.find(';',n);
    if(pos==-1)
      pos=GetLastPosition();
    else
      ++pos;
  }
  else
    pos = fname.rfind(';',n);
  if(pos!=-1)
  {
    SetInsertionPoint(pos);
    OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
    wxString nxtname = SelectFilename();
    if(nxtname.IsEmpty())
      return false;
    frame->DisplayInFile(nxtname);
    return true;
  }
  else
    return false;
}
//***********************************************


