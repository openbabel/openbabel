/**********************************************************************
OBGUI.h -  Cross-platform Graphical User Interface

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
#ifndef OB_GUI_H
#define OB_GUI_H

#include "stdwx.h"
#include <wx/dnd.h>
#include "optswx.h"
#include <sstream>

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

// Define a new application type, each program should derive a class from wxApp
class OBGUIApp : public wxApp
{
public:
    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows one to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit();
    wxString HelpFile;
};

//*******************************************
///Class for input filenames textctrl
class CFilenames : public wxTextCtrl
{
public:
  CFilenames(wxWindow* parent, wxWindowID id, const wxString& value = _T(""),
    const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize)
    : wxTextCtrl(parent,id,value,pos,size,wxTE_NOHIDESEL|wxTE_PROCESS_TAB)
  {
    SetToolTip(_T("With multiple files, display any of them below using double click, tab or mousewheel"));
  }
  void OnDblClick(wxMouseEvent& event);
  void OnKeyPress(wxKeyEvent& event);
  bool ToNextFile(int delta);
  wxString SelectFilename();
  int  Expand(std::vector<std::string>& filelist);

private:
  DECLARE_EVENT_TABLE()
  wxString nameWithWildcard;
};


//*******************************************
/// The main window
class OBGUIFrame : public wxFrame
{
public:
  OBGUIFrame(const wxString& title, wxPoint position, wxSize size);

  // event handlers
  void OnQuit(wxCommandEvent& event);
  void OnSaveInputText(wxCommandEvent& event);
  void OnCopyToInput(wxCommandEvent& event);
  void OnAbout(wxCommandEvent& event);
  void OnHelp(wxCommandEvent& event);
  void OnGetInputFile(wxCommandEvent& event);
  void OnGetOutputFile(wxCommandEvent& event);
  void OnInFormatInfo(wxCommandEvent& event);
  void OnOutFormatInfo(wxCommandEvent& event);
  void OnOutFileNameUpdate(wxUpdateUIEvent& event);
  void OnInFileNameUpdate(wxUpdateUIEvent& event);
  void OnConvert(wxCommandEvent& event);
  void OnChangeInputHere(wxCommandEvent& event);
  void ChangeInputHere(bool chk);
  void OnChangeFormat(wxCommandEvent& event);
  void OnClose(wxCloseEvent& event);
  void OnMouseWheel(wxMouseEvent& event);
  void OnSelectFormats(wxCommandEvent& event);
  void OnRestrictFormats(wxCommandEvent& event);
  void OnSetDisplayFile(wxCommandEvent& event);
  void OnClickPlugin(wxCommandEvent& event);
  void OnSetMinSize(wxCommandEvent& event);
  void OnSaveConfig(wxCommandEvent& event);
  void OnShowDataDir(wxCommandEvent& event);

  void DisplayInFile(wxString filename);
  wxString GetInFileBasePath(){ return m_InFileBasePath;}
  void DisplayMessage(wxString& message){m_pMessages->SetValue(message);}
  void DisplayInputFiles(wxArrayString filepatharray);
  void SetInitialFocus();

private:
  // any class wishing to process wxWidgets events must use this macro
  DECLARE_EVENT_TABLE()

  wxMenu* fileMenu;
  wxMenu* listMenu;
  wxMenu* viewMenu;
  wxMenu* helpMenu;

  wxStaticText* m_pInPath;
  wxChoice*   m_pInFormat;
  wxChoice*   m_pOutFormat;
  wxCheckBox* m_pForceInFormat;
  wxCheckBox* m_pNoOutFile;
  wxCheckBox* m_pDisplay;
  wxCheckBox* m_pInputHere;
  wxButton*   m_pInFiles;
  wxButton*   m_pOutFiles;
  CFilenames* m_pInFilename;
  wxTextCtrl* m_pOutFilename;
  wxTextCtrl* m_pInText;
  wxTextCtrl* m_pOutText;
  wxSplitterWindow* m_pSplitter;
  wxButton*   m_pInInfo;
  wxButton*   m_pOutInfo;
  wxButton*   m_pConvert;
  wxTextCtrl* m_pMessages;

  DynOptionswx* m_pGenOptsPanel;
  DynOptionswx* m_pAPIOptsPanel;
  DynOptionswx* m_pConvOptsPanel;
  DynOptionswx* m_pInOptsPanel;
  DynOptionswx* m_pOutOptsPanel;
  wxScrolledWindow* m_pOptsWindow;

  wxBoxSizer *topSizer;
  wxBoxSizer *InSizer;
  wxBoxSizer *OutSizer;
  wxBoxSizer* CenterSizer;
  wxBoxSizer* OptionsSizer;

  wxString InputFilterString, OutputFilterString;
  wxString    m_InFileBasePath;
  wxString    m_DisplayFile, m_DisplayCmd;
  wxFont* m_pfixedFont;
  ActiveFormats m_ActiveFormats;

private:
  void DoOptions(OpenBabel::OBConversion& Conv);

  /// Returns a file path name shortened (using /.../) so that will fit in a specified window
  /// The width in pixels can optionally be explicitly specified.
  wxString ShortenedPath(const wxString& path, const wxWindow& wnd, int wndwidth=-1);

  void GetAvailableFormats();
  wxString GetFilter(wxChoice* pChoice);

  bool SetChoice(wxChoice* pChoice, const wxString& FileName);
  void MakeBold(wxWindow* pWnd);
  void MakePluginsMenu();
};

class MyDialog : public wxDialog
{
public:
  MyDialog(wxWindow *parent, const wxString &title );
};

class DnD : public wxFileDropTarget
{
public:
  DnD(OBGUIFrame* parent) : m_pParent(parent){};
  virtual bool OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
  {
    m_pParent->DisplayInputFiles(filenames);
    return true;
  }
private:
  OBGUIFrame* m_pParent;
};




// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// IDs for the controls and the menu commands
enum
{
    ID_BUTTON = wxID_HIGHEST+1,
    ID_VIEW,
    ID_SHOWCONVOPTIONS,ID_SHOWAPIOPTIONS,ID_SHOWOBJOPTIONS1,ID_SHOWOBJOPTIONS2,
    ID_SHOWINOPTIONS,ID_SHOWOUTOPTIONS,
    ID_INWRAPPED,ID_OUTWRAPPED,
    ID_TEXTCTRL, ID_INFO, ID_INTEXT, ID_OUTTEXT,
    ID_INFILENAME,ID_INGETFILES,
    ID_OUTFILENAME,ID_OUTGETFILES,
    ID_INFORMAT,ID_OUTFORMAT,ID_ININFO,ID_OUTINFO,
    ID_INFORCEFORMAT,ID_OUTFORCEFORMAT,ID_NOOUTFILE,ID_CONVERT,
    ID_MESSAGES,ID_INPUTHERE,ID_RESTRICTFORMATS,ID_SELFORMATS,
    ID_COPYTOINPUT,ID_SAVECONFIG,ID_HINT,ID_SHOWDATADIR,
    ID_PLUGINS,ID_DISPLAY,ID_SETDISPLAYFILE,ID_MINSIZE,ID_PLUGININFO

};

#endif
