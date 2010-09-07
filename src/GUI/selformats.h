#ifndef SELFORMATS_DIALOG_H
#define SELFORMATS_DIALOG_H

// Class for dialog to select active formats
class SelFormatsDialog : public wxDialog
{
private:
	DECLARE_EVENT_TABLE()
	wxCheckListBox* pList;
	wxString& active;
public:
	SelFormatsDialog(wxArrayString& AllFormatsArray, wxString& ActiveFormatsString)
		: wxDialog(NULL, wxID_ANY, _T("Set Active Formats"),
		wxPoint(100,60), wxDefaultSize, wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER),
		active(ActiveFormatsString)
	{
		//Cannot get SetSize() to work here so use fixed item height
		size_t height = 5 + 15 * AllFormatsArray.GetCount();
		size_t maxheight = static_cast<size_t> (0.8*wxSystemSettings::GetMetric(wxSYS_SCREEN_Y));
		if(height > maxheight)
			height = maxheight;
		pList = new wxCheckListBox(this, wxID_ANY, wxDefaultPosition,
			wxSize(-1,height), AllFormatsArray);

		wxBoxSizer* topSizer = new wxBoxSizer(wxHORIZONTAL);
		topSizer->Add(pList,1,wxEXPAND);
		wxBoxSizer* buttonSizer = new wxBoxSizer(wxVERTICAL);
		buttonSizer->Add(new wxButton(this, wxID_OK, wxT("OK")),0,wxALL,10);
		buttonSizer->Add(new wxButton(this,wxID_CANCEL),0,wxALL,10);

		topSizer->Add(buttonSizer);
		SetSizer(topSizer);
		topSizer->Fit(this);
		topSizer->SetSizeHints(this);

		unsigned i;
		for(i=0;i<AllFormatsArray.GetCount();++i)
		{
			int pos = ActiveFormatsString.find(AllFormatsArray[i]+_T(';'));
			pList->Check(i, pos!=-1);
		}
	}
	void OnOK(wxCommandEvent& WXUNUSED(event))
	{
		active.clear();
		unsigned i;
		for(i=0;i< pList->GetCount();++i)
		{
			if(pList->IsChecked(i))
				active.Append(pList->GetString(i) + _T(';'));
		}
		EndModal(wxID_OK);
	}
};
BEGIN_EVENT_TABLE(SelFormatsDialog,wxDialog)
	EVT_BUTTON(wxID_OK,SelFormatsDialog::OnOK)
END_EVENT_TABLE()

class ActiveFormats
{
private:
	wxString ActiveFormatsString;
	wxArrayString AllFormatsArray;
	wxArrayInt    ActiveFormatsFlags;
public:
	bool ReadConfig(wxConfig& config) //Called in OBGUIFrame constructor
	{
		bool chk=false;
		config.Read(_T("ActiveFormats"),&ActiveFormatsString,wxEmptyString);
		if(!ActiveFormatsString.IsEmpty())
			config.Read(_T("UseRestrictedFormats"), &chk, false);
		return chk;
	}
	void Clear() //Called in OBGUIFrame::GetAvailableFormats
	{
		AllFormatsArray.clear();
	}

	bool Add(wxString& format)//Called in OBGUIFrame::GetAvailableFormats
	{
		AllFormatsArray.Add(format);
		return (ActiveFormatsString.find(format+_T(';')) != -1);
	}

	void WriteConfig(wxConfig& config) //Called in OBGUIFrame::OnClose()
	{
		config.Write(_T("ActiveFormats"), ActiveFormatsString);
	}

	bool SelectFormats() // Called in OBGUIFrame::OnSelectFormats
	{
		SelFormatsDialog dialog(AllFormatsArray, ActiveFormatsString);
//		dialog.Connect(wxID_OK, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(SelFormatsDialog::OnOK));
		return(dialog.ShowModal()==wxID_OK);
	}
};

#endif

