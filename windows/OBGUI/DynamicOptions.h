// DynamicOptions.h: interface for the CDynamicOptions class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DYNAMICOPTIONS_H__21C0F058_352D_4FF2_9E19_F3F7FA9CF970__INCLUDED_)
#define AFX_DYNAMICOPTIONS_H__21C0F058_352D_4FF2_9E19_F3F7FA9CF970__INCLUDED_
#pragma warning (disable : 4786)

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

typedef multimap<char,CWnd*> OMapType;
class CDynamicOptions  
{
public:
	CDynamicOptions();
	void Construct(const char* OptionTxt, CWnd* pWnd, CRect& CtlRect );
	virtual ~CDynamicOptions();
	void InsertText(const char* txt, CWnd* pWnd, CRect& CtlRect );
	const char* GetOptions();
private:
	Clear();
	char* strcasestr(const char* haystack, const char* needle);
	
	OMapType OptionMap;
	string opts;
	UINT nID;
};
#endif // !defined(AFX_DYNAMICOPTIONS_H__21C0F058_352D_4FF2_9E19_F3F7FA9CF970__INCLUDED_)
