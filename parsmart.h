/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef PARSMART_H
#define PARSMART_H

#include "smarts.h"

namespace OpenBabel
{

//
//SMARTS Parser
//

class OESmartsParser
{
	const char                    *_ptr;
	int                            _stereo;
	int                            _vb;
	OENode                        *_prev;
	vector<OENode*>                _vprev;
	vector<pair<OEEdgeBase*,int> > _vclose;
public:
	int         GetVectorBinding();
	void        ReportError() {}
	void        AddClosure(OEEdgeBase*,int);
	bool        Parse(OESmartsPattern&,const char*);
	bool        Parse(OESmartsPattern&,string&);
	OEExprBase *ParseSimpleAtomPrimitive();
	OEExprBase *ParseComplexAtomPrimitive();
	OEExprBase *ParseBondPrimitive();
	OEExprBase *ParseAtomExpr(int);
	OEExprBase *ParseBondExpr(int);
	OEEdgeBase *GetClosure(int);
};

#define ELEMMAX 104
#define OE_CLOCK   1
#define OE_ACLOCK  2

} //namespace OpenBabel

#endif //PARSMART_H
