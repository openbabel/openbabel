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

#include "smarts.h"
#include "parsmart.h"

//
//Code for SMARTS parser
//

namespace OpenBabel
{

bool OESmartsParser::Parse(OESmartsPattern &sp,string &s)
{
	return(Parse(sp,s.c_str()));
}

bool OESmartsParser::Parse(OESmartsPattern &sp,const char *buf)
{
	//reset closure info if necessary

	int         vb;
	int         idx;
	OENode     *node;
	OEEdgeBase *cb;
	OEExprBase *eexpr=NULL;
	OEExprBase *vexpr=NULL;

	_prev = NULL;
	_vprev.clear();

	for (_ptr=buf;*_ptr;_ptr++)
      switch(*_ptr)
	  {
		//disconnected structure
        case '.': 
			if (!_vprev.empty())
			{
				ReportError();
				return(false);
			}
			_prev=NULL; 
			break;

		//start of branch
		case '(': 
			if (_prev == NULL)
			{
				ReportError();
				return(false);
			}

			_vprev.push_back(_prev); 
			break;

		//end of branch
		case ')': 
			if (_vprev.empty())
			{
				ReportError();
				return(false);
			}
			_prev = _vprev.back(); 
			_vprev.pop_back(); 
			break;

		//complex atom expression
		case '[': 

			_ptr++; //skip open bracket
			vexpr = ParseAtomExpr(0);
			if (!vexpr) return(false);
			if (*_ptr == ':') vb = GetVectorBinding();
			else              vb = 0;

			if (*_ptr != ']' || !vexpr)
			{
				if (vexpr) delete vexpr;
				ReportError();
				return(false);
			}

			node = sp.NewNode(vexpr);
			node->SetVectorBinding(vb);
			if (_prev) //bond to previous
			{
				if (eexpr) sp.NewEdge(_prev,node,eexpr);
				else       sp.NewEdge(_prev,node,new OEDefaultEdgeExpr);
			}

			_prev = node;
			eexpr = (OEExprBase*)NULL;
			break; 

		//bond expression
		case '-':	case '=':	case '#':  
		case ':':	case '/':	case '\\':
		case '!':	case '~':	case '@':
		eexpr = ParseBondExpr(0); 
		if (!eexpr) //cleanup and bail
			return(false);
		break;

		//double ring closure digit
		case '%': 
			break;

		//single ring closure digit
		case '0': case '1': 
		case '2': case '3': case '4': case '5':
		case '6': case '7': case '8': case '9': //ring closure

			if (!_prev)
			{
				ReportError();
				return(false);
			}
			idx = (int)((*_ptr)-'0');

			if ((cb = GetClosure(idx))) //closure already exists
			{
				cb->SetEnd(_prev);
				_prev->AddEdge(cb);
				if (eexpr) ((OEEdge*)cb)->ReplaceExpr(eexpr);
			}
			else //add closure bonds
			{
				if (eexpr) cb = sp.NewEdge(_prev,NULL,eexpr);
				else       cb = sp.NewEdge(_prev,NULL,new OEDefaultEdgeExpr);
				AddClosure(cb,idx);
			}
			eexpr = (OEExprBase*)NULL;
			break;

		//simple atom expression
		default: 

			if (!(vexpr = ParseSimpleAtomPrimitive()))
			{
				ReportError();
				return(false);
			}

			node = sp.NewNode(vexpr);
			if (_prev) //bond to previous
			{
				if (eexpr) sp.NewEdge(_prev,node,eexpr);
				else       sp.NewEdge(_prev,node,new OEDefaultEdgeExpr);
			}

			_prev = node;
			eexpr = (OEExprBase*)NULL;
		} // end switch

	sp.PrepForMatch();
	sp.SetSMARTS(buf);
	return(true);
}

OEExprBase *OESmartsParser::ParseAtomExpr( int level )
{
    const char *prev;
    OEExprBase *expr1;
    OEExprBase *expr2;

    switch( level )
    {
	case(0): /* Low Precedence Conjunction */

		if( !(expr1=ParseAtomExpr(1)) ) return((OEExprBase*)NULL);

		while( *_ptr == ';' )
		{
			_ptr++;
			if( !(expr2=ParseAtomExpr(1)) )
			{
				delete expr1;
				return((OEExprBase*)NULL);
			}
			expr1 = new OEAndExpr(expr1,expr2);
		}
		return(expr1);

	case(1): /* Disjunction */

		if( !(expr1=ParseAtomExpr(2)) ) return((OEExprBase*)NULL);

		while( *_ptr == ',' )
		{
			_ptr++;
			if( !(expr2=ParseAtomExpr(2)) )
			{
				delete expr1;
				return((OEExprBase*)NULL);
			}
			expr1 = new OEOrExpr(expr1,expr2);
		}
		return(expr1);

	case(2): /* High Precedence Conjunction */

		if( !(expr1=ParseAtomExpr(3)) ) return((OEExprBase*)NULL);
			
		while( (*_ptr!=']') && (*_ptr!=';') && *_ptr != ',' && *_ptr )
		{   
			if(*_ptr == '&') _ptr++;
			prev = _ptr;
			if( !(expr2=ParseAtomExpr(3)) )
			{
				if( prev != _ptr )
				{
					delete expr1;
					return((OEExprBase*)NULL);
				} 
				else return(expr1);
			}
			expr1 = new OEAndExpr(expr1,expr2);
		}
		return(expr1);

	case(3): /* Negation or Primitive */
		if(*_ptr == '!')
		{
			_ptr++;
			if( !(expr1=ParseAtomExpr(3)) )
				return( (OEExprBase*)NULL);
			return(new OENotExpr(expr1));
		}
		return(ParseComplexAtomPrimitive());
    }
    return((OEExprBase*)NULL);
}

OEExprBase *OESmartsParser::ParseSimpleAtomPrimitive()
{
	if (islower(*_ptr))
	{
		switch( *_ptr)
		{
		case 'a':  return(new OEAromaticExpr(true));
		case 'c':  return(new OEAromElemExpr(6,true));
		case 'n':  return(new OEAromElemExpr(7,true));
		case 'o':  return(new OEAromElemExpr(8,true));
		case 'p':  return(new OEAromElemExpr(15,true));
		case 's':  return(new OEAromElemExpr(16,true));
		}
	}
	else
	{
		const char *next = _ptr; next++;

		switch( *_ptr)
		{
			case 'C':  if( *next == 'l' )
					   {
						   _ptr++;
						   return(new OEElementExpr(17));
					   }
					   return(new OEAromElemExpr(6,false));

			case 'N':  return(new OEAromElemExpr(7,false));
			case 'O':  return(new OEAromElemExpr(8,false));
			case 'S':  return(new OEAromElemExpr(16,false));
			case 'P':  return(new OEAromElemExpr(15,false));
			case '*':  return(new OEConstExpr);
			case 'A':  return(new OEAromaticExpr(false));
			case 'B':  if( *next == 'r' )
					   {
						   _ptr++;
						   return(new OEElementExpr(35));
					   }
					   return(new OEElementExpr(5));

			case 'F':  return(new OEElementExpr(9 ));
			case 'I':  return(new OEElementExpr(53));
		}
	}

    return(NULL);
}

OEExprBase *OESmartsParser::ParseComplexAtomPrimitive()
{
    int index;

	if (isalpha(*_ptr))
	{
		if (islower(*_ptr))
		{
			switch(*_ptr++)
			{
			case 'a':
				if( *_ptr == 's' )
				{
					_ptr++;
					return(new OEAromElemExpr(33,true));
				}
				return(new OEAromaticExpr(true));

			case 'c': return(new OEAromElemExpr(6,true));

			case 'h':
				if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
				}
				else index = 1;
				return(new OEImplicitExpr(index));

			case 'n':  return(new OEAromElemExpr(7,true));
			case 'o':  return(new OEAromElemExpr(8,true));
			case 'p':  return(new OEAromElemExpr(15,true));

			case 'r':
				if( isdigit(*_ptr) )
				{
					index = 0;
					while(isdigit(*_ptr))
						index = index*10 + ((*_ptr++)-'0');
					if( index == 0 )
						return(new OERingExpr(0));
					return(new OESizeExpr(index));
				}
				return(new OERingExpr);

			case 's':  
				if(*_ptr == 'i')
				{
					_ptr++;
					return(new OEAromElemExpr(14,true));
				}
				return(new OEAromElemExpr(16,true));

			case 'v':
				if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OEValenceExpr(index));
				}
				return((OEExprBase*)NULL);
			} //switch
		} //if (islower())
		else
		{
			switch (*_ptr++)
			{
			case 'C':  
				switch( *_ptr++ )
				{
				case 'a':  return(new OEElementExpr(20));
				case 'd':  return(new OEElementExpr(48));
				case 'e':  return(new OEElementExpr(58));
				case 'f':  return(new OEElementExpr(98));
				case 'l':  return(new OEElementExpr(17));
				case 'm':  return(new OEElementExpr(96));
				case 'o':  return(new OEElementExpr(27));
				case 'r':  return(new OEElementExpr(24));
				case 's':  return(new OEElementExpr(55));
				case 'u':  return(new OEElementExpr(29));
				}
				_ptr--;
				return(new OEAromElemExpr(6,false));
				
			case 'H':
				switch(*_ptr++)
				{
				case 'e': return(new OEElementExpr(2));
				case 'f': return(new OEElementExpr(72));
				case 'g': return(new OEElementExpr(80));
				case 'o': return(new OEElementExpr(67));
				}
				_ptr--;
				
				if(isdigit(*_ptr))
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OEHCountExpr(index));
				}	
				return(new OEHCountExpr(1));

			case 'A':
				switch(*_ptr++)
				{
				case 'c':  return(new OEElementExpr(89));
				case 'g':  return(new OEElementExpr(47));
				case 'l':  return(new OEElementExpr(13));
				case 'm':  return(new OEElementExpr(95));
				case 'r':  return(new OEElementExpr(18));
				case 's':  return(new OEElementExpr(33));
				case 't':  return(new OEElementExpr(85));
				case 'u':  return(new OEElementExpr(79));
				}
				_ptr--;
                return(new OEAromaticExpr(false));

			case 'B':
				switch( *_ptr++ )
				{
				case 'a':  return(new OEElementExpr(56));
				case 'e':  return(new OEElementExpr( 4));
				case 'i':  return(new OEElementExpr(83));
				case 'k':  return(new OEElementExpr(97));
				case 'r':  return(new OEElementExpr(35));
				}
				_ptr--;
				return (new OEElementExpr(5));

			case 'D':
				if( *_ptr == 'y' )
				{
					_ptr++;
					return(new OEElementExpr(66));
				}
				else if(isdigit(*_ptr))
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OEDegreeExpr(index));
				}
				return((OEExprBase*)NULL);

			case 'E':
				switch(*_ptr++)
				{
				case 'r': return(new OEElementExpr(68));
				case 's': return(new OEElementExpr(99));
				case 'u': return(new OEElementExpr(63));
				}
				return((OEExprBase*)NULL);

			case 'F':
				switch(*_ptr++)
				{
				case 'e': return(new OEElementExpr(26));
				case 'm': return(new OEElementExpr(100));
				case 'r': return(new OEElementExpr(87));
				}
				_ptr--;
				return (new OEElementExpr(9));

			case 'G':
				switch(*_ptr++)
				{
				case 'a': return(new OEElementExpr(31));
				case 'd': return(new OEElementExpr(64));
				case 'e': return(new OEElementExpr(32));
				}
				return((OEExprBase*)NULL);

			case 'I':
				switch(*_ptr++)
				{
				case 'n': return(new OEElementExpr(49));
				case 'r': return(new OEElementExpr(77));
				}				
				_ptr--;
				return(new OEElementExpr(53));

			case 'K':  
				if( *_ptr++ == 'r' )
					return(new OEElementExpr(36));
				_ptr--;
				return(new OEElementExpr(19));

			case 'L':  
				switch(*_ptr++)
				{
				case 'a': return(new OEElementExpr(57));
				case 'i': return(new OEElementExpr(3));
				case 'r': return(new OEElementExpr(103));
				case 'u': return(new OEElementExpr(71));
				}
				return((OEExprBase*)NULL);
			
			case 'M':
				switch(*_ptr++)
				{
				case 'd': return(new OEElementExpr(101));
				case 'g': return(new OEElementExpr(12));
				case 'n': return(new OEElementExpr(25));
				case 'o': return(new OEElementExpr(42));
				}
				return((OEExprBase*)NULL);

			case 'N':
				switch( *_ptr++ )
				{
				case 'a':  return(new OEElementExpr( 11));
				case 'b':  return(new OEElementExpr( 41));
				case 'd':  return(new OEElementExpr( 60));
				case 'e':  return(new OEElementExpr( 10));
				case 'i':  return(new OEElementExpr( 28));
				case 'o':  return(new OEElementExpr(102));
				case 'p':  return(new OEElementExpr( 93));
				}
				_ptr--;
				return(new OEAromElemExpr(7,false));

			case 'O':
				if(*_ptr == 's')
				{
					_ptr++;
					return(new OEElementExpr(76));
				}
                return(new OEAromElemExpr(8,false));

			case 'P':
				switch(*_ptr++)
				{
				case 'a':  return(new OEElementExpr(91));
				case 'b':  return(new OEElementExpr(82));
				case 'd':  return(new OEElementExpr(46));
				case 'm':  return(new OEElementExpr(61));
				case 'o':  return(new OEElementExpr(84));
				case 'r':  return(new OEElementExpr(59));
				case 't':  return(new OEElementExpr(78));
				case 'u':  return(new OEElementExpr(94));
				}
				_ptr--;
				return(new OEElementExpr(15));

			case 'R':
				switch( *_ptr++ )
				{
				case 'a':  return(new OEElementExpr(88));
				case 'b':  return(new OEElementExpr(37));
				case 'e':  return(new OEElementExpr(75));
				case 'h':  return(new OEElementExpr(45));
				case 'n':  return(new OEElementExpr(86));
				case 'u':  return(new OEElementExpr(44));
				}
				_ptr--;
                    
				if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
				}
				else index = -1;
				return(new OERingExpr(index));

			case 'S':  
				switch( *_ptr++ )
				{   
				case 'b':  return(new OEElementExpr(51));
				case 'c':  return(new OEElementExpr(21));
				case 'e':  return(new OEElementExpr(34));
				case 'i':  return(new OEElementExpr(14));
				case 'm':  return(new OEElementExpr(62));
				case 'n':  return(new OEElementExpr(50));
				case 'r':  return(new OEElementExpr(38));
				}
				_ptr--;
				return(new OEAromElemExpr(16,false));

			case 'T':
				switch( *_ptr++ )
				{
				case 'a':  return(new OEElementExpr(73));
				case 'b':  return(new OEElementExpr(65));
				case 'c':  return(new OEElementExpr(43));
				case 'e':  return(new OEElementExpr(52));
				case 'h':  return(new OEElementExpr(90));
				case 'i':  return(new OEElementExpr(22));
				case 'l':  return(new OEElementExpr(81));
				case 'm':  return(new OEElementExpr(69));
				}
				_ptr--;
				return((OEExprBase*)NULL);

			case 'U':  return(new OEElementExpr(92));
			case 'V':  return(new OEElementExpr(23));
			case 'W':  return(new OEElementExpr(74));

			case 'X':
				if( *_ptr == 'e' )
				{
					_ptr++;  return(new OEElementExpr(54));
				}
				else if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OEConnectExpr(index));
				}
				return((OEExprBase*)NULL);

			case 'Y':
				if( *_ptr == 'b' )
				{
					_ptr++;
					return(new OEElementExpr(70));
				}
				return(new OEElementExpr(39));

			case 'Z':
				switch( *_ptr++ )
				{
				case 'n':  return(new OEElementExpr(30));
				case 'r':  return(new OEElementExpr(40));
				}
				_ptr--;
				return((OEExprBase*)NULL);
			} //switch
		} //if (islower)
	}
	else //isalpha
	{
		switch( *_ptr++)
		{   
		case '#':  
			if(!isdigit(*_ptr)) return((OEExprBase*)NULL);
			index = 0;
			while(isdigit(*_ptr)) index = index*10 + ((*_ptr++)-'0');
			
			if( index > ELEMMAX )
			{
				_ptr--;
				return((OEExprBase*)NULL);
			}
			else if( !index )
				return((OEExprBase*)NULL);
			return(new OEElementExpr(index));

		case '$':
			if( *_ptr != '(' ) 
				return((OEExprBase*)NULL);
			else
			{
				int count,size=0;
				const char *start = _ptr; start++;
				//advance _ptr to end of recursive smarts section
				for (count=1,size=0;*_ptr++ != '\0';size++)
				{
					if (*_ptr == '(') count++;
					if (*_ptr == ')') count--;
					if (count == 0) {_ptr++;break;}
				}
				if (*_ptr == '\0') return((OEExprBase*)NULL);
				if (size <= 0) return((OEExprBase*)NULL);
				//make a temp copy of the recursive section
				char *tmp = new char [size+1];
				memset(tmp,'\0',sizeof(char)*size+1);
				memcpy(tmp,start,sizeof(char)*size);

				OESmartsParser sp;
				OESmartsPattern *pat = new OESmartsPattern;
				if (!sp.Parse(*pat,tmp))
				{
					delete pat;
					return((OEExprBase*)NULL);
				}
				delete [] tmp;
				return(new OERecursExpr(pat));
			}
			break;

        case '*':  
			return(new OEConstExpr);

        case '+':
			if( isdigit(*_ptr) )
			{
				index = 0;
				while( isdigit(*_ptr) )
					index = index*10 + ((*_ptr++)-'0');
			} 
			else
			{
				index = 1;
				while( *_ptr == '+' )
				{
					_ptr++;
					index++;
				}
			}
			return(new OEPosChargeExpr(index));

        case('-'):
			if( isdigit(*_ptr) )
			{
				index = 0;
				while( isdigit(*_ptr) )
					index = index*10 + ((*_ptr++)-'0');
			}
			else
			{
				index = 1;
				while( *_ptr == '-' )
				{
					_ptr++;
					index++;
				}
			}
			return(new OENegChargeExpr(index));

        case '@':
			if (*_ptr == '@')
			{
				_ptr++;
				//_stereo = OE_ACLOCK;
			}
			else
				//_stereo = OE_CLOCK;
			if (*_ptr == 'H') _ptr++;

			return(new OENotExpr(new OEConstExpr));

        case '^': 
			if (isdigit(*_ptr))
			{
				index = 0;
				while( isdigit(*_ptr) )
					index = index*10 + ((*_ptr++)-'0');
				return(new OEHybExpr(index));
			}
			else
				return(new OEHybExpr(1));

        case('0'): case('1'): case('2'): case('3'): case('4'):
        case('5'): case('6'): case('7'): case('8'): case('9'):
			index = _ptr[-1]-'0';
			while( isdigit(*_ptr) )
				index = index*10 + ((*_ptr++)-'0');
			return(new OEMassExpr(index));
		}
	}

    _ptr--;
    return((OEExprBase*)NULL);
}

int OESmartsParser::GetVectorBinding()
{
  int vb=0;

  _ptr++; //skip colon
  if(isdigit(*_ptr))
    {
      vb = 0;
      while( isdigit(*_ptr) )
		  vb = vb*10 + ((*_ptr++)-'0');
    }

  return(vb);
}

OEExprBase *OESmartsParser::ParseBondExpr( int level )
{
    const char *prev;
    OEExprBase *expr1;
    OEExprBase *expr2;

    switch( level )
    {
	case 0: /* Low Precedence Conjunction */

		if( !(expr1=ParseBondExpr(1)) ) return(NULL);

		while( *_ptr == ';' )
		{   
			_ptr++;
			if( !(expr2=ParseBondExpr(1)) )
			{   
				delete expr1;
				return((OEExprBase*)NULL);
			}
			expr1 = new OEAndExpr(expr1,expr2);
		}
		_ptr--;
		return(expr1);

	case 1: /* Disjunction */

		if( !(expr1=ParseBondExpr(2)) ) return(NULL);

		while( *_ptr == ',' )
		{
			_ptr++;
			if( !(expr2=ParseBondExpr(2)) )
			{
				delete expr1; 
				return((OEExprBase*)NULL);
			}
			expr1 = new OEOrExpr(expr1,expr2);
		}
		return(expr1);

	case 2: /* High Precedence Conjunction */

		if( !(expr1=ParseBondExpr(3)) ) return((OEExprBase*)NULL);

		while(*_ptr !=']' && *_ptr !=';' && *_ptr !=',' && *_ptr )
		{
			if( *_ptr == '&' ) _ptr++;
			prev = _ptr;
			if( !(expr2=ParseBondExpr(3)) )
			{   
				if( prev != _ptr )
				{
					delete expr1;
					return((OEExprBase*)NULL);
				} 
				else return(expr1);
			}
			expr1 = new OEAndExpr(expr1,expr2); 
		}
		return(expr1);

        case(3): /* Negation or Primitive */
                 
			if( *_ptr == '!' )
			{   
				_ptr++;
				if( !(expr1=ParseBondExpr(3)) ) return((OEExprBase*)NULL);
				return(new OENotExpr(expr1));
			}
			return(ParseBondPrimitive());
    }

    return((OEExprBase*)NULL);
}

OEExprBase *OESmartsParser::ParseBondPrimitive()
{
	switch(*_ptr++)
	{
	case '-': return(new OESingleExpr);
	case '=': return(new OEDoubleExpr);
	case '#': return(new OETripleExpr);
	case ':': return(new OEAromaticExpr(true));
	case '@': return(new OERingExpr);
	case '~': return(new OEConstExpr);
	case '/':
		if (*_ptr == '?')
		{
			_ptr++;
			return(new OEUpUnspecExpr);
		}
		return(new OEUpExpr);
	case '\\':
		if (*_ptr)
		{
			_ptr++;
			return(new OEDownUnspecExpr);
		}
		return(new OEDownExpr);
	}

	_ptr--;
	return((OEExprBase*)NULL);
}


OEEdgeBase *OESmartsParser::GetClosure(int idx)
{
	OEEdgeBase *edge;
	vector<pair<OEEdgeBase*,int> >::iterator i;
	for (i = _vclose.begin();i != _vclose.end();i++)
		if (idx == i->second)
		{
			edge = i->first;
			_vclose.erase(i);
			return(edge);
		}

	return(NULL);
}

void OESmartsParser::AddClosure(OEEdgeBase *edge,int idx)
{
	_vclose.push_back(pair<OEEdgeBase*,int> (edge,idx));
}

} //namespace OpenBabel

