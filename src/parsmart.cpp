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

using namespace std;

namespace OpenBabel
{

OBSmartsParser::~OBSmartsParser()
{
}

bool OBSmartsParser::Parse(OBSmartsPattern &sp,string &s)
{
	return(Parse(sp,s.c_str()));
}

bool OBSmartsParser::Parse(OBSmartsPattern &sp,const char *buf)
{
	//reset closure info if necessary

	int         vb;
	int         idx;
	OBNode     *node;
	OBEdgeBase *cb;
	OBExprBase *eexpr=NULL;
	OBExprBase *vexpr=NULL;
	
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
			  if (vexpr) {delete vexpr; vexpr = NULL;}
			  ReportError();
			  return(false);
			}
			node = sp.NewNode(vexpr);
			node->SetVectorBinding(vb);
			if (_prev) //bond to previous
			{
				if (eexpr) sp.NewEdge(_prev,node,eexpr);
				else       sp.NewEdge(_prev,node,new OBDefaultEdgeExpr);
			}

			_prev = node;
			eexpr = (OBExprBase*)NULL;
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
				if (eexpr) ((OBEdge*)cb)->ReplaceExpr(eexpr);
			}
			else //add closure bonds
			{
				if (eexpr) cb = sp.NewEdge(_prev,NULL,eexpr);
				else       cb = sp.NewEdge(_prev,NULL,new OBDefaultEdgeExpr);
				AddClosure(cb,idx);
			}
			eexpr = (OBExprBase*)NULL;
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
				else       sp.NewEdge(_prev,node,new OBDefaultEdgeExpr);
			}

			_prev = node;
			eexpr = (OBExprBase*)NULL;
		} // end switch

	sp.PrepForMatch();
	sp.SetSMARTS(buf);
	return(true);
}

OBExprBase *OBSmartsParser::ParseAtomExpr( int level )
{
    const char *prev;
    OBExprBase *expr1;
    OBExprBase *expr2;

    switch( level )
    {
	case(0): /* Low Precedence Conjunction */

		if( !(expr1=ParseAtomExpr(1)) ) return((OBExprBase*)NULL);

		while( *_ptr == ';' )
		{
			_ptr++;
			if( !(expr2=ParseAtomExpr(1)) )
			{
			  delete expr1; expr1 = NULL;
				return((OBExprBase*)NULL);
			}
			expr1 = new OBAndExpr(expr1,expr2);
		}
		return(expr1);

	case(1): /* Disjunction */

		if( !(expr1=ParseAtomExpr(2)) ) return((OBExprBase*)NULL);

		while( *_ptr == ',' )
		{
			_ptr++;
			if( !(expr2=ParseAtomExpr(2)) )
			{
			  delete expr1; expr1 = NULL;
				return((OBExprBase*)NULL);
			}
			expr1 = new OBOrExpr(expr1,expr2);
		}
		return(expr1);

	case(2): /* High Precedence Conjunction */

		if( !(expr1=ParseAtomExpr(3)) ) return((OBExprBase*)NULL);
			
		while( (*_ptr!=']') && (*_ptr!=';') && *_ptr != ',' && *_ptr )
		{   
			if(*_ptr == '&') _ptr++;
			prev = _ptr;
			if( !(expr2=ParseAtomExpr(3)) )
			{
				if( prev != _ptr )
				{
				  delete expr1; expr1 = NULL;
					return((OBExprBase*)NULL);
				} 
				else return(expr1);
			}
			expr1 = new OBAndExpr(expr1,expr2);
		}
		return(expr1);

	case(3): /* Negation or Primitive */
		if(*_ptr == '!')
		{
			_ptr++;
			if( !(expr1=ParseAtomExpr(3)) )
				return( (OBExprBase*)NULL);
			return(new OBNotExpr(expr1));
		}
		return(ParseComplexAtomPrimitive());
    }
    return((OBExprBase*)NULL);
}

OBExprBase *OBSmartsParser::ParseSimpleAtomPrimitive()
{
	if (islower(*_ptr))
	{
		switch( *_ptr)
		{
		case 'a':  return(new OBAromaticExpr(true));
		case 'c':  return(new OBAromElemExpr(6,true));
		case 'n':  return(new OBAromElemExpr(7,true));
		case 'o':  return(new OBAromElemExpr(8,true));
		case 'p':  return(new OBAromElemExpr(15,true));
		case 's':  return(new OBAromElemExpr(16,true));
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
						   return(new OBElementExpr(17));
					   }
					   return(new OBAromElemExpr(6,false));

			case 'N':  return(new OBAromElemExpr(7,false));
			case 'O':  return(new OBAromElemExpr(8,false));
			case 'S':  return(new OBAromElemExpr(16,false));
			case 'P':  return(new OBAromElemExpr(15,false));
			case '*':  return(new OBConstExpr);
			case 'A':  return(new OBAromaticExpr(false));
			case 'B':  if( *next == 'r' )
					   {
						   _ptr++;
						   return(new OBElementExpr(35));
					   }
					   return(new OBElementExpr(5));

			case 'F':  return(new OBElementExpr(9 ));
			case 'I':  return(new OBElementExpr(53));
		}
	}

    return(NULL);
}

OBExprBase *OBSmartsParser::ParseComplexAtomPrimitive()
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
					return(new OBAromElemExpr(33,true));
				}
				return(new OBAromaticExpr(true));

			case 'c': return(new OBAromElemExpr(6,true));

			case 'h':
				if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
				}
				else index = 1;
				return(new OBImplicitExpr(index));

			case 'n':  return(new OBAromElemExpr(7,true));
			case 'o':  return(new OBAromElemExpr(8,true));
			case 'p':  return(new OBAromElemExpr(15,true));

			case 'r':
				if( isdigit(*_ptr) )
				{
					index = 0;
					while(isdigit(*_ptr))
						index = index*10 + ((*_ptr++)-'0');
					if( index == 0 )
						return(new OBRingExpr(0));
					return(new OBSizeExpr(index));
				}
				return(new OBRingExpr);

			case 's':  
				if(*_ptr == 'i')
				{
					_ptr++;
					return(new OBAromElemExpr(14,true));
				}
				return(new OBAromElemExpr(16,true));

			case 'v':
				if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OBValenceExpr(index));
				}
				return((OBExprBase*)NULL);
			} //switch
		} //if (islower())
		else
		{
			switch (*_ptr++)
			{
			case 'C':  
				switch( *_ptr++ )
				{
				case 'a':  return(new OBElementExpr(20));
				case 'd':  return(new OBElementExpr(48));
				case 'e':  return(new OBElementExpr(58));
				case 'f':  return(new OBElementExpr(98));
				case 'l':  return(new OBElementExpr(17));
				case 'm':  return(new OBElementExpr(96));
				case 'o':  return(new OBElementExpr(27));
				case 'r':  return(new OBElementExpr(24));
				case 's':  return(new OBElementExpr(55));
				case 'u':  return(new OBElementExpr(29));
				}
				_ptr--;
				return(new OBAromElemExpr(6,false));
				
			case 'H':
				switch(*_ptr++)
				{
				case 'e': return(new OBElementExpr(2));
				case 'f': return(new OBElementExpr(72));
				case 'g': return(new OBElementExpr(80));
				case 'o': return(new OBElementExpr(67));
				}
				_ptr--;
				
				if(isdigit(*_ptr))
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OBHCountExpr(index));
				}	
				return(new OBHCountExpr(1));

			case 'A':
				switch(*_ptr++)
				{
				case 'c':  return(new OBElementExpr(89));
				case 'g':  return(new OBElementExpr(47));
				case 'l':  return(new OBElementExpr(13));
				case 'm':  return(new OBElementExpr(95));
				case 'r':  return(new OBElementExpr(18));
				case 's':  return(new OBElementExpr(33));
				case 't':  return(new OBElementExpr(85));
				case 'u':  return(new OBElementExpr(79));
				}
				_ptr--;
                return(new OBAromaticExpr(false));

			case 'B':
				switch( *_ptr++ )
				{
				case 'a':  return(new OBElementExpr(56));
				case 'e':  return(new OBElementExpr( 4));
				case 'i':  return(new OBElementExpr(83));
				case 'k':  return(new OBElementExpr(97));
				case 'r':  return(new OBElementExpr(35));
				}
				_ptr--;
				return (new OBElementExpr(5));

			case 'D':
				if( *_ptr == 'y' )
				{
					_ptr++;
					return(new OBElementExpr(66));
				}
				else if(isdigit(*_ptr))
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OBDegreeExpr(index));
				}
				return((OBExprBase*)NULL);

			case 'E':
				switch(*_ptr++)
				{
				case 'r': return(new OBElementExpr(68));
				case 's': return(new OBElementExpr(99));
				case 'u': return(new OBElementExpr(63));
				}
				return((OBExprBase*)NULL);

			case 'F':
				switch(*_ptr++)
				{
				case 'e': return(new OBElementExpr(26));
				case 'm': return(new OBElementExpr(100));
				case 'r': return(new OBElementExpr(87));
				}
				_ptr--;
				return (new OBElementExpr(9));

			case 'G':
				switch(*_ptr++)
				{
				case 'a': return(new OBElementExpr(31));
				case 'd': return(new OBElementExpr(64));
				case 'e': return(new OBElementExpr(32));
				}
				return((OBExprBase*)NULL);

			case 'I':
				switch(*_ptr++)
				{
				case 'n': return(new OBElementExpr(49));
				case 'r': return(new OBElementExpr(77));
				}				
				_ptr--;
				return(new OBElementExpr(53));

			case 'K':  
				if( *_ptr++ == 'r' )
					return(new OBElementExpr(36));
				_ptr--;
				return(new OBElementExpr(19));

			case 'L':  
				switch(*_ptr++)
				{
				case 'a': return(new OBElementExpr(57));
				case 'i': return(new OBElementExpr(3));
				case 'r': return(new OBElementExpr(103));
				case 'u': return(new OBElementExpr(71));
				}
				return((OBExprBase*)NULL);
			
			case 'M':
				switch(*_ptr++)
				{
				case 'd': return(new OBElementExpr(101));
				case 'g': return(new OBElementExpr(12));
				case 'n': return(new OBElementExpr(25));
				case 'o': return(new OBElementExpr(42));
				}
				return((OBExprBase*)NULL);

			case 'N':
				switch( *_ptr++ )
				{
				case 'a':  return(new OBElementExpr( 11));
				case 'b':  return(new OBElementExpr( 41));
				case 'd':  return(new OBElementExpr( 60));
				case 'e':  return(new OBElementExpr( 10));
				case 'i':  return(new OBElementExpr( 28));
				case 'o':  return(new OBElementExpr(102));
				case 'p':  return(new OBElementExpr( 93));
				}
				_ptr--;
				return(new OBAromElemExpr(7,false));

			case 'O':
				if(*_ptr == 's')
				{
					_ptr++;
					return(new OBElementExpr(76));
				}
                return(new OBAromElemExpr(8,false));

			case 'P':
				switch(*_ptr++)
				{
				case 'a':  return(new OBElementExpr(91));
				case 'b':  return(new OBElementExpr(82));
				case 'd':  return(new OBElementExpr(46));
				case 'm':  return(new OBElementExpr(61));
				case 'o':  return(new OBElementExpr(84));
				case 'r':  return(new OBElementExpr(59));
				case 't':  return(new OBElementExpr(78));
				case 'u':  return(new OBElementExpr(94));
				}
				_ptr--;
				return(new OBElementExpr(15));

			case 'R':
				switch( *_ptr++ )
				{
				case 'a':  return(new OBElementExpr(88));
				case 'b':  return(new OBElementExpr(37));
				case 'e':  return(new OBElementExpr(75));
				case 'h':  return(new OBElementExpr(45));
				case 'n':  return(new OBElementExpr(86));
				case 'u':  return(new OBElementExpr(44));
				}
				_ptr--;
                    
				if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
				}
				else index = -1;
				return(new OBRingExpr(index));

			case 'S':  
				switch( *_ptr++ )
				{   
				case 'b':  return(new OBElementExpr(51));
				case 'c':  return(new OBElementExpr(21));
				case 'e':  return(new OBElementExpr(34));
				case 'i':  return(new OBElementExpr(14));
				case 'm':  return(new OBElementExpr(62));
				case 'n':  return(new OBElementExpr(50));
				case 'r':  return(new OBElementExpr(38));
				}
				_ptr--;
				return(new OBAromElemExpr(16,false));

			case 'T':
				switch( *_ptr++ )
				{
				case 'a':  return(new OBElementExpr(73));
				case 'b':  return(new OBElementExpr(65));
				case 'c':  return(new OBElementExpr(43));
				case 'e':  return(new OBElementExpr(52));
				case 'h':  return(new OBElementExpr(90));
				case 'i':  return(new OBElementExpr(22));
				case 'l':  return(new OBElementExpr(81));
				case 'm':  return(new OBElementExpr(69));
				}
				_ptr--;
				return((OBExprBase*)NULL);

			case 'U':  return(new OBElementExpr(92));
			case 'V':  return(new OBElementExpr(23));
			case 'W':  return(new OBElementExpr(74));

			case 'X':
				if( *_ptr == 'e' )
				{
					_ptr++;  return(new OBElementExpr(54));
				}
				else if( isdigit(*_ptr) )
				{
					index = 0;
					while( isdigit(*_ptr) )
						index = index*10 + ((*_ptr++)-'0');
					return(new OBConnectExpr(index));
				}
				return((OBExprBase*)NULL);

			case 'Y':
				if( *_ptr == 'b' )
				{
					_ptr++;
					return(new OBElementExpr(70));
				}
				return(new OBElementExpr(39));

			case 'Z':
				switch( *_ptr++ )
				{
				case 'n':  return(new OBElementExpr(30));
				case 'r':  return(new OBElementExpr(40));
				}
				_ptr--;
				return((OBExprBase*)NULL);
			} //switch
		} //if (islower)
	}
	else //isalpha
	{
		switch( *_ptr++)
		{   
		case '#':  
			if(!isdigit(*_ptr)) return((OBExprBase*)NULL);
			index = 0;
			while(isdigit(*_ptr)) index = index*10 + ((*_ptr++)-'0');
			
			if( index > ELEMMAX )
			{
				_ptr--;
				return((OBExprBase*)NULL);
			}
			else if( !index )
				return((OBExprBase*)NULL);
			return(new OBElementExpr(index));

		case '$':
			if( *_ptr != '(' ) 
				return((OBExprBase*)NULL);
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
				if (*_ptr == '\0') return((OBExprBase*)NULL);
				if (size <= 0) return((OBExprBase*)NULL);
				//make a temp copy of the recursive section
				char *tmp = new char [size+1];
				memset(tmp,'\0',sizeof(char)*size+1);
				memcpy(tmp,start,sizeof(char)*size);

				OBSmartsParser sp;
				OBSmartsPattern *pat = new OBSmartsPattern;
				if (!sp.Parse(*pat,tmp))
				{
				  delete pat; pat = NULL;
					return((OBExprBase*)NULL);
				}
				delete [] tmp;
				return(new OBRecursExpr(pat));
			}
			break;

        case '*':  
			return(new OBConstExpr);

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
			return(new OBPosChargeExpr(index));

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
			return(new OBNegChargeExpr(index));

        case '@':
			if (*_ptr == '@')
			{
				_ptr++;
				//_stereo = OB_ACLOCK;
			}
			else
				//_stereo = OB_CLOCK;
			if (*_ptr == 'H') _ptr++;

			return(new OBNotExpr(new OBConstExpr));

        case '^': 
			if (isdigit(*_ptr))
			{
				index = 0;
				while( isdigit(*_ptr) )
					index = index*10 + ((*_ptr++)-'0');
				return(new OBHybExpr(index));
			}
			else
				return(new OBHybExpr(1));

        case('0'): case('1'): case('2'): case('3'): case('4'):
        case('5'): case('6'): case('7'): case('8'): case('9'):
			index = _ptr[-1]-'0';
			while( isdigit(*_ptr) )
				index = index*10 + ((*_ptr++)-'0');
			return(new OBMassExpr(index));
		}
	}

    _ptr--;
    return((OBExprBase*)NULL);
}

int OBSmartsParser::GetVectorBinding()
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

OBExprBase *OBSmartsParser::ParseBondExpr( int level )
{
    const char *prev;
    OBExprBase *expr1;
    OBExprBase *expr2;

    switch( level )
    {
	case 0: /* Low Precedence Conjunction */

		if( !(expr1=ParseBondExpr(1)) ) return(NULL);

		while( *_ptr == ';' )
		{   
			_ptr++;
			if( !(expr2=ParseBondExpr(1)) )
			{   
			  delete expr1; expr1 = NULL;
				return((OBExprBase*)NULL);
			}
			expr1 = new OBAndExpr(expr1,expr2);
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
			  delete expr1; expr1 = NULL;
				return((OBExprBase*)NULL);
			}
			expr1 = new OBOrExpr(expr1,expr2);
		}
		return(expr1);

	case 2: /* High Precedence Conjunction */

		if( !(expr1=ParseBondExpr(3)) ) return((OBExprBase*)NULL);

		while(*_ptr !=']' && *_ptr !=';' && *_ptr !=',' && *_ptr )
		{
			if( *_ptr == '&' ) _ptr++;
			prev = _ptr;
			if( !(expr2=ParseBondExpr(3)) )
			{   
				if( prev != _ptr )
				{
				  delete expr1; expr1 = NULL;
					return((OBExprBase*)NULL);
				} 
				else return(expr1);
			}
			expr1 = new OBAndExpr(expr1,expr2); 
		}
		return(expr1);

        case(3): /* Negation or Primitive */
                 
			if( *_ptr == '!' )
			{   
				_ptr++;
				if( !(expr1=ParseBondExpr(3)) ) return((OBExprBase*)NULL);
				return(new OBNotExpr(expr1));
			}
			return(ParseBondPrimitive());
    }

    return((OBExprBase*)NULL);
}

OBExprBase *OBSmartsParser::ParseBondPrimitive()
{
	switch(*_ptr++)
	{
	case '-': return(new OBSingleExpr);
	case '=': return(new OBDoubleExpr);
	case '#': return(new OBTripleExpr);
	case ':': return(new OBAromaticExpr(true));
	case '@': return(new OBRingExpr);
	case '~': return(new OBConstExpr);
	case '/':
		if (*_ptr == '?')
		{
			_ptr++;
			return(new OBUpUnspecExpr);
		}
		return(new OBUpExpr);
	case '\\':
		if (*_ptr)
		{
			_ptr++;
			return(new OBDownUnspecExpr);
		}
		return(new OBDownExpr);
	}

	_ptr--;
	return((OBExprBase*)NULL);
}


OBEdgeBase *OBSmartsParser::GetClosure(int idx)
{
	OBEdgeBase *edge;
	vector<pair<OBEdgeBase*,int> >::iterator i;
	for (i = _vclose.begin();i != _vclose.end();i++)
		if (idx == i->second)
		{
			edge = i->first;
			_vclose.erase(i);
			return(edge);
		}

	return(NULL);
}

void OBSmartsParser::AddClosure(OBEdgeBase *edge,int idx)
{
	_vclose.push_back(pair<OBEdgeBase*,int> (edge,idx));
}

} //namespace OpenBabel

