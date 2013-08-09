/**********************************************************************
Copyright (C) 2013 by Matt Swain <m.swain@me.com>

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

#include <openbabel/json/json.h>
#include <openbabel/json/customwriter.h>

namespace Json {

CustomWriter::CustomWriter( std::string opencurly,
                            std::string closecurly,
                            std::string opensquare,
                            std::string closesquare,
                            std::string colon,
                            std::string comma,
                            std::string indent,
                            int maxWidth)
   : opencurly_( opencurly )
   , closecurly_( closecurly )
   , opensquare_( opensquare )
   , closesquare_( closesquare )
   , colon_( colon )
   , comma_( comma )
   , indent_( indent )
   , maxWidth_( maxWidth )
{
}


std::string 
CustomWriter::write( const Value &root )
{
   document_ = "";
   indentString_ = "";
   writeValue( root, document_, false );
   document_ += "\n";
   return document_;
}


void 
CustomWriter::writeValue( const Value &value, std::string &doc, bool forceSingleLine )
{
   switch ( value.type() )
   {
   case nullValue:
      doc += "null";
      break;
   case intValue:
      doc += valueToString( value.asLargestInt() );
      break;
   case uintValue:
      doc += valueToString( value.asLargestUInt() );
      break;
   case realValue:
      doc += valueToString( value.asDouble() );
      break;
   case stringValue:
      doc += valueToQuotedString( value.asCString() );
      break;
   case booleanValue:
      doc += valueToString( value.asBool() );
      break;
   case arrayValue:
      {
         bool isMulti = false;
         if (!forceSingleLine)
         {
            std::string valLine = "";
            writeValue( value, valLine, true);
            if (valLine.length() > maxWidth_)
            {
               isMulti = true;
            }
            else
            {
               doc += valLine;
               break;
            }
         }
         doc += opensquare_;
         if (isMulti)
            indent();
         for ( int index =0; index < value.size(); ++index )
         {
            if (isMulti)
            {
               doc += "\n";
               doc += indentString_;
            }
            writeValue( value[index], doc, false );
            if ( index < value.size()-1 )
               doc += comma_;
         }
         if (isMulti)
         {
            unindent();
            doc += "\n";
            doc += indentString_;
         }
         doc += closesquare_;
      }
      break;
   case objectValue:
      {
         bool isMulti = false;
         if (!forceSingleLine)
         {
            std::string valLine = "";
            writeValue( value, valLine, true);
            if (valLine.length() > maxWidth_)
            {
               isMulti = true;
            }
            else
            {
               doc += valLine;
               break;
            }
         }
         Value::Members members( value.getMemberNames() );
         doc += opencurly_;
         if (isMulti)
            indent();
         for ( Value::Members::iterator it = members.begin(); 
               it != members.end(); 
               ++it )
         {
            if (isMulti)
            {
               doc += "\n";
               doc += indentString_;
               
            }
            const std::string &name = *it;
            doc += valueToQuotedString( name.c_str() );
            doc += colon_;
            writeValue( value[name], doc, forceSingleLine );
            if ( !(it + 1 == members.end()) )
               doc += comma_;
         }
         if (isMulti)
         {
            unindent();
            doc += "\n";
            doc += indentString_;
         }
         doc += closecurly_;
      }
      break;
   }
}


void 
CustomWriter::indent()
{
   indentString_ += indent_;
}


void 
CustomWriter::unindent()
{
   int idSize = int(indent_.size());
   int idsSize = int(indentString_.size());
   if (idsSize >= idSize)
      indentString_.resize (idsSize - idSize);
}

}
