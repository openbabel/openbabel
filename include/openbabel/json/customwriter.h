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

#include <string>

namespace Json {

   /** \brief Writes a Value in <a HREF="http://www.json.org">JSON</a> format with custom formatting.
    *
    * The JSON document is written according to the rules specified in the constructor. Objects and
    * arrays are printed on a single line if they are below a certain length, otherwise they are 
    * indented. It is possible to output invalid json if the customizable parameters are specified
    * incorrectly. Set maxWidth to 0 to print output on a single line. 
    *
    * \sa Reader, Value
    */
   class JSON_API CustomWriter : public Writer
   {
   public:
      CustomWriter( std::string opencurly = "{",
                    std::string closecurly = "}",
                    std::string opensquare = "[",
                    std::string closesquare = "]",
                    std::string colon = ":",
                    std::string comma = ",",
                    std::string indent = "  ",
                    int maxWidth = 74);
      virtual ~CustomWriter(){}

   public: // overridden from Writer
      virtual std::string write( const Value &root );

   private:
      void writeValue( const Value &value, std::string &doc, bool forceSingleLine );
      void indent();
      void unindent();

      std::string document_;
      std::string indentString_;
      std::string opencurly_;
      std::string closecurly_;
      std::string opensquare_;
      std::string closesquare_;
      std::string colon_;
      std::string comma_;
      std::string indent_;
      int maxWidth_;
   };
   
}
