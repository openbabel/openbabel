/**********************************************************************
descriptor.h - Base class for molecular descriptors
 
Copyright (C) 2007 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_DESCRIPTOR_H
#define OB_DESCRIPTOR_H

#include <string>
#include <sstream>
#include <limits>

#include <openbabel/plugin.h>

namespace OpenBabel
{
  class OBBase; //Forward declaration; used only as pointer.

  // Class introduction in descriptor.cpp
  class OBAPI OBDescriptor : public OBPlugin
  {
    MAKE_PLUGIN(OBDescriptor)

    public:
      /** 
       * @return The plugin type. 
       */
      const char* TypeID() { return "descriptors"; }
      /** 
       * @param pOb The molecule.
       * 
       * @return The value of a numeric descriptor.
       */
      virtual double Predict(OBBase* pOb) { return std::numeric_limits<double>::quiet_NaN(); }
      /** 
       * @param pOb The molecule.
       * 
       * @return The value of the descriptor and adds it to the object's OBPairData.
       */
      double PredictAndSave(OBBase* pOb);
      /** 
       * @brief Provides a string value for non-numeric descriptors and returns 
       * NaN, or a string representation and returns a numeric value.
       * 
       * @param pOb The molecule.
       * @param svalue String representation for the descriptor.
       * 
       * @return NaN (for non-numeric descriptors) or a numeric value.
       */
      virtual double GetStringValue(OBBase* pOb, std::string& svalue);
      /** 
       * @brief Parses the filter stream for a relational expression and returns 
       * its result when applied to the chemical object.
       * 
       * @param pOb The molecule.
       * @param ss
       * @param noEval
       * 
       * @return 
       */
      virtual bool Compare(OBBase* pOb, std::istream& ss, bool noEval);
      /** 
       * @brief Write information on a plugin class to the string txt.
       * 
       * @param txt
       * @param param
       * @param ID
       * 
       * @return 
       *
       * If the parameter is a descriptor ID, displays the verbose description 
       * for that descriptor only  (e.g. babel -L descriptors HBA1).
       */
      virtual bool Display(std::string &txt, const char* param, const char* ID = NULL);
      /** 
       * @brief  Interprets the --filter option string and returns the combined 
       * result of all the comparisons it contains.
       * 
       * @param pOb The molecule.
       * @param ss
       * @param noEval
       * 
       * @return 
       */
      static bool FilterCompare(OBBase* pOb, std::istream& ss, bool noEval);
      /** 
       * @brief Reads list of descriptor IDs and calls PredictAndSave() for each.
       * 
       * @param pOb The molecule.
       * @param DescrList
       */
      static void AddProperties(OBBase* pOb, const std::string& DescrList);
      /** 
       * @brief Deletes all the OBPairDatas whose attribute names are in the 
       * list (if they exist). 
       * 
       * @param pOb The molecule.
       * @param DescrList OBPairData attribute names to delete.
       */
      static void DeleteProperties(OBBase* pOb, const std::string& DescrList);
      /** 
       * @brief Reads list of descriptor IDs and OBPairData names and returns a 
       * list of values, each precede by a space or the first character in the 
       * list if it is whitespace or punctuation.
       * 
       * @param pOb The molecule.
       * @param DescrList
       * 
       * @return 
       */
      static std::string GetValues(OBBase* pOb, const std::string& DescrList);
    protected:
      /** 
       * @brief Read an identifier from the filter string.
       * 
       * @param optionText
       * 
       * @return The identifier. 
       */
      static std::string GetIdentifier(std::istream& optionText);
      /** 
       * @brief 
       * 
       * @param optionText
       * @param ch1
       * @param ch2
       * @param svalue
       * 
       * @return 
       * @todo docs
       */
      static double ParsePredicate(std::istream& optionText, char& ch1, char& ch2, std::string& svalue);
      /** 
       * @brief Reads a string from the filter stream, optionally preceded by = or !=.
       * 
       * @param ss
       * @param result
       * 
       * @return False if != operator found, and true otherwise.
       */
      static bool ReadStringFromFilter(std::istream& ss, std::string& result);
      /** 
       * @brief Makes a comparison using the operator and a string read from 
       * the filter stream with a provided string.
       * 
       * @param optionText
       * @param s
       * @param noEval
       * @param NoCompOK
       * 
       * @return The result of the comparison and true if NoCompOK==true and 
       * there is no comparison operator.
       */
      static bool CompareStringWithFilter(std::istream& optionText, std::string& s, bool noEval, bool NoCompOK=false);
      /** 
       * @brief Treats _ as not a punctuation character. 
       * 
       * @param ch The char to test.
       * 
       * @return True if @p ch is a punctuation character.
       */
      static bool ispunctU(char ch) { return ispunct(ch) && ch != '_'; }
      /** 
       * @param pOb The molecule.
       * @param s The OBPairData attribute.
       * 
       * @return True if s (with or without _ replaced by spaces) is a PairData 
       * attribute. On return s is the form which matches.
       */
      static bool MatchPairData(OBBase* pOb, std::string& s);
  };

  template <class T>
  static bool DoComparison(char ch1, char ch2, T& val, T& filterval)
  {
    switch(ch1)
    {
      case (0):  //no comparison operator is same as =
      case('='):
        return val==filterval; //**needs a better floating point comparison**	
      case('!'):
        return val!=filterval; //**needs a better floating point comparison**	
      case('>'):              
        if(ch2=='=')
          return val>=filterval;
        else
          return val>filterval;
      case('<'):              
        if(ch2=='=')
          return val<=filterval;
        else
          return val<filterval;
    }
    
    return false;
  }

}//namespace

#endif

//! @file descriptor.h
//! @brief Base class for molecular descriptors
