/**********************************************************************
opsplit.cpp - A OBOp to write each molecule to a file with a name from the molecule, usually thetitle.

Copyright(C) 2011 by Chris Morley

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
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <fstream>
#include "deferred.h"
#include <openbabel/descriptor.h>

namespace OpenBabel
{

class OpSplit : public OBOp
{
public:
  OpSplit(const char* ID) : OBOp(ID, false){};
  const char* Description()
  { 
    return "[param] split, each filename from param;default title\n"
           "param is an optional property or descriptor.\n"
           "    obabel  in.sdf  -omol2  --split\n"
           "will output each molecule in a file named ``<title>.mol2``\n"
           "in the current directory.\n"
           "If an output filename is specified it is used to define the format\n"
           "and the path of the output files, but its name is ignored.\n"
           "If the molecule has no title or if it contains illegal characters,\n"
           "the input count (starting at 1) is used as the filename.\n"
           "The existing command:\n"
           "    obabel  in.sdf  -O .mol2  -m\n"
           "will give the same output as this, whatever the molecule titles.\n \n"

           "A descriptor or property, if present as a parameter, can be used\n"
           "to define the filename, e.g.\n"
           "    obabel  chembl_02.sdf.gz  -omol  --split chebi_id  -l8\n"
           "will produce a file for each of the first 8 molecules with filename\n"
           "the same as its sdf property ``chebi_id``.\n \n"

           "Filtering of the molecules or minor manipulation (e.g. make H explicit)\n"
           "can be done at the same time as splitting by adding further options\n"
           "to the command line.\n \n"

           "For files with Unix line-endings, copy format can be used to avoid\n"
           "uncertainties of interpretation. With SDF files, it can also be much\n"
           "faster if the -aT or -aP options are used when the title or\n"
           "properties respectively are used to derive the filename or in filtering.\n"
           "The following produces (in a subdirectory) an sdf file for each molecule,\n"
           "with a filename derived from its title:\n"
           "    obabel  in.sdf  -O sub/out.sdf  --split  -aT  -ocopy\n \n"

           "If a format specification like -osdf is used without a filename, the\n"
           "extension is chosen by OpenBabel from the alternatives for the format, and\n"
           "cannot be changed. So if you don't like .smiles, use a dummy filename.\n\n"
           ;
  }
  virtual bool WorksWith(OBBase* pOb)const { return true; } //all OBBase objects
  virtual bool Do(OBBase* pOb, const char* OptionText=nullptr, OpMap* pOptions=nullptr, OBConversion* pConv=nullptr);
private:
  int _inputCount;
  OBFormat* _realOutFormat;
  std::string _optionText, _outExt, _outPath;
  OBDescriptor* _pDesc;
};

/////////////////////////////////////////////////////////////////
OpSplit theOpSplit("split"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpSplit::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  if(!strcmp(OptionText, "inactive"))
  {
    ++_inputCount;
    return true;
  }

  if(!pConv)
    return false;
  
    if(pConv->IsFirstInput())
  {
    _inputCount=0;
    _optionText = OptionText; //because gets overwritten by "inactive"
    _pDesc = *OptionText ? OBDescriptor::FindType(OptionText) : nullptr;
    _realOutFormat = pConv->GetOutFormat();

    // If there is an output file specified, delete the file,close and invalidate the outstream so OBConversion is not confused.
    std::ofstream* oldfs = dynamic_cast<std::ofstream*>(pConv->GetOutStream());
    if(oldfs && oldfs->is_open())
    {
      oldfs->close();
      oldfs->setstate(std::ios::failbit);
      remove(pConv->GetOutFilename().c_str());
    }

    // If there is an output file name, use its path and its extension for the output files.
    // Otherwise use current directory and the ID of the output format.
    _outExt = _outPath = pConv->GetOutFilename(); //This may not be present yet in OBConversion
    std::string::size_type pos = _outPath.find_last_of("\\/");
    if(pos!=std::string::npos)
      _outPath.erase(pos+1);
    else
      _outPath.clear();
    pos = _outExt.rfind('.');
    if(pos!=std::string::npos)
      _outExt.erase(0, pos+1);
    else
      _outExt = _realOutFormat->GetID();

    /* Need to call this op after all other options because it outputs directly
       rather than letting OBConversion do the writing. Use the trick of
       the op deactivating itself and setting output to DeferredFormat. When this
       is written to, it immediately calls back this function, which writes to
       an individual file of the real output format.
    */
    pConv->AddOption(GetID(),OBConversion::GENOPTIONS,"inactive");//removing messes up DoOps()
    new DeferredFormat(pConv, this, true); //it will delete itself
    return true;
  }

  //The following is called from Deferred Format

  std::ofstream ofs;
  std::stringstream filename;
  filename.str("");
  filename << _outPath;
  std::string name;

  if(!_pDesc && _optionText.empty()) //no param, use title
  {
    const char* pFilename = pOb->GetTitle();
    name = pFilename;
  }
  else if(_pDesc) //use descriptor
  {
    std::string s;
    _pDesc->GetStringValue(pOb, s);
    name = s;
  }
  else //use OBPairData property
  {
    if(pOb->HasData(_optionText))
      name = pOb->GetData(_optionText)->GetValue();
  }

  bool filenameok = !name.empty() && name.find_first_of("/\\:*|?\"") == std::string::npos;
  if(filenameok)
  {
    filename << name << '.' << _outExt;
    std::string fname(filename.str()); //DEBUG
    ofs.open(fname.c_str());
  }
  if(!filenameok || !ofs)
  {
    //filename is empty or contain an illegal character or failed to open. Use "n.xxx"
    obErrorLog.ThrowError(__FUNCTION__,
      "The fallback filename, based on input index, has been used for at least one object.", obWarning, onceOnly);
    std::stringstream ss;
    ss << _outPath << _inputCount << '.' << _outExt;
    ofs.clear();
    if(ofs.is_open())
      ofs.close();
    ofs.open(ss.str().c_str());
    if(!ofs)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open file " + ss.str(), obError);
       pConv->SetOneObjectOnly(); //stop all conversion
      return false;
    }
  }
  pConv->SetOutStream(&ofs);
  //Output the object now - do not use the normal queue of two.
  _realOutFormat->WriteChemObject(pConv);
  //Correct for second call to GetChemObject(), which increments index, in above.
  //GetChemObject() already called in DeferredFormat
  pConv->SetOutputIndex(pConv->GetOutputIndex()-1);

  ofs.close();
  return false; //do not store anything in DeferredFormat
}

}//namespace
