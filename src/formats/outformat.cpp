/**********************************************************************
Copyright (C) 2001-2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include <sstream>

using namespace std;
namespace OpenBabel
{

  // The ".out" format:
  // Detect GAMESS, Q-Chem, Gaussian, or MOPAC output files
  class OutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    OutputFormat()
    {
      OBConversion::RegisterFormat("out", this);
      OBConversion::RegisterFormat("output", this);
      OBConversion::RegisterFormat("log", this);
      OBConversion::RegisterFormat("dat", this);
    }

    virtual const char* Description() //required
    {
      return
        "Generic Output file format\n"
        "Read ADF, Gaussian, GAMESS, Q-Chem, MOPAC, etc. file.out"
        "files by detecting contents\n"
        "Read Options e.g. -as\n"
        "   s  Output single bonds only\n"
        "   b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.sourceforge.net/wiki/Output";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  OutputFormat theOutputFormat;

  /////////////////////////////////////////////////////////////////
  bool OutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    // so we want to read through the file until we can figure out
    // what program actually created it
    // if we get to the end, emit a warning
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    OBFormat *pFormat = NULL;

    // the detection strings are from the Chemical MIME project
    // http://chemical-mime.sourceforge.net/chemical-mime-data.html
    while (ifs.getline(buffer,BUFF_SIZE)) {
      if (strstr(buffer,"GAMESS execution script") != NULL) {
        // GAMESS output
        pFormat = pConv->FindFormat("gamout");
        break;
      } else if (strstr(buffer,"Gaussian, Inc.") != NULL) {
        // Gaussian output
        pFormat = pConv->FindFormat("g03");
        break;
      } else if (strstr(buffer,"GENERAL UTILITY LATTICE PROGRAM") != NULL) {
        // GULP output -- not currently supported
        continue;
      } else if (strstr(buffer,"MOPAC") != NULL) {
        // MOPAC output
        pFormat = pConv->FindFormat("mopout");
        break;
      } else if (strstr(buffer,"Welcome to Q-Chem") != NULL) {
        // Q-Chem output
        pFormat = pConv->FindFormat("qcout");
        break;
      } else if (strstr(buffer,"Amsterdam Density Functional") != NULL) {
        // ADF output
        pFormat = pConv->FindFormat("adfout");
        break;
      } else if (strstr(buffer,"Northwest Computational Chemistry") != NULL) {
        // NWChem output
        pFormat = pConv->FindFormat("nwo");
        break;
      } else if (strstr(buffer,"MPQC: Massively Parallel Quantum Chemistry") != NULL) {
        // MPQC output
        pFormat = pConv->FindFormat("mpqc");
        break;
      } else if (strstr(buffer,"PROGRAM SYSTEM MOLPRO") != NULL) {
        // MOLPRO output
        pFormat = pConv->FindFormat("mpo");
        break;
      } else if ((strstr(buffer,"Schrodinger, Inc.") != NULL) &&
                 (strstr(buffer,"Jaguar") != NULL)) {
        // Jaguar
        pFormat = pConv->FindFormat("jout");
        break;
      }
    }

    if (pFormat) {
      ifs.seekg (0, ios::beg); // reset the stream to the beginning
      return pFormat->ReadMolecule(pOb, pConv);
    }

    obErrorLog.ThrowError(__FUNCTION__,
                          "Problems reading an output file: Could not determine the format of this file. Please report it to the openbabel-discuss @ lists.sourceforge.net mailing list.", obError);
    return(false); // we couldn't figure out the format
  }

} //namespace OpenBabel
