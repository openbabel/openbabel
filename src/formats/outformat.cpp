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
#include <openbabel/generic.h>

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

#include <sstream>

using namespace std;
namespace OpenBabel
{

  // The ".out" format:
  // Detect GAMESS, GAMESS-UK, Q-Chem, PWSCF, Gaussian,  or MOPAC output files
  // also detect ORCA output files now     --- added 14.03.2014 by D.Lenk
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
        "Automatically detect and read computational chemistry output files\n\n"
        "This format can be used to read ADF, Gaussian, GAMESS, PWSCF, Q-Chem,\n"
        "MOPAC, ORCA etc. output files by automatically detecting the file type.\n\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.org/wiki/Output";}; //optional

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
    std::string formatName;

    // the detection strings are from the Chemical MIME project
    // http://chemical-mime.sourceforge.net/chemical-mime-data.html
    while (ifs.getline(buffer,BUFF_SIZE)) {
      if ((strstr(buffer,"GAMESS execution script") != NULL) ||
          (strstr(buffer,"PC GAMESS") != NULL) ||
          (strstr(buffer,"GAMESS VERSION") != NULL)) {
        // GAMESS output
        formatName = "gamout";
        break;
      } else if (strstr(buffer,"===  G A M E S S - U K    === ") != NULL) {
        // GAMESS-UK output
        formatName = "gukout";
        break;
      } else if (strstr(buffer,"Gaussian, Inc") != NULL) {
        // Gaussian output
        formatName = "g03";
        break;
      } else if (strstr(buffer,"GENERAL UTILITY LATTICE PROGRAM") != NULL) {
        // GULP output -- not currently supported
        break;
      } else if (strstr(buffer,"MOPAC") != NULL) {
        // MOPAC output
        formatName = "mopout";
        break;
      } else if (strstr(buffer,"Program PWSCF") != NULL) {
        // PWSCF
        formatName = "pwscf";
        break;
      } else if (strstr(buffer,"Welcome to Q-Chem") != NULL) {
        // Q-Chem output
        formatName = "qcout";
        break;
      } else if (strstr(buffer,"Amsterdam Density Functional") != NULL) {
        // ADF output
        // Determine the kind of ADF output
        while (ifs.getline(buffer, BUFF_SIZE)) {
          if (strstr(buffer, "|     A D F     |") != NULL) {
            formatName = "adfout";
            break;
          } else if (strstr(buffer, "|     B A N D     |") != NULL) {
            formatName = "adfband";
            break;
          } else if (strstr(buffer, "|     D F T B     |") != NULL) {
            formatName = "adfdftb";
            break;
          } else if (strstr(buffer, "DFTB Engine") != NULL) {
            // "|     D F T B     |" is no longer printed in ADF 2018
            // Hopefully, "DFTB Engine" will work fine...
            formatName = "adfdftb";
            break;
          }
        }
        break;
      } else if (strstr(buffer,"Northwest Computational Chemistry") != NULL) {
        // NWChem output
        formatName = "nwo";
        break;
      } else if (strstr(buffer,"MPQC: Massively Parallel Quantum Chemistry") != NULL) {
        // MPQC output
        formatName = "mpqc";
        break;
      } else if (strstr(buffer,"PROGRAM SYSTEM MOLPRO") != NULL) {
        // MOLPRO output
        formatName = "mpo";
        break;
      } else if ((strstr(buffer,"Schrodinger, Inc.") != NULL) &&
                 (strstr(buffer,"Jaguar") != NULL)) {
        // Jaguar
        formatName = "jout";
        break;
      } else if (strstr(buffer, "ABINIT") != NULL) {
        // Abinit
        formatName = "abinit";
        break;
      } else if (strstr(buffer, "ACES2") != NULL) {
        // ACESII
        formatName = "acesout";
        break;
      } else if (strstr(buffer, "CRYSTAL06") != NULL ||
                 strstr(buffer, "CRYSTAL09") != NULL) {
        // CRYSTAL09
        formatName = "c09out";
        break;
      } else if (strstr(buffer, "* O   R   C   A *") != NULL) {
        // ORCA
        formatName = "orca";
        break;
      } else if (strstr(buffer, "WELCOME TO SIESTA") != NULL) {
        // SIESTA
        formatName = "siesta";
        break;
      }
    }

    // if we assigned something above, let's try to find it
    if (formatName.length())
      pFormat = pConv->FindFormat(formatName);

    if (pFormat) {
      ifs.seekg (0, ios::beg); // reset the stream to the beginning
      ifs.clear();
      bool success = pFormat->ReadMolecule(pOb, pConv);

      // Tag the molecule with the format (e.g., if a program wants to know the kind of "out" or "log" file)
      // We have to do this *after* ReadMolecule returns, or the data might be cleared
      if (pOb) {
        OBPairData *dp = new OBPairData;
        dp->SetAttribute("File Format");
        dp->SetValue(formatName);
        dp->SetOrigin(fileformatInput);
        pOb->SetData(dp);
      }

      return success;
    }

    obErrorLog.ThrowError(__FUNCTION__,
                          "Problems reading an output file: Could not determine the format of this file. Please report it to the openbabel-discuss @ lists.sourceforge.net mailing list.", obError);
    return(false); // we couldn't figure out the format
  }

} //namespace OpenBabel
