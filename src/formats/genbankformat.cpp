/**********************************************************************
genbankformat.cpp - Conversion from genbank format.
Copyright (C) Scarlet Line 2008

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
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>

#include <iostream>

using namespace std;
namespace OpenBabel
{
  class GenBankFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GenBankFormat()
    { // Copied from the Chemical MIME Page at http://www.ch.ic.ac.uk/chemime/
      OBConversion::RegisterFormat("gen", this, "chemical/x-genbank");
      OBConversion::RegisterFormat("embl", this);
      OBConversion::RegisterFormat("ddbj", this);

      OBConversion::RegisterOptionParam("s", this);
      OBConversion::RegisterOptionParam("b", this);
    }

    virtual const char* Description() //required
    {
      return
        "GenBank, DDBJ, EMBL Flat File format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    }

    virtual const char* SpecificationURL()
    { return "http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html";} //optional
    // Genbank itself is here: http://www.ncbi.nlm.nih.gov/Genbank/
    // European EMBL-Bank is here: http://www.ebi.ac.uk/embl/
    // Japanese DDBJ is here: http://www.ddbj.nig.ac.jp/
    // EMBL format http://www.ebi.ac.uk/embl/Documentation/User_manual/usrman.html

    virtual const char* GetMIMEType()
    { return "chemical/x-genbank"; }

    virtual unsigned int Flags()
    { return NOTWRITABLE | READONEONLY; }
    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

    enum SequenceType
      {
      UnknownSequence,
      ProteinSequence,
      DNASequence,
      RNASequence,
      MAXSequence
      };
  };

  //Make an instance of the format class
  GenBankFormat theGenBankFormat;

  // from formats/fastaformat.cpp
  bool ReadFASTASequence(OBMol * pmol, int seq_type, std::istream * in, bool create_bonds, bool bond_orders,
                         bool singleStrand, const char *turns = 0);
  /////////////////////////////////////////////////////////////////
  bool GenBankFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == 0)
      return false;


    std::istream * in = pConv->GetInStream();
    pmol->BeginModify();

    std::string line;

    SequenceType sequence_type = UnknownSequence;
    getline( * in, line);
    while (!in->eof())
      { // We have a new line in 'line'
      if (!line.compare(0, 6, "LOCUS", 6) || !line.compare(0, 2, "ID ", 2))
        {
        if (sequence_type == GenBankFormat::UnknownSequence)
          { // attempt to determine the sequence type
          if (line.find("RNA") != std::string::npos)
            sequence_type = GenBankFormat::RNASequence;
          else if (line.find("DNA") != std::string::npos)
            sequence_type = GenBankFormat::DNASequence;
          }
        }
      else if (!line.compare(0, 6, "DEFINITION", 6) || !line.compare(0, 2, "DE ", 2))
        {
        if (sequence_type == GenBankFormat::UnknownSequence)
          { // attempt to determine the sequence type
          if (line.find("RNA") != std::string::npos)
            sequence_type = GenBankFormat::RNASequence;
          else if (line.find("DNA") != std::string::npos)
            sequence_type = GenBankFormat::DNASequence;
          else if (line.find("gene") != std::string::npos)
            sequence_type = GenBankFormat::DNASequence;
          }
        if (pmol->GetTitle()[0] == 0)
          {
          std::string::size_type fc = line.find(' ');
          while (fc != std::string::npos && strchr(" \t\n\r", line[fc]))
            ++ fc;
          if (fc != std::string::npos)
            pmol->SetTitle( & (line.c_str()[fc]) );
          }
        }
      else if (!line.compare(0, 6, "ORIGIN", 6) || !line.compare(0, 2, "SQ ", 2))
        break;

      getline( * in, line);
      }
    if (sequence_type == GenBankFormat::UnknownSequence)
      sequence_type = GenBankFormat::DNASequence;

    bool rv = ReadFASTASequence(pmol, sequence_type, in,
        !pConv->IsOption("b",OBConversion::INOPTIONS), !pConv->IsOption("s",OBConversion::INOPTIONS), false);
	  pmol->EndModify();
    return rv;
  }

} //namespace OpenBabel
