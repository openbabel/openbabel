/**********************************************************************
Copyright (C) 2006 by Sangwoo Shim
 
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

using namespace std;
namespace OpenBabel
{

  class FASTAFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID in the constructor
    FASTAFormat() {
      OBConversion::RegisterFormat("fasta", this, "chemical/x-fasta");
      OBConversion::RegisterFormat("fa", this);
      OBConversion::RegisterFormat("fsa", this);
      OBConversion::RegisterOptionParam("n", this);
    }

    virtual const char* Description() //required
    {
      return
        "FASTA format\n"
        "A file format used to exchange information between\n"
        "genetic sequence databases\n"
        "Write Options e.g. -xn \n"
        "  n  Omit title and comments\n";
    };

    virtual const char* SpecificationURL() {
      return "http://www.ebi.ac.uk/help/formats_frame.html";
    };

    virtual unsigned int Flags() {
      return NOTREADABLE | WRITEONEONLY;
    };

    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    string conv_3to1(string three);
  };  

  FASTAFormat theFASTAFormat;

  bool FASTAFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    string seq;
    OBMol* pmol;
    //   OBResidue *res;
	
    pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol == NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    FOR_RESIDUES_OF_MOL(res,pmol) {
      if (res->GetAtoms().size() > 3)
        seq.append(conv_3to1(res->GetName()));
    }
    if(!pConv->IsOption("n")) {
      if (strlen(pmol->GetTitle()) > 0)
        ofs << ">" << pmol->GetTitle();
      else
        ofs << ">Unknown molecule";
      ofs << " " << pmol->NumResidues() << " bp";
      ofs << "; generated with OpenBabel " << BABEL_VERSION << endl;
    }
    ofs << seq << endl;
    return true;
  }

  string
  FASTAFormat::conv_3to1(string three)
  {
    const char *aa_tbl[][3] = {
	    {"alanine", "ALA", "A"}, 
	    {"arginine", "ARG", "R"}, 
	    {"asparagine", "ASN", "N"}, 
	    {"aspartic acid", "ASP", "D"}, 
	    {"asparagine or aspartic acid", "ASX", "B"}, 
	    {"cysteine", "CYS", "C"}, 
	    {"glutamic acid", "GLU", "E"}, 
	    {"glutamine", "GLN", "Q"}, 
	    {"glutamine or glutamic acid", "GLX", "Z"}, 
	    {"glycine", "GLY", "G"}, 
	    {"histidine", "HIS", "H"}, 
	    {"isoleucine", "ILE", "I"}, 
	    {"leucine", "LEU", "L"}, 
	    {"lysine", "LYS", "K"}, 
	    {"methionine", "MET", "M"}, 
	    {"phenylalanine", "PHE", "F"}, 
	    {"proline", "PRO", "P"}, 
	    {"serine", "SER", "S"}, 
	    {"threonine", "THR", "T"}, 
	    {"tryptophan", "TRP", "W"}, 
	    {"tyrosine", "TYR", "Y"}, 
	    {"valine", "VAL", "V"}, 
	    {NULL, NULL, NULL}
    };
    int i;

    for (i = 0; aa_tbl[i][1] != NULL; i++) {
      if (strncasecmp(three.c_str(), aa_tbl[i][1], 3) == 0)
        return (string)aa_tbl[i][2];
    }
    return "X";
  }

} //namespace OpenBabel
