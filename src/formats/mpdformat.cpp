/***********************************************************************
mpdformat.cpp - Write only format to produce descriptors of molecules

Copyright (C) 2005 Nick England

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
// Output format is #Origatomtype;#layer-#frequency-#atomtype;#l-#f-#aty;...<tab>Next atom<newline>next molecule

#include <openbabel/babelconfig.h>

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/data.h>
#include <cstdlib>


#define LAYER_DEPTH 2 // cannot increase past 2 without adding more *nbr atom pointers and loops
#define LAYER_SIZE 184 // number of types needed for types system used
#define SEP_0 ";"  // separator between types
#define SEP_1 "-"  // separator for data layer-freq-type
#define SEP_2 '\t' // separator for atoms

using namespace std;
namespace OpenBabel
{
  class MPDFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MPDFormat()
    {
      OBConversion::RegisterFormat("mpd",this);
      OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("c", this);
      OBConversion::RegisterOptionParam("i", this);
    }

    virtual const char* Description() //required
    {
      return
        "MolPrint2D format\n"
        "An implementation of the circular fingerprint MolPrint2D\n"
        "MolPrint2D is an atom-environment fingerprint developed by Bender et al [bmg2004]_\n"
        "which has been used in QSAR studies and for measuring molecular similarity.\n\n"

        "The format of the output is as follows::\n\n"
        "   [Molec_name]\\t[atomtype];[layer]-[frequency]-[neighbour_type];\n\n"
        "Example for the SMILES string ``CC(=O)Cl``::\n\n"
        "   acid chloride   1;1-1-2;2-1-9;2-1-15;   2;1-1-1;1-1-9;1-1-15;\n"
        "                   9;1-1-2;2-1-1;2-1-15;   15;1-1-2;2-1-1;2-1-9;\n\n"

".. [bmg2004] Andreas Bender, Hamse Y. Mussa, and Robert C. Glen. **Molecular\n"
"             Similarity Searching Using Atom Environments, Information-Based\n"
"             Feature Selection, and a Naive Bayesian Classifier.**\n"
"             *J. Chem. Inf. Comput. Sci.* **2004**, *44*, 170-178.\n"
"             [`Link <https://doi.org/10.1021/ci034207y>`_]\n\n"

           " Write Options: e.g. -xnc\n"
           "  n prefix molecule names with name of file \n"
           "  c use XML style separators instead \n"
           "  i use IDX atom types of babel internal \n\n";
    };

    virtual const char* SpecificationURL()
    {
      return "https://doi.org/10.1021/ci034207y";
    }; //optional


    virtual unsigned int Flags() //Flags() can return be any the following combined by | or be omitted if none apply
    {                            // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
      return NOTREADABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    void ClearLayer(int a[][LAYER_SIZE]);
    void PrintLayer(int a[][LAYER_SIZE],ostream &ofs);
    void PrintXML(int layer_a[][LAYER_SIZE],ostream &ofs);
    int MyType(string a);
  };
  //***

  //Make an instance of the format class
  MPDFormat theMPDFormat;

  void MPDFormat::ClearLayer(int layer_a[][LAYER_SIZE])
  {
    for(int n=0;n<LAYER_DEPTH;n++)
      {
        for(int m=0;m<LAYER_SIZE;m++)
          {
            layer_a[n][m]=0;
          }
      }
  }

  void MPDFormat::PrintLayer(int layer_a[][LAYER_SIZE],ostream &ofs)
  {
    int freq=0;
    for(int n=0;n<LAYER_DEPTH;n++)
      {
        for(int m=0;m<LAYER_SIZE;m++)
          {
            freq=layer_a[n][m];
            if (freq == 0) continue;
            ofs << n+1 << SEP_1 << freq << SEP_1 << m << SEP_0;
            layer_a[n][m]=0;
          }
      }
    ofs << SEP_2;
  }
  void MPDFormat::PrintXML(int layer_a[][LAYER_SIZE],ostream &ofs)
  {
    int freq=0;
    string outType;
    for(int n=0;n<LAYER_DEPTH;n++)
      {
        for(int m=0;m<LAYER_SIZE;m++)
          {
            freq=layer_a[n][m];
            if (freq == 0) continue;
            ofs << "<layer depth=\"" << n+1 << "\" "
                << "frequency=\"" << freq <<"\" "<<"type=\""<< m <<"\"/>";
            layer_a[n][m]=0;
          }
      }
    ofs << "</atom>";
  }
  /*int MPDFormat::MyType(string a)
    {
    int o=0;
    if (strcmp("C.3",a.c_str())==0) o=1;
    else if(strcmp("C.2",a.c_str())==0) o=2;
    else if(strcmp( "C.1",a.c_str())==0) o=4;
    else if(strcmp( "C.ar",a.c_str())==0) o=3;
    else if(strcmp( "C.cat",a.c_str())==0) o=33;
    else if(strcmp( "N.3",a.c_str())==0) o=5;
    else if(strcmp( "N.2",a.c_str())==0) o=6;
    else if(strcmp( "N.1",a.c_str())==0) o=7;
    else if(strcmp( "N.ar",a.c_str())==0) o=11;
    else if(strcmp( "N.am",a.c_str())==0) o=28;
    else if(strcmp( "N.pl3",a.c_str())==0) o=19;
    else if(strcmp( "N.4",a.c_str())==0) o=31;
    else if(strcmp( "O.3",a.c_str())==0) o=8;
    else if(strcmp( "O.2",a.c_str())==0) o=9;
    else if(strcmp( "O.co2",a.c_str())==0) o=32;
    else if(strcmp( "O.spc",a.c_str())==0) o=8;
    else if(strcmp( "O.t3p",a.c_str())==0) o=8;
    else if(strcmp( "S.3",a.c_str())==0) o=10;
    else if(strcmp( "S.2",a.c_str())==0) o=18;
    else if(strcmp( "S.o",a.c_str())==0) o=29;
    else if(strcmp( "S.o2",a.c_str())==0) o=30;
    else if(strcmp( "P.3",a.c_str())==0) o=12;
    else if(strcmp( "H",a.c_str())==0) o=13;
    else if(strcmp( "H.spc",a.c_str())==0) o=13;
    else if(strcmp( "H.t3p",a.c_str())==0) o=13;
    else if(strcmp( "F",a.c_str())==0) o=16;
    else if(strcmp( "Cl",a.c_str())==0) o=15;
    else if(strcmp( "Br",a.c_str())==0) o=14;
    else if(strcmp( "I",a.c_str())==0) o=17;
    else if(strcmp( "Si",a.c_str())==0) o=27;
    else if(strcmp( "LP",a.c_str())==0) o=20;
    else if(strcmp( "Du",a.c_str())==0) o=26;
    else if(strcmp( "Na",a.c_str())==0) o=21;
    else if(strcmp( "K",a.c_str())==0) o=22;
    else if(strcmp( "Ca",a.c_str())==0) o=23;
    else if(strcmp( "Li",a.c_str())==0) o=24;
    else if(strcmp( "Al",a.c_str())==0) o=25;
    else o=26;
    return (o);
    }
  */


  ///////////////////////////////////////////////////
  /* Now the Write molecule code */

  bool MPDFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;



    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    OBAtom *atom,*nbr,*nbr2; // define atom and neghbour atom pointers
    string str,src,name;     // str used for output, src for handling
    unsigned int orig,otyp;  // orig holds first index for removal from layer 2, otype for output
    //    char buffer[BUFF_SIZE];
    bool xml_true=false, pre_true=false, idx_true=false;
    ttab.SetFromType("INT");
    ttab.SetToType("SBN");
    int layer[LAYER_DEPTH][LAYER_SIZE]; // layer stores the frequencies of each atom type
    ClearLayer(layer);

		if(pConv->IsOption("n")) // appending file name to molecule names
      {
        name = pConv->GetInFilename();     // string name holds the filename for appending
        unsigned int dotpos=name.find(".");         // removes the extension(s) from the filename
        if (dotpos < name.length())name.erase(dotpos);
        pre_true = true;
      }

		if(pConv->IsOption("c")) // outputting in XML format
			xml_true=true;

		if(pConv->IsOption("i")) // using IDX not SBN
      {
        idx_true=true;
        ttab.SetToType("IDX");
      }

    str = mol.GetTitle();
    if(xml_true==true) // <xml>
      {
        ofs << "<molecule id=\"";
        if(pre_true==true)ofs << name;
        if (str.empty())
          {
            ofs << pConv->GetOutputIndex() << "\">";
          }
        else ofs << str << pConv->GetOutputIndex() << "\">";
      } // </xml>
    else{
      if (str.empty())
        {
          if (pre_true==true) {ofs << name << SEP_1;}
          ofs << "***" << pConv->GetOutputIndex()<< SEP_2;
        }
      else
        { if (pre_true==true){ofs << name << SEP_1;}
        ofs << str << SEP_2;
        }
    }
    vector<OBAtom*>::iterator i; // iterate over all atoms
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        src = atom->GetType();
        ttab.Translate(str,src);
        // if (idx_true==true){
        otyp = atoi(str.c_str());
        //}
        //  else {otyp=MyType(str);}
        orig = atom->GetIdx();
        if(xml_true==true){ ofs << "<atom type=\"" << otyp << "\">";}
        else ofs << otyp << SEP_0;

        vector<OBBond*>::iterator j; // iterate over its neighbours
        for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
          {
            src = nbr->GetType();
            ttab.Translate(str,src);
            // if (idx_true==true){
            otyp = atoi(str.c_str());
            //}
            //  else {otyp=MyType(str);}
            layer[0][otyp]=layer[0][otyp]+1;

            vector<OBBond*>::iterator k; // iterate again over neighbours
            for (nbr2 = nbr->BeginNbrAtom(k);nbr2;nbr2 = nbr->NextNbrAtom(k))
              {
                if (nbr2->GetIdx()==orig) continue;
                src = nbr2->GetType();
                ttab.Translate(str,src);
                // if (idx_true==true){
                otyp = atoi(str.c_str());
                //}
                //  else {otyp=MyType(str);}
                layer[1][otyp]=layer[1][otyp]+1;
              } // end k
          } // end j
        if(xml_true==true)PrintXML(layer,ofs);
        else PrintLayer(layer,ofs);
        //ClearLayer(layer);
      } // end i
    if(xml_true==true)ofs << "</molecule>";
    ofs << endl;
    return(true);
  } // writemolecule

} // namespace openbabel
