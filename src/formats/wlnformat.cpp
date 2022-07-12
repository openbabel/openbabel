/**********************************************************************
Copyright (C) 2019 by NextMove Software
Some portions Copyright (C) 2022 by Michael Blakey

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
#include <string>
#include <openbabel/base.h>
#include <openbabel/mol.h>

bool NMReadWLN(const char *ptr, OpenBabel::OBMol* mol);
bool MBWriterWLN(OpenBabel::OBMol *mol, std::string &buffer);

using namespace std;
namespace OpenBabel
{

  class WLNFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    WLNFormat()
    {
      OBConversion::RegisterFormat("wln", this);
    }

    virtual const char* Description() //required
    {
      return
        "Wiswesser Line Notation\n"
	    "A chemical line notation developed by Wiswesser\n\n"

        "WLN was invented in 1949, by William J. Wiswesser, as one of the first attempts\n"
        "to codify chemical structure as a line notation, enabling collation on punched\n"
        "cards using automatic tabulating machines and early electronic computers. WLN\n"
        "was a forerunner to the SMILES notation used in modern cheminformatics systems,\n"
        "which attempted to simplify the complex rules used in WLN encoding (at the\n"
        "expense of brevity) to come up with an algorithmic system more suitable for\n"
        "implementation on computers, where historically WLN was typically encoded\n"
        "by hand by trained registrars.\n\n"

        "WLN encoding makes use of uppercase letters, digits, spaces and punctuation:\n\n"

        "- E       Bromine atom\n"
        "- F       Fluorine atom\n"
        "- G       Chlorine atom\n"
        "- H       Hydrogen atom\n"
        "- I       Iodine atom\n"
        "- Q       Hydroxyl group, -OH\n"
        "- R       Benzene ring\n"
        "- S       Sulfur atom\n"
        "- U       Double bond\n"
        "- UU      Triple bond\n"
        "- V       Carbonyl, -C(=O)-\n"
        "- C       Unbranched carbon multiply bonded to non-carbon atom\n"
        "- K       Nitrogen atom bonded to more than three other atoms\n"
        "- L       First symbol of a carbocyclic ring notation\n"
        "- M       Imino or imido -NH-group\n"
        "- N       Nitrogen atom, hydrogen free, bonded to fewer than 4 atoms\n"
        "- O       Oxygen atom, hydrogen-free\n"
        "- T       First symbol of a heterocyclic ring notation\n"
        "- W       Non-linear dioxo group, as in -NO2 or -SO2-\n"
        "- X       Carbon attached to four atoms other than hydrogen\n"
        "- Y       Carbon attached to three atoms other then hydrogen\n"
        "- Z       Amino and amido NH2 group\n"
        "- <digit> Digits '1' to '9' denote unbranched alkyl chains\n"
        "- &       Sidechain terminator or, after a space, a component separator\n\n"

        "For a more complete description of the grammar see Smith's book [1], which more\n"
        "accurately reflects the WLN commonly encountered than Wiswesser's book [2].\n"
        "Additional WLN dialects include inorganic salts, and methyl contractions.\n\n"

        "Here are some examples of WLN strings along with a corresponding SMILES string:\n\n"

        "- WN3        [O-][N+](=O)CCC\n"
        "- G1UU1G     ClC#CCl\n"
        "- VH3        O=CCCC\n"
        "- NCCN       N#CC#N\n"
        "- ZYZUM      NC(=N)N\n"
        "- QY         CC(C)O\n"
        "- OV1 &-NA-  CC(=O)[O-].[Na+]\n"
        "- RM1R       c1ccccc1NCc2ccccc2\n"
        "- T56 BMJ B D - DT6N CNJ BMR BO1 DN1 & 2N1 & 1 EMV1U1   (osimertinib)\n"
        "  Cn1cc(c2c1cccc2)c3ccnc(n3)Nc4cc(c(cc4OC)N(C)CCN(C)C)NC(=O)C=C\n\n"

        "The reader was contributed by Roger Sayle (NextMove Software). The text of\n"
        "this description was taken from his Bio-IT World poster [3]. Note that not\n"
        "all of WLN is currently supported; initially 76% of the WLN strings\n"
        "found in PubChem can be interpreted. Contributions by Michael Blakey \n"
        "has added functionality for poly/peri fused rings, taking this percentage to"
        "97%, with a paper currently in review discussing read/write rules\n\n"

        "The Writer was contributed by Michael Blakey (University of Southampton). Not\n"
        "all of WLN can be converted at this time, with limited support for complex ring\n"
        "systems containing multiple peri and poly fused cycles. This still in active\n"
        "development, publication is currently in review stage.\n\n"

        "1. Elbert G. Smith, \"The Wiswesser Line-Formula Chemical Notation\",\n"
        "   McGraw-Hill Book Company publishers, 1968.\n"
        "2. William J. Wiswesser, \"A Line-Formula Chemical Notation\", Thomas Crowell\n"
        "   Company publishers, 1954.\n"
        "3. Roger Sayle, Noel O'Boyle, Greg Landrum, Roman Affentranger. \"Open\n"
        "   sourcing a Wiswesser Line Notation (WLN) parser to facilitate electronic\n"
        "   lab notebook (ELN) record transfer using the Pistoia Alliance's UDM\n"
        "   (Unified Data Model) standard.\" BioIT World. Apr 2019.\n"
        "   https://www.nextmovesoftware.com/posters/Sayle_WisswesserLineNotation_BioIT_201904.pdf\n"
        ;
    };

    // Now writable so does the flag go?
    //virtual unsigned int Flags()
    //{
      //return NOTWRITABLE;
    //}

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
    virtual const char* TargetClassDescription(){return OBMol::ClassDescription();};


  };
  //***

  //Make an instance of the format class
  WLNFormat theWLNFormat;

  /////////////////////////////////////////////////////////////////
  bool WLNFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == nullptr)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    const char* title = pConv->GetTitle();
    char buffer[BUFF_SIZE];

    if (!ifs.getline(buffer,BUFF_SIZE))
      return false;

    NMReadWLN(buffer, pmol);

    return true;
  }

  bool WLNFormat::WriteMolecule(OBBase *pOb, OBConversion *pConv) {
      // Following formats given for Smiles and Inchi
      OBMol* pmol = dynamic_cast<OBMol*>(pOb);
      if (pmol == nullptr)
          return false;


      ostream &ofs = *pConv->GetOutStream();

      if(!pConv->IsOption("n")) //OBConversion::OUTOPTIONS is the default
          ofs << "Title = " << pmol->GetTitle() << endl;

      std::string buffer;
      buffer.reserve(1000);

      if (!MBWriterWLN(pmol, buffer))
          return false;

      // I think this should do it?
      ofs << buffer;
      return true;
  }

  ////////////////////////////////////////////////////////////////

} //namespace OpenBabel
