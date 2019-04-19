/**********************************************************************
Copyright (C) 2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004-2008 by Chris Morley

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
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/internalcoord.h>
#include <openbabel/generic.h>

#include <cstdlib>

//Possible replacement for strcasestr. See end of file
const char *_strcasestr(const char *s, const char *pattern);

#ifndef strcasestr
#define strcasestr _strcasestr
#endif

using namespace std;
namespace OpenBabel
{

  class GaussianZMatrixInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GaussianZMatrixInputFormat()
    {
      OBConversion::RegisterFormat("gzmat",this, "chemical/x-gaussian-input");
    }

    virtual const char* Description() //required
    {
      return
        "Gaussian Z-Matrix Input\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n"
        "Write Options e.g. -xk\n"
        "  k  \"keywords\" Use the specified keywords for input\n"
        "  f    <file>     Read the file specified for input keywords\n\n";
    };

    virtual const char* SpecificationURL()
    { return "https://www.gaussian.com/zmat/"; };

    virtual const char* GetMIMEType()
    { return "chemical/x-gaussian-input"; };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  GaussianZMatrixInputFormat theGaussianZMatrixInputFormat;

  ////////////////////////////////////////////////////////////////

  bool GaussianZMatrixInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordsEnable = pConv->IsOption("k",OBConversion::GENOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);
    string defaultKeywords = "!Put Keywords Here, check Charge and Multiplicity.\n#";

    if(keywords) {
      defaultKeywords = keywords;
    }

    if (keywordsEnable) {
      string model;
      string basis;
      string method;

      OBPairData *pd = (OBPairData *) pmol->GetData("model");
      if(pd)
        model = pd->GetValue();

      pd = (OBPairData *) pmol->GetData("basis");
      if(pd)
        basis = pd->GetValue();

      pd = (OBPairData *) pmol->GetData("method");
      if(pd)
        method = pd->GetValue();

      if(method == "optimize") {
        method = "opt";
      }

      if(model != "" && basis != "" && method != "") {
        ofs << model << "/" << basis << "," << method << "\n";
      }
      else {
        ofs << "#Unable to translate keywords!\n";
        ofs << defaultKeywords << "\n";
      }
    }
    else if (keywordFile) {
      ifstream kfstream(keywordFile);
      string keyBuffer;
      if (kfstream) {
        while (getline(kfstream, keyBuffer))
          ofs << keyBuffer << "\n";
      }
    }
    else {
      ofs << defaultKeywords << "\n";
    }
    ofs << "\n"; // blank line after keywords
    ofs << " " << mol.GetTitle() << "\n\n";

    snprintf(buffer, BUFF_SIZE, "%d  %d",
             mol.GetTotalCharge(),
             mol.GetTotalSpinMultiplicity());
    ofs << buffer << "\n";

		// OK, now that we set the keywords, charge, etc. we generate the z-matrix

    OBAtom *a,*b,*c;//, *atom;

    vector<OBInternalCoord*> vic;
    vic.push_back((OBInternalCoord*)NULL);
		FOR_ATOMS_OF_MOL(atom, mol)
      vic.push_back(new OBInternalCoord);

    CartesianToInternal(vic,mol);

    double r,w,t;
		string type;
		FOR_ATOMS_OF_MOL(atom, mol)
      {
        a = vic[atom->GetIdx()]->_a;
        b = vic[atom->GetIdx()]->_b;
        c = vic[atom->GetIdx()]->_c;

				type = OBElements::GetSymbol(atom->GetAtomicNum());
				if (atom->GetIsotope() != 0) {
					snprintf(buffer, BUFF_SIZE, "(Iso=%d)", atom->GetIsotope());
					type += buffer;
				}

				switch(atom->GetIdx())
          {
					case 1:
            snprintf(buffer, BUFF_SIZE, "%-s\n",type.c_str());
            ofs << buffer;
            continue;
            break;

					case 2:
            snprintf(buffer, BUFF_SIZE, "%-s  %d  r%d\n",
                     type.c_str(), a->GetIdx(), atom->GetIdx());
            ofs << buffer;
            continue;
            break;

					case 3:
            snprintf(buffer, BUFF_SIZE, "%-s  %d  r%d  %d  a%d\n",
                     type.c_str(), a->GetIdx(), atom->GetIdx(), b->GetIdx(), atom->GetIdx());
            ofs << buffer;
            continue;
            break;

					default:
            snprintf(buffer, BUFF_SIZE, "%-s  %d  r%d  %d  a%d  %d  d%d\n",
                     type.c_str(), a->GetIdx(), atom->GetIdx(),
                     b->GetIdx(), atom->GetIdx(), c->GetIdx(), atom->GetIdx());
            ofs << buffer;
          }
      }

		ofs << "Variables:\n";

		FOR_ATOMS_OF_MOL(atom, mol)
      {
        r = vic[atom->GetIdx()]->_dst;
        w = vic[atom->GetIdx()]->_ang;
				if (w < 0.0)
					w += 360.0;
        t = vic[atom->GetIdx()]->_tor;
				if (t < 0.0)
					t += 360.0;

				switch(atom->GetIdx())
          {
					case 1:
            continue;
            break;

					case 2:
            snprintf(buffer, BUFF_SIZE, "r2= %6.4f\n", r);
            ofs << buffer;
            continue;
            break;

					case 3:
            snprintf(buffer, BUFF_SIZE, "r3= %6.4f\na3= %6.2f\n", r, w);
            ofs << buffer;
            continue;
            break;

					default:
            snprintf(buffer, BUFF_SIZE, "r%d= %6.4f\na%d= %6.2f\nd%d= %6.2f\n",
                     atom->GetIdx(), r, atom->GetIdx(), w, atom->GetIdx(), t);
            ofs << buffer;
          }
      }

    // file should end with a blank line
    ofs << "\n";
    return(true);
  }

	// This is based on the Babel 1.6 code. It may not be ideal.
  // If you have problems (or examples of input files that do not work), please contact
  // the openbabel-discuss@lists.sourceforge.net mailing list and/or post a bug
  bool GaussianZMatrixInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];

    OBAtom *atom;
		vector<OBInternalCoord*> vic;
	  vic.push_back((OBInternalCoord*)NULL); // OBMol indexed from 1 -- potential atom index problem

    vector<string> vs;
    int charge = 0;
    unsigned int spin = 1;
		unsigned int blankLines = 0;
		bool readVariables = false; // are we reading variables yet?
    bool foundVariables = false; // have we seen actual text for a variable line?
		map<string, double> variables; // map from variable name to value
    vector<string> atomLines; // save atom lines to re-parse after reading variables

    mol.BeginModify();

    //@todo: Read keywords from this line
    while (ifs.getline(buffer,BUFF_SIZE))
      if (strncmp(buffer,"#", 1) == 0) // begins with '#'
			  break;

		while (ifs.getline(buffer,BUFF_SIZE)) { // blank line, title, blank, then we begin
			if (strlen(buffer) == 0)
				blankLines++;
			if (blankLines == 2)
				break;
		}

    ifs.getline(buffer,BUFF_SIZE); // Charge = # Multiplicty = #
    tokenize(vs, buffer, " \t\n");
    if (vs.size() == 2)
      {
        charge = atoi(vs[0].c_str());
        spin = atoi(vs[1].c_str());
      }

		// We read through atom lines and cache them (into atomLines)
		while (ifs.getline(buffer,BUFF_SIZE)) {
			if (strlen(buffer) == 0) {
        if (foundVariables)
          break; // blank line after variable section
        else {
          readVariables = true;
          continue;
        }
      }

			if (strcasestr(buffer, "VARIABLE") != NULL) {
				readVariables = true;
        continue;
			}

			if (readVariables) {
			  tokenize(vs, buffer, "= \t\n");
			  if (vs.size() >= 2) {
          variables[vs[0]] = atof(vs[1].c_str());
          foundVariables = true;
          //          cerr << "var: " << vs[0] << " " << vs[1] << endl;
        }
			}
			else {
        atomLines.push_back(buffer);
        vic.push_back(new OBInternalCoord);
      }
		} // end while

    char *endptr;
    double temp;
		if (atomLines.size() > 0) {
      unsigned int i, j;
      for (i = 0; i < atomLines.size(); ++i) {
        j = i+1;
        tokenize(vs, atomLines[i]);
        atom = mol.NewAtom();
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));

        if (j == 1) {
          continue; // first atom, just create it
        }

        if (j >= 2) {
          if (vs.size() < 3) {return false;}
          vic[j]->_a = mol.GetAtom(atoi(vs[1].c_str()));

          temp = strtod((char*)vs[2].c_str(), &endptr);
          if (endptr != (char*)vs[2].c_str())
            vic[j]->_dst = temp;
          else
            vic[j]->_dst = variables[vs[2].c_str()];
        }

        if (j >= 3) {
          if (vs.size() < 5) {return false;}
          vic[j]->_b = mol.GetAtom(atoi(vs[3].c_str()));

          temp = strtod((char*)vs[4].c_str(), &endptr);
          if (endptr != (char*)vs[4].c_str())
            vic[j]->_ang = temp;
          else
            vic[j]->_ang = variables[vs[4].c_str()];
        }

        if (j >= 4) {
          if (vs.size() < 7) {return false;}
          vic[j]->_c = mol.GetAtom(atoi(vs[5].c_str()));

          temp = strtod((char*)vs[6].c_str(), &endptr);
          if (endptr != (char*)vs[6].c_str()) {
            vic[j]->_tor = temp;
          } else {
            const char* tor_str = vs[6].c_str();
            if (tor_str[0] == '-')
              vic[j]->_tor = -1 * variables[tor_str+1];
            else
              vic[j]->_tor = variables[tor_str];
          }
        }
      }
		}

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

		InternalToCartesian(vic,mol);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    mol.SetTotalCharge(charge);
    mol.SetTotalSpinMultiplicity(spin);

    mol.SetTitle(title);
    return(true);
  }

} //namespace OpenBabel

//*****************************************************************************
// Adapted from the ht://Dig package   <http://www.htdig.org/>
// Copyright (c) 1999, 2000, 2001 The ht://Dig Group
// For copyright details, see the file COPYING in your distribution
// or the GNU General Public License version 2 or later
// <http://www.gnu.org/copyleft/gpl.html>

const char *
_strcasestr(const char *s, const char *pattern)
{
    int		length = strlen(pattern);

    while (*s)
    {
	if (strncasecmp(s, pattern, length) == 0)
	    return s;
	s++;
    }
    return 0;
}

