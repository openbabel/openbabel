/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2003-2008 Geoffrey R. Hutchison
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
#include <openbabel/mol.h>

// Parse Atom records in PDB or PQR formats

namespace OpenBabel
{

  ////////////////////////////////////////////////////////////////
  bool ParseAtomRecord(char *buffer, OBMol &mol,int chainNum)
  /* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3)" */
  {
    string sbuf = &buffer[6];
    if (sbuf.size() < 48)
      return(false);

    bool hetatm = (EQn(buffer,"HETATM",6)) ? true : false;

    /* serial number */
    string serno = sbuf.substr(0,5);
    //SerialNum(the_atom) = atoi(tmp_str);

    /* atom name */
    string atmid = sbuf.substr(6,4);

    /* element */
    string element;
    if (sbuf.size() > 71)
      element = sbuf.substr(70,2);
    else
      element = "  ";

    //trim spaces on the right and left sides
    while (!atmid.empty() && atmid[0] == ' ')
      atmid = atmid.substr(1,atmid.size()-1);

    while (!atmid.empty() && atmid[atmid.size()-1] == ' ')
      atmid = atmid.substr(0,atmid.size()-1);

    /* residue name */

    string resname = sbuf.substr(11,3);
    if (resname == "   ")
      resname = "UNK";
    else
      {
        while (!resname.empty() && resname[0] == ' ')
          resname = resname.substr(1,resname.size()-1);

        while (!resname.empty() && resname[resname.size()-1] == ' ')
          resname = resname.substr(0,resname.size()-1);
      }

    string type;
    if (EQn(buffer,"ATOM",4))
      {
        type = atmid.substr(0,2);
        if (isdigit(type[0])) {
          // sometimes non-standard files have, e.g 11HH
          if (!isdigit(type[1])) type = atmid.substr(1,1);
          else type = atmid.substr(2,1); 
        } else if (sbuf[6] == ' ' &&
                   strncasecmp(type.c_str(), "Zn", 2) != 0 &&
                   strncasecmp(type.c_str(), "Fe", 2) != 0)
          type = atmid.substr(0,1);     // one-character element
        

        if (resname.substr(0,2) == "AS" || resname[0] == 'N')
          {
            if (atmid == "AD1")
              type = "O";
            if (atmid == "AD2")
              type = "N";
          }
        if (resname.substr(0,3) == "HIS" || resname[0] == 'H')
          {
            if (atmid == "AD1" || atmid == "AE2")
              type = "N";
            if (atmid == "AE1" || atmid == "AD2")
              type = "C";
          }
        if (resname.substr(0,2) == "GL" || resname[0] == 'Q')
          {
            if (atmid == "AE1")
              type = "O";
            if (atmid == "AE2")
              type = "N";
          }
      }
    else //must be hetatm record
      {
        if (isalpha(element[1]) && (isalpha(element[0]) || (element[0] == ' ')))
          {
            if (isalpha(element[0]))
              type = element.substr(0,2);
            else
              type = element.substr(1,1);
            if (type.size() == 2)
              type[1] = tolower(type[1]);
          }
        else
          {
            if (isalpha(atmid[0])) {
              
              if (atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' '))
                type = atmid.substr(0,2);
              else if (atmid[0] == 'A') // alpha prefix
                type = atmid.substr(1, atmid.size() - 1);
              else
                type = atmid.substr(0,1);
            }
            else if (atmid[0] == ' ')
              type = atmid.substr(1,1); // one char element
            else
              type = atmid.substr(1,2);

            // Some cleanup steps
            if (atmid == resname)
              {
                type = atmid;
                if (type.size() == 2)
                  type[1] = tolower(type[1]);
              }
            else
              if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
                  resname == "GPG" || resname == "NAD" || resname == "NAL" ||
                  resname == "NDP" || resname == "ABA")
                {
                  if (type.size() > 1)
                    type = type.substr(0,1);
                  //type.erase(1,type.size()-1);
                }
              else
                if (isdigit(type[0]))
                  {
                    type = type.substr(1,1);
                  }
                else
                  if (type.size() > 1 && isdigit(type[1]))
                    type = type.substr(0,1);
                  else
                    if (type.size() > 1 && isalpha(type[1])) {
                      if (type[0] == 'O' && type[1] == 'H')
                        type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
                      else if(isupper(type[1]))
                        {
                          type[1] = tolower(type[1]);
                        }
                    }
          }
        
      }

    OBAtom atom;
    /* X, Y, Z */
    string xstr = sbuf.substr(24,8);
    string ystr = sbuf.substr(32,8);
    string zstr = sbuf.substr(40,8);
    vector3 v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
    atom.SetVector(v);

    // useful for debugging unknown atom types (e.g., PR#1577238)
    //    cout << mol.NumAtoms() + 1 << " " << atmid << " type: " << type << endl;
    atom.SetAtomicNum(etab.GetAtomicNum(type.c_str()));

    /* residue sequence number */
    string resnum = sbuf.substr(16,4);
    OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
    if (res == NULL || res->GetName() != resname 
        || res->GetNumString() != resnum)
      {
        vector<OBResidue*>::iterator ri;
        for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
          if (res->GetName() == resname 
              && res->GetNumString() == resnum
              && static_cast<int>(res->GetChainNum()) == chainNum)
            break;

        if (res == NULL)
          {
            res = mol.NewResidue();
            res->SetChainNum(chainNum);
            res->SetName(resname);
            res->SetNum(resnum);
          }
      }

    if (!mol.AddAtom(atom))
      return(false);
    else
      {
        OBAtom *atom = mol.GetAtom(mol.NumAtoms());

        res->AddAtom(atom);
        res->SetSerialNum(atom, atoi(serno.c_str()));
        res->SetAtomID(atom, atmid);
        res->SetHetAtom(atom, hetatm);

        return(true);
      }
  }

} // end namespace Open Babel
