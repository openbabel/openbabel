/**********************************************************************
Copyright (C) 2013, 2018 by Matt Swain <m.swain@me.com>

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

#include <map>
#include <openbabel/babelconfig.h>
#include <openbabel/json.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>

using namespace std;
namespace OpenBabel
{

class ChemDoodleJSONFormat : public OBMoleculeFormat
{
  public:
    ChemDoodleJSONFormat()
    {
      OBConversion::RegisterFormat("cdjson", this);
    }

    virtual const char* Description()
    {
      return
      "ChemDoodle JSON\n"
      "The native way to present data to the ChemDoodle Web Components\n\n"

      "Read Options, e.g. -ac\n"
      " c<num>  coordinate multiplier (default: 20)\n\n"

      "Write Options, e.g. -xc\n"
      " c<num>  coordinate multiplier (default: 20)\n"
      " m  minified output formatting, with no line breaks or indents\n"
      " v  verbose output (include default values)\n"
      " w  use wedge/hash bonds from input instead of perceived stereochemistry\n\n"
      ;
    };

    virtual const char* SpecificationURL()
    { return "http://web.chemdoodle.com/docs/chemdoodle-json-format"; };

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    rapidjson::Document inRoot;
    rapidjson::Document outRoot;
    int currentMolIndex;

  };
  
  ChemDoodleJSONFormat theChemDoodleJSONFormat;

  bool ChemDoodleJSONFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol *pmol = pOb->CastAndClear<OBMol>();
    if (pmol == NULL) return false;
    istream &ifs = *pConv->GetInStream();

    if (!ifs.good())
      return false;

    map<OBBond *, OBStereo::BondDirection> updown;

    pmol->BeginModify();

    // Parse entire file into memory once, then reuse inRoot for subsequent molecules
    // (It's really tricky to stream json)
    if (ifs.peek() != EOF) {
      rapidjson::IStreamWrapper isw(ifs);
      inRoot.ParseStream(isw);
      if (inRoot.HasParseError()) {
        stringstream msg;
        msg << "JSON parse error at offset " << inRoot.GetErrorOffset() << ": "
            << rapidjson::GetParseError_En(inRoot.GetParseError());
        obErrorLog.ThrowError("ChemDoodleJSONFormat", msg.str(), obError);
        return false;
      }
      // Clear ifs flags so it is "good" and any subsequent mols are read, but leave it at EOF position.
      // Therefore, when not at EOF position we know to parse next file and reset currentMolIndex
      ifs.clear();
      currentMolIndex = 0;
    }

    if (!inRoot.IsObject()) {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "JSON file should be a single object", obError);
      return false;
    }

    // Get the root level of the molecule
    rapidjson::Value molRoot;
    if (inRoot.HasMember("m") && inRoot["m"].IsArray()) {
      // File contains an array of molecules, iterate over them
      if (inRoot["m"].Size() > currentMolIndex) {
        molRoot = inRoot["m"][currentMolIndex];
        currentMolIndex++;
      } else {
        // Finished the last molecule
        ifs.seekg(0, ios::end);
        return false;
      }
    } else {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "JSON file must contain an 'm' molecules array", obError);
      return false;
    }
    
    // Set dimension to 2 unless z coordinates are found
    unsigned short dim = 2;
    
    // Atoms
    if (!molRoot.HasMember("a") || !molRoot["a"].IsArray()) {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "Atoms must be specified in an array", obError);
      return false;
    }
    const rapidjson::Value &atoms = molRoot["a"];
    pmol->ReserveAtoms(atoms.Size());
    double c = 20;
    if (pConv->IsOption("c")) {
      c = atof(pConv->IsOption("c"));
    }
    for (rapidjson::SizeType i = 0; i < atoms.Size(); i++) {
      const rapidjson::Value &atom = atoms[i];
      if (!atom.IsObject()) {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid atom", obWarning);
        continue;
      }
      // Coordinates
      double x = 0, y = 0, z = 0;
      OBAtom *patom = pmol->NewAtom();
      if (atom.HasMember("x") && atom["x"].IsNumber()) {
        x = atom["x"].GetDouble() / c;
      } else {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Atom missing x coordinate", obWarning);
      }
      if (atom.HasMember("y") && atom["y"].IsNumber()) {
        y = atom["y"].GetDouble() / c;
      } else {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Atom missing y coordinate", obWarning);
      }
      if (atom.HasMember("z") && atom["z"].IsNumber()) {
        z = atom["z"].GetDouble() / c;
        dim = 3;
      }
      patom->SetVector(x, y, z);
      // Element
      string symbol = "C";
      if (atom.HasMember("l") && atom["l"].IsString()) {
        symbol = atom["l"].GetString();
      }
      if (atom.HasMember("q")) {
        // q indicates a query atom
        patom->SetAtomicNum(0);
        // TODO: Parse query?
      } else {
        patom->SetAtomicNum(OBElements::GetAtomicNum(symbol.c_str()));
      }
      // Charge
      if (atom.HasMember("c") && atom["c"].IsInt()) {
        patom->SetFormalCharge(atom["c"].GetInt());
      }
      // Mass
      if (atom.HasMember("m") && atom["m"].IsInt()) {
        int mass = atom["m"].GetInt();
        if (mass > -1) {
          patom->SetIsotope((unsigned int) mass);
        }
      }
      // Radicals
      if (atom.HasMember("r") && atom["r"].IsInt()) {
        unsigned int rad = atom["r"].GetUint();
        unsigned short sm = 0;
        if (rad == 0) {
          sm = 0;
        } else if (rad == 1) {
          sm = 2;   // Monovalent
        } else if (rad == 2) {
          sm = 3;   // Divalent (Carbene, Nitrene)
        } else {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid atom radical", obWarning);
        }
        patom->SetSpinMultiplicity(sm);
      }
      // Hydrogen count
      if (atom.HasMember("h") && atom["h"].IsInt()) {
        int hcount = atom["h"].GetInt();
        if (hcount > -1) {
          patom->SetImplicitHCount((unsigned int) hcount);
        }
      }
      // Lone pairs
      if (atom.HasMember("p") && atom["p"].IsInt()) {
        OBPairInteger *lp = new OBPairInteger;
        lp->SetAttribute("p");
        lp->SetValue(atom["p"].GetInt());
        lp->SetOrigin(fileformatInput);
        patom->SetData(lp);
      }
      // Atom identifier string
      if (atom.HasMember("i") && atom["i"].IsString()) {
        OBPairData *id = new OBPairData;
        id->SetAttribute("i");
        id->SetValue(atom["i"].GetString());
        id->SetOrigin(fileformatInput);
        patom->SetData(id);
      }
    }

    pmol->SetDimension(dim);

    // Bonds
    if (molRoot.HasMember("b") && molRoot["b"].IsArray()) {
      const rapidjson::Value &bonds = molRoot["b"];
      for (rapidjson::SizeType i = 0; i < bonds.Size(); i++) {
        const rapidjson::Value &bond = bonds[i];
        if (!bond.IsObject()) {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond", obWarning);
          continue;
        }
        int begin = 0, end = 0, flag = 0;
        if (bond.HasMember("b") && bond["b"].IsInt()) {
          begin = bond["b"].GetInt() + 1;  // +1 to index as OB index starts at 1
        } else {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond missing begin atom (\"b\")", obWarning);
          continue;
        }
        if (bond.HasMember("e") && bond["e"].IsInt()) {
          end = bond["e"].GetInt() + 1;    // +1 to index as OB index starts at 1
        } else {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond missing end atom (\"e\")", obWarning);
          continue;
        }
        int order = 1;
        if (bond.HasMember("o") && bond["o"].IsNumber()) {
          double bo = bond["o"].GetDouble();
          if (fabs(bo - 1) < 0.01) {
            order = 1;
          } else if (fabs(bo - 2) < 0.01) {
            order = 2;
          } else if (fabs(bo - 3) < 0.01) {
            order = 3;
          } else if (fabs(bo - 4) < 0.01) {
            order = 4;
          } else if (fabs(bo - 0) < 0.01) {
            order = 0;
          } else if (fabs(bo - 0.5) < 0.01) {
            obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond order 0.5 not supported, using 1", obWarning);
            order = 1;
          } else if (fabs(bo - 1.5) < 0.01) {
            obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond order 1.5 not supported, using 2", obWarning);
            order = 2;
          } else {
            obErrorLog.ThrowError("ChemDoodleJSONFormat", "Unsupported bond order, skipping bond", obWarning);
            continue;
          }
        } else if (bond.HasMember("o")) {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond order", obWarning);
          continue;
        }
        if (begin == end || (unsigned) begin > pmol->NumAtoms() ||
            (unsigned) end > pmol->NumAtoms() || pmol->GetBond(begin, end)) {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond", obWarning);
          continue;
        }
        // Bond stereo
        string stereo = "none";
        if (bond.HasMember("s") && bond["s"].IsString()) {
          stereo = bond["s"].GetString(); // 'none', 'protruding', 'recessed', 'ambiguous'
        }
        if (stereo == "protruding") {
          flag |= OBBond::Wedge;
        } else if (stereo == "recessed") {
          flag |= OBBond::Hash;
        } else if (stereo == "ambiguous") {
          if (order == 1) {
            flag |= OBBond::WedgeOrHash;  // Single bond
          } else {
            flag |= OBBond::CisOrTrans;   // Double bond
          }
        } else if (stereo != "none") {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond stereo", obWarning);
        }
        OBBond *pbond = pmol->NewBond();
        pbond->Set(pbond->GetIdx(), pmol->GetAtom(begin), pmol->GetAtom(end), order, flag);
        pmol->GetAtom(begin)->AddBond(pbond);
        pmol->GetAtom(end)->AddBond(pbond);
        // Bond identifier string
        if (bond.HasMember("i") && bond["i"].IsString()) {
          OBPairData *id = new OBPairData;
          id->SetAttribute("i");
          id->SetValue(bond["i"].GetString());
          id->SetOrigin(fileformatInput);
          pbond->SetData(id);
        }
      }
    }

    // Set up the updown map we are going to use to derive stereo info
    FOR_BONDS_OF_MOL(pbond, pmol) {
      OBStereo::BondDirection bd = OBStereo::NotStereo;
      unsigned int flag = pbond->GetFlags();
      if (flag & OBBond::Wedge)
        bd = OBStereo::UpBond;
      if (flag & OBBond::Hash)
        bd = OBStereo::DownBond;
      if (flag & OBBond::WedgeOrHash)
        bd = OBStereo::UnknownDir;
      if (flag & OBBond::CisOrTrans && pbond->GetBondOrder() == 2)
        bd = OBStereo::UnknownDir;
      if (bd != OBStereo::NotStereo)
        updown[&*pbond] = bd;
    }

    // TODO: Do we need to do SetImplicitValence for each atom?

    // Automatically determine spin multiplicity for atoms with hydrogens specified
    pmol->AssignSpinMultiplicity();
    pmol->EndModify();

    if (pmol->Has3D()) {
      // Use 3D coordinates to determine stereochemistry
      StereoFrom3D(pmol);

      // For unspecified cis/trans stereos, set their Configs to unspecified
      map<OBBond *, OBStereo::BondDirection>::const_iterator bd_it;
      OpenBabel::OBStereoFacade facade(pmol);
      for (bd_it = updown.begin(); bd_it != updown.end(); ++bd_it) {
        OBBond *bond = bd_it->first;
        if (bond->GetBondOrder() != 2 || bd_it->second != OBStereo::UnknownDir)
          continue; // Only continue for those double bonds with UnknownDir
        OBCisTransStereo *ct = facade.GetCisTransStereo(bond->GetId());
        if (ct) {
          OBCisTransStereo::Config config = ct->GetConfig();
          config.specified = false;
          ct->SetConfig(config);
        }
      }
    } else if (pmol->Has2D()) {
      // Use 2D coordinates + hash/wedge to determine stereochemistry
      StereoFrom2D(pmol, &updown);
    }
    
    return true;
  }
  
  bool ChemDoodleJSONFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol *pmol = dynamic_cast<OBMol *>(pOb);
    if (pmol == NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    if (pmol->GetDimension() == 0) {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "No 2D or 3D coordinates exist. "
          "To generate 2D or 3D coordinates use --gen2D or --gen3D.", obError);
      return false;
    }

    PerceiveStereo(pmol);

    // Set up all the stereochemistry information
    set<OBBond *> unspec_ctstereo = GetUnspecifiedCisTrans(*pmol);
    map<OBBond *, OBStereo::BondDirection> updown;
    map<OBBond *, OBStereo::Ref> from;
    map<OBBond *, OBStereo::Ref>::const_iterator from_cit;
    if (!pConv->IsOption("w", pConv->OUTOPTIONS))
      TetStereoToWedgeHash(*pmol, updown, from);

    // Whether to include default values
    bool verbose = pConv->IsOption("v", pConv->OUTOPTIONS) != NULL;

    // Must always pass an allocator when memory may need to be allocated
    rapidjson::Document::AllocatorType &al = outRoot.GetAllocator();

    rapidjson::Value doc(rapidjson::kObjectType);  // Root of molecule JSON

    // Atoms
    rapidjson::Value atoms(rapidjson::kArrayType);
    double c = 20;
    if (pConv->IsOption("c")) {
      c = atof((const char *) pConv->IsOption("c"));
    }
    FOR_ATOMS_OF_MOL(patom, pmol) {
      rapidjson::Value atom(rapidjson::kObjectType);
      // Coordinates
      // TODO: An option to round coordinates to n decimal places?
      atom.AddMember("x", rapidjson::Value(patom->GetX() * c).Move(), al);
      atom.AddMember("y", rapidjson::Value(patom->GetY() * c).Move(), al);
      if (pmol->GetDimension() == 3) {
        atom.AddMember("z", rapidjson::Value(patom->GetZ() * c).Move(), al);
      }
      // Element
      if (patom->GetAtomicNum()) {
        if (patom->GetAtomicNum() != 6 || verbose) {
          atom.AddMember("l", rapidjson::Value(patom->GetAtomicNum()).Move(), al);
        }
      } else {
        // No atomic number, is a query atom
        atom.AddMember("q", rapidjson::Value(rapidjson::kObjectType).Move(), al);
      }
      // Charge
      if (patom->GetFormalCharge() != 0 || verbose) {
        atom.AddMember("c", rapidjson::Value(patom->GetFormalCharge()).Move(), al);
      }
      // Mass
      int m = patom->GetIsotope();
      if (m != 0 || verbose) {
        atom.AddMember("m", rapidjson::Value((m == 0) ? -1 : m).Move(), al);
      }
      // Radicals
      int sm = patom->GetSpinMultiplicity();
      if (sm != 0 || verbose) {
        if (sm == 0) {
          atom.AddMember("r", rapidjson::Value(0).Move(), al);
        } else if (sm == 2) {
          atom.AddMember("r", rapidjson::Value(1).Move(), al);
        } else if (sm == 3 || sm == 1) {
          atom.AddMember("r", rapidjson::Value(2).Move(), al);
        }
      }
      // Lone pairs
      if (patom->HasData("p")) {
        OBPairInteger *lp = dynamic_cast<OBPairInteger*>(patom->GetData("p"));
        atom.AddMember("p", rapidjson::Value(lp->GetGenericValue()).Move(), al);
      } else if (verbose) {
        atom.AddMember("p", rapidjson::Value(0).Move(), al);
      }
      // Atom identifier string
      if (patom->HasData("i")) {
        OBPairData *id = dynamic_cast<OBPairData *>(patom->GetData("i"));
        rapidjson::Value idValue(rapidjson::kStringType);
        idValue.SetString(id->GetValue().c_str(), al);
        atom.AddMember("i", idValue, al);
      }
      atoms.PushBack(atom, al);
    }
    doc.AddMember("a", atoms, al);


    // Bonds
    if (pmol->NumBonds() > 0) {
      rapidjson::Value bonds(rapidjson::kArrayType);
      FOR_BONDS_OF_MOL(pbond, pmol) {
        rapidjson::Value bond(rapidjson::kObjectType);

        bond.AddMember("b", rapidjson::Value(pbond->GetBeginAtomIdx() - 1).Move(), al);
        bond.AddMember("e", rapidjson::Value(pbond->GetEndAtomIdx() - 1).Move(), al);

        // Order
        int order = pbond->GetBondOrder();
        if (order != 1 || verbose) {
          bond.AddMember("o", rapidjson::Value(order).Move(), al);
        }

        // Stereochemistry
        string stereo = "none";
        if (pConv->IsOption("w", pConv->OUTOPTIONS)) {
          // Option w means just use input wedge/hash/ambiguous bonds
          if (pbond->IsWedge()) {
            stereo = "protruding";
          } else if (pbond->IsHash()) {
            stereo = "recessed";
          } else if (pbond->IsWedgeOrHash() || pbond->IsCisOrTrans()) {
            stereo = "ambiguous";
          }
        } else {
          // No option w means use calculated stereochemistry

          // Swap start and end atom if necessary
          from_cit = from.find(&*pbond);
          if (from_cit != from.end() && from_cit->second == pbond->GetEndAtom()->GetId()) {
            bond["b"].Swap(bond["e"]);
          }

          // Unspecified cis-trans stereo
          if (unspec_ctstereo.find(&*pbond) != unspec_ctstereo.end()) {
            stereo = "ambiguous";
          }

          if (updown.find(&*pbond) != updown.end()) {
            if (updown[&*pbond] == 1) {
              stereo = "protruding";
            } else if (updown[&*pbond] == 4) {
              stereo = "ambiguous";
            } else if (updown[&*pbond] == 6) {
              stereo = "recessed";
            }
          }
        }
        if (stereo != "none" || verbose) {
          rapidjson::Value stereoValue(rapidjson::kStringType);
          stereoValue.SetString(stereo.c_str(), al);
          bond.AddMember("s", stereoValue, al);
        }
        // Bond identifier string
        if (pbond->HasData("i")) {
          OBPairData *id = dynamic_cast<OBPairData *>(pbond->GetData("i"));
          rapidjson::Value idValue(rapidjson::kStringType);
          idValue.SetString(id->GetValue().c_str(), al);
          bond.AddMember("i", idValue, al);
        }
        bonds.PushBack(bond, al);
      }
      doc.AddMember("b", bonds, al);
    }

    // Create root object and m array if this is the first molecule in the file
    if (!outRoot.IsObject() || !outRoot.HasMember("m")) {
      outRoot.SetObject();
      outRoot.AddMember("m", rapidjson::Value(rapidjson::kArrayType).Move(), al);
    }

    // Add molecule to m array
    outRoot["m"].PushBack(doc, al);

    // Write json to output stream if this is the last molecule in the file
    if (pConv->IsLast()) {
      rapidjson::OStreamWrapper osw(ofs);
      if (pConv->IsOption("m", pConv->OUTOPTIONS)) {
        rapidjson::Writer<rapidjson::OStreamWrapper> writer(osw);
        outRoot.Accept(writer);
      } else {
        rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(osw);
        writer.SetIndent(' ', 2);
        outRoot.Accept(writer);
      }
      // Clear outRoot so it can be re-used
      rapidjson::Document(rapidjson::kObjectType).Swap(outRoot);
    }
    return true;
  }
  
}

