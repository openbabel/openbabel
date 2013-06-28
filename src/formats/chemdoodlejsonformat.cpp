/**********************************************************************
Copyright (C) 2013 by Matt Swain <m.swain@me.com>

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
#include <openbabel/obmolecformat.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>

#include "jsoncpp/json/json.h"
#include "jsoncpp/json/customwriter.h"

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
      "ChemDoodle JSON format is the native way to present data to \n"
      "the ChemDoodle Web Components.\n\n"
      
      "Write Options, e.g. -xd\n"
      " m  minified output formatting, with no line breaks or indents\n"
      " v  verbose output (include default values)\n"
      " w  use wedge/hash bonds from input instead of perceived stereochemistry\n\n"
      ;
    };
  
    virtual const char* SpecificationURL()
    { return "http://web.chemdoodle.com/docs/chemdoodle-json-format"; };
  
    virtual unsigned int Flags()
    { return ZEROATOMSOK; };

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    Json::Value inRoot;
    Json::Value outRoot;
    int currentMolIndex;
    
  };
  
  ChemDoodleJSONFormat theChemDoodleJSONFormat;
  
  bool ChemDoodleJSONFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == NULL) return false;
    istream& ifs = *pConv->GetInStream();
    
    if ( !ifs.good() || ifs.peek() == EOF )
      return false;
    
    map<OBBond*, OBStereo::BondDirection> updown;
    
    pmol->BeginModify();
    
    // Parse entire file into memory once, then reuse inRoot for subsequent molecules
    // (It's really tricky to stream json)
    if (inRoot.empty()) {
      Json::Reader reader;
      if (!reader.parse(ifs, inRoot)) {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", reader.getFormattedErrorMessages(), obError);
        return false;
      }
      ifs.clear();
      ifs.seekg(0, ios::beg);
      currentMolIndex = 0;
    }
    
    // Get the root level of the molecule
    Json::Value molRoot;
    if (!inRoot["m"].empty() && inRoot["m"].isArray()) {
      // File contains an array of molecules, iterate over them
      if (inRoot["m"].size() > currentMolIndex) {
        molRoot = inRoot["m"][currentMolIndex];
        currentMolIndex++;
      } else {
        // Finished the last molecule, reset inRoot in case multiple input files
        inRoot.clear();
        ifs.seekg(0, ios::end);
        return false;
      }
    } else if (!inRoot["a"].empty()) {
      // File is a single molecule, copy and reset inRoot in case multiple input files
      molRoot = Json::Value(inRoot);
      inRoot.clear();
      ifs.seekg(0, ios::end);
    } else {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "No molecules found", obError);
      return false;
    }
    
    // Set dimension to 2 unless z coordinates are found
    int dim = 2;
    
    // Atoms
    Json::Value atoms = molRoot["a"];
    if (!atoms.empty() && !atoms.isArray()) {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "Atoms must be specified in an array", obError);
      return false;
    }
    pmol->ReserveAtoms(atoms.size());
    for(Json::ArrayIndex i = 0; i < atoms.size(); i++) {
      Json::Value atom = atoms.get(i, 0);
      if (!atom.isObject()) {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid atom", obWarning);
        continue;
      }
      // Coordinates
      double x,y,z = 0;
      OBAtom* patom = pmol->NewAtom();
      if (atom["x"].isNumeric()) {
        x = atom["x"].asDouble();
      } else {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Atom missing x coordinate", obWarning);
      }
      if (atom["y"].isNumeric()) {
        y = atom["y"].asDouble();
      } else {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Atom missing y coordinate", obWarning);
      }
      if (atom["z"].isNumeric()) {
        z = atom["z"].asDouble();
        dim = 3;
      }
      patom->SetVector(x,y,z);
      // Element
      string symbol = "C";
      if (atom["l"].isString()) {
        symbol = atom["l"].asString();
      }
      if ((!atom["q"].empty() && atom["q"].asBool()) || (!atom["rg"].empty() && atom["rg"].asInt() != -1)) {
        // q and rg indicate an "any" atom or an rgroup
        patom->SetAtomicNum(0);
      } else {
        int isotope = 0;
        patom->SetAtomicNum(etab.GetAtomicNum(symbol, isotope));
        if (isotope != 0) {
          // isotope gets set for e.g. 'D' or 'T' symbol
          patom->SetIsotope(isotope);
        }
      }
      // Charge
      if (atom["c"].isIntegral()) {
        patom->SetFormalCharge(atom["c"].asInt());
      }
      // Mass
      if (atom["m"].isIntegral() && atom["m"].asInt() != -1) {
        patom->SetIsotope(atom["m"].asInt());
      }
      // Radicals
      if (atom["r"].isIntegral()) {
        int sm = 0;
        if (atom["r"].asInt() == 0) {
          sm = 0;
        } else if (atom["r"].asInt() == 1) {
          sm = 2;   // Monovalent
        } else if (atom["r"].asInt() == 2) {
          sm = 3;   // Divalent (Carbene, Nitrene)
        } else {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid atom radical", obWarning);
        }
        patom->SetSpinMultiplicity(sm);        
      }
      // Lone pairs
      if (atom["p"].isIntegral()) {
        OBPairInteger *lp = new OBPairInteger;
        lp->SetAttribute("p");
        lp->SetValue(atom["p"].asInt());
        lp->SetOrigin(fileformatInput);
        patom->SetData(lp);
      }
      // Atom identifier string
      if (atom["i"].isString()) {
        OBPairData *id = new OBPairData;
        id->SetAttribute("i");
        id->SetValue(atom["i"].asString());
        id->SetOrigin(fileformatInput);
        patom->SetData(id);
      }
      // R group
      if (atom["rg"].isIntegral()) {
        OBPairInteger *rg = new OBPairInteger;
        rg->SetAttribute("rg");
        rg->SetValue(atom["rg"].asInt());
        rg->SetOrigin(fileformatInput);
        patom->SetData(rg);
      }
      // "Any" atom
      if (atom["q"].isBool()) {
        OBPairBool *q = new OBPairBool;
        q->SetAttribute("q");
        q->SetValue(atom["q"].asBool());
        q->SetOrigin(fileformatInput);
        patom->SetData(q);
      }
    }
    pmol->SetDimension(dim);
    
    // Bonds
    Json::Value bonds = molRoot["b"];
    if (!bonds.empty() && !bonds.isArray()) {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bonds must be specified in an array", obError);
      return false;
    }
    for(Json::ArrayIndex i = 0; i < bonds.size(); i++) {
      Json::Value bond = bonds.get(i, 0);
      if (!bond.isObject()) {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond", obWarning);
        continue;
      }
      int begin, end, flag = 0;
      if (bond["b"].isIntegral()) {
        begin = bond["b"].asInt() + 1;  // +1 to index as OB index starts at 1
      } else {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond missing begin atom (\"b\")", obWarning);
        continue;
      }
      if (bond["e"].isIntegral()) {
        end = bond["e"].asInt() + 1;    // +1 to index as OB index starts at 1
      } else {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond missing end atom (\"e\")", obWarning);
        continue;
      }
      int order = 1;
      if (bond["o"].isNumeric()) {
        if (fabs(bond["o"].asDouble() - 1) < 0.01) {
          order = 1;
        } else if (fabs(bond["o"].asDouble() - 2) < 0.01) {
          order = 2;
        } else if (fabs(bond["o"].asDouble() - 3) < 0.01) {
          order = 3;
        } else if (fabs(bond["o"].asDouble() - 4) < 0.01) {
          order = 4;
        } else if (fabs(bond["o"].asDouble() - 0) < 0.01) {
          order = 0;
        } else if (fabs(bond["o"].asDouble() - 0.5) < 0.01) {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond order 0.5 not supported, using 1", obWarning);
          order = 1;
        } else if (fabs(bond["o"].asDouble() - 1.5) < 0.01) {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Bond order 1.5 not supported, using 2", obWarning);
          order = 2;
        } else {
          obErrorLog.ThrowError("ChemDoodleJSONFormat", "Unsupported bond order, skipping bond", obWarning);
          continue;
        }
      } else if (!bond["o"].empty()) {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond order", obWarning);
        continue;
      }
      if (begin == end || (unsigned)begin > pmol->NumAtoms() ||
          (unsigned)end > pmol->NumAtoms() || pmol->GetBond(begin, end)) {
        obErrorLog.ThrowError("ChemDoodleJSONFormat", "Invalid bond", obWarning);
        continue;
      }
      // Bond stereo
      string stereo = "none";
      if (bond["s"].isString()) {
        stereo = bond["s"].asString(); // 'none', 'protruding', 'recessed', 'ambiguous'
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
      OBBond* pbond = pmol->NewBond();
      pbond->Set(pbond->GetIdx(), pmol->GetAtom(begin), pmol->GetAtom(end), order, flag);
      pmol->GetAtom(begin)->AddBond(pbond);
      pmol->GetAtom(end)->AddBond(pbond);
      // Bond identifier string
      if (bond["i"].isString()) {
        OBPairData *id = new OBPairData;
        id->SetAttribute("i");
        id->SetValue(bond["i"].asString());
        id->SetOrigin(fileformatInput);
        pbond->SetData(id);
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
      map<OBBond*, OBStereo::BondDirection>::const_iterator bd_it;
      OpenBabel::OBStereoFacade facade(pmol);
      for(bd_it = updown.begin(); bd_it != updown.end(); ++bd_it) {
        OBBond* bond = bd_it->first;
        if (bond->GetBondOrder() != 2 || bd_it->second != OBStereo::UnknownDir)
          continue; // Only continue for those double bonds with UnknownDir
        OBCisTransStereo* ct = facade.GetCisTransStereo(bond->GetId());
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
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == NULL)
      return false;
    ostream& ofs = *pConv->GetOutStream();
    
    if (pmol->GetDimension() == 0) {
      obErrorLog.ThrowError("ChemDoodleJSONFormat", "No 2D or 3D coordinates exist. "
                 "To generate 2D or 3D coordinates use --gen2D or --gen3D.", obError);
      return false;
    }
    
    PerceiveStereo(pmol);
    
    // Kekulize any untyped aromatic bonds (5)
    FOR_BONDS_OF_MOL(bond, pmol) {
      if (bond->GetBondOrder() == 5) {
        pmol->Kekulize();
        break;
      }
    }
    
    // Set up all the stereochemistry information
    set<OBBond*> unspec_ctstereo = GetUnspecifiedCisTrans(*pmol);
    map<OBBond*, OBStereo::BondDirection> updown;
    map<OBBond*, OBStereo::Ref> from;
    map<OBBond*, OBStereo::Ref>::const_iterator from_cit;
    if (!pConv->IsOption("w", pConv->OUTOPTIONS))
        TetStereoToWedgeHash(*pmol, updown, from);
    
    // Construct the JSON objects
    Json::Value doc(Json::objectValue);
    Json::Value atoms(Json::arrayValue);
    Json::Value bonds(Json::arrayValue);
    
    // Atoms
    FOR_ATOMS_OF_MOL(patom, pmol) {
      Json::Value atom(Json::objectValue);
      // Coordinates
      // TODO: An option to round coordinates to n decimal places?
      atom["x"] = patom->GetX();
      atom["y"] = patom->GetY();
      if (pmol->GetDimension() == 3) {
        atom["z"] = patom->GetZ();
      }
      // Element
      if (patom->GetAtomicNum()) {
        if (patom->GetAtomicNum() != 6 || pConv->IsOption("v", pConv->OUTOPTIONS)) {
          atom["l"] = etab.GetSymbol(patom->GetAtomicNum());
        }
      } else {
        // No atomic number, could be an r group or a * atom
        if (patom->HasData("rg")) {
          OBPairInteger *rg = dynamic_cast<OBPairInteger*>(patom->GetData("rg"));
          atom["rg"] = rg->GetGenericValue();
        } else if (patom->HasData("q")) {
          OBPairBool *q = dynamic_cast<OBPairBool*>(patom->GetData("q"));
          atom["q"] = q->GetGenericValue();
        }
      }
      // Charge
      if (patom->GetFormalCharge() != 0 || pConv->IsOption("v", pConv->OUTOPTIONS)) {
        atom["c"] = patom->GetFormalCharge();
      }
      // Mass
      int m = patom->GetIsotope();
      if (m != 0 || pConv->IsOption("v", pConv->OUTOPTIONS)) {
        atom["m"] = (m == 0) ? -1 : m;
      }
      // Radicals
      int sm = patom->GetSpinMultiplicity();
      if (sm != 0 || pConv->IsOption("v", pConv->OUTOPTIONS)) {
        if (sm == 0) {
          atom["r"] = 0;
        } else if (sm == 2) {
          atom["r"] = 1;
        } else if (sm == 3 || sm == 1) {
          atom["r"] = 2;
        }
      }
      // Lone pairs
      if (patom->HasData("p")) {
        OBPairInteger *lp = dynamic_cast<OBPairInteger*>(patom->GetData("p"));
        atom["p"] = lp->GetGenericValue();
      } else if (pConv->IsOption("v", pConv->OUTOPTIONS)) {
        atom["p"] = 0;
      }
      // Atom identifier string
      if (patom->HasData("i")) {
        OBPairData *id = dynamic_cast<OBPairData*>(patom->GetData("i"));
        atom["i"] = id->GetValue();
      }
      atoms.append(atom);
    }
    
    // Bonds
    FOR_BONDS_OF_MOL(pbond, pmol) {
      Json::Value bond(Json::objectValue);
      
      bond["b"] = pbond->GetBeginAtomIdx()-1;
      bond["e"] = pbond->GetEndAtomIdx()-1;
      
      // Order
      int order = pbond->GetBondOrder();
      if (order != 1 || pConv->IsOption("v", pConv->OUTOPTIONS)) {
        bond["o"] = order;
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
          int tmp = bond["b"].asInt();
          bond["b"] = bond["e"];
          bond["e"] = tmp;
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
      if (stereo != "none" || pConv->IsOption("v", pConv->OUTOPTIONS)) {
        bond["s"] = stereo;
      }
      // Bond identifier string
      if (pbond->HasData("i")) {
        OBPairData *id = dynamic_cast<OBPairData*>(pbond->GetData("i"));
        bond["i"] = id->GetValue();
      }
      bonds.append(bond);
    }
    
    doc["a"] = atoms;
    doc["b"] = bonds;
    outRoot["m"].append(doc);
    
    if (pConv->IsLast()) {
      // Different JSON output formatting
      if (pConv->IsOption("m", pConv->OUTOPTIONS)) {
        Json::FastWriter fwriter;
        ofs << fwriter.write(outRoot);
      } else {
        Json::CustomWriter cwriter = Json::CustomWriter("{", "}", "[", "]", ": ", ",", "  ", 74);
        ofs << cwriter.write(outRoot);
      }
      outRoot.clear();  // Clear in case multiple output files
    }
    return true;
  }
  
}

