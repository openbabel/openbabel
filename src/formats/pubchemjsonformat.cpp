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
#include <algorithm>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/json/json.h>
#include <openbabel/json/customwriter.h>

using namespace std;
namespace OpenBabel
{

typedef OBPairTemplate< vector<string> > AnnotationData;

class PubChemJSONFormat : public OBMoleculeFormat
{
  public:
    PubChemJSONFormat()
    {
      OBConversion::RegisterFormat("pcjson",this);
    }

    virtual const char* Description()
    {
      return
      "PubChem JSON\n"
      "This is the format returned by the PubChem PUG REST service when\n"
      "requesting JSON. It closely resembles PubChem's internal data structure.\n\n"
      
      "Read Options, e.g. -as\n"
      " s  disable stereo perception and just read stereo information from input\n\n"
      
      "Write Options, e.g. -xm\n"
      " m  minified output, with no line breaks or indents\n"
      " w  use bond styles from input instead of perceived stereochemistry\n\n"
      ;
    };
  
    virtual const char* SpecificationURL()
    { return "http://www.ncbi.nlm.nih.gov/data_specs/asn/pcsubstance.asn"; };
    // http://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html also useful
    
    virtual unsigned int Flags()
    { return ZEROATOMSOK; };

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    
  private:
    Json::Value inRoot;
    Json::Value outRoot;
    int currentMolIndex;
    
  };
  
  PubChemJSONFormat thePubChemJSONFormat;
  
  bool PubChemJSONFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
        obErrorLog.ThrowError("PubChemJSONFormat", reader.getFormattedErrorMessages(), obError);
        return false;
      }
      ifs.clear();
      ifs.seekg(0, ios::beg);
      currentMolIndex = 0;
    }
    
    // Get the root level of the molecule
    Json::Value molRoot;
    if (!inRoot["PC_Compounds"].empty() && inRoot["PC_Compounds"].isArray()) {
      // File is a PC_Compounds array
      if (inRoot["PC_Compounds"].size() > currentMolIndex) {
        molRoot = inRoot["PC_Compounds"][currentMolIndex];
        currentMolIndex++;
      } else {
        // Finished the last molecule, reset inRoot in case multiple input files
        inRoot.clear();
        ifs.seekg(0, ios::end);
        return false;
      }
    } else if (!inRoot["PC_Substances"].empty() && inRoot["PC_Substances"].isArray()) {
      // File is a PC_Substances array
      if (inRoot["PC_Substances"].size() > currentMolIndex) {
        // We are assuming the first item in compound array is the deposited entry
        molRoot = inRoot["PC_Substances"][currentMolIndex]["compound"][0];
        currentMolIndex++;
      } else {
        // Finished the last molecule, reset inRoot in case multiple input files
        inRoot.clear();
        ifs.seekg(0, ios::end);
        return false;
      }
    } else if (!inRoot["atoms"].empty() && inRoot["atoms"].isObject()) {
      // File is a single molecule, copy and reset inRoot in case multiple input files
      molRoot = Json::Value(inRoot);
      inRoot.clear();
      ifs.seekg(0, ios::end);
    } else {
      obErrorLog.ThrowError("PubChemJSONFormat", "Cannot locate molecule root", obError);
      return false;
    }
    
    // CID or SID
    if (molRoot["id"]["id"]["cid"].isInt()) {
      // PC_Compound
      OBPairInteger *cid = new OBPairInteger;
      cid->SetAttribute("cid");
      cid->SetValue(molRoot["id"]["id"]["cid"].asInt());
      cid->SetOrigin(fileformatInput);
      pmol->SetData(cid);
      ostringstream s;
      s << molRoot["id"]["id"]["cid"].asInt();
      std::string title(s.str());
      pmol->SetTitle(title);
    } else if (inRoot["PC_Substances"][currentMolIndex-1]["sid"]["id"].isInt()) {
      // PC_Substace
      OBPairInteger *sid = new OBPairInteger;
      sid->SetAttribute("sid");
      sid->SetValue(inRoot["PC_Substances"][currentMolIndex-1]["sid"]["id"].asInt());
      sid->SetOrigin(fileformatInput);
      pmol->SetData(sid);
      ostringstream s;
      s << inRoot["PC_Substances"][currentMolIndex-1]["sid"]["id"].asInt();
      std::string title(s.str());
      pmol->SetTitle(title);
    }
    
    // Atom elements
    Json::Value eAids = molRoot["atoms"]["aid"];
    Json::Value elements = molRoot["atoms"]["element"];
    pmol->ReserveAtoms(eAids.size());
    for(Json::ArrayIndex i = 0; i < eAids.size(); i++) {
      if (eAids[i].isInt() && elements[i].isString()) {
        string elementstring = elements[i].asString();
        OBAtom* patom = pmol->NewAtom((unsigned long)eAids[i].asInt());
        if (elementstring == "a" || elementstring == "d" || elementstring == "r" || elementstring == "lp") {
          patom->SetAtomicNum(0);
        } else {
          int isotope = 0;
          patom->SetAtomicNum(etab.GetAtomicNum(elements[i].asString(), isotope));
          if (isotope != 0) {
            // isotope gets set for e.g. 'D' or 'T' symbol
            patom->SetIsotope(isotope);
          }
        }
        
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom", obWarning);
      }
    }
    
    // Atom charges
    Json::Value charges = molRoot["atoms"]["charge"];
    for(Json::ArrayIndex i = 0; i < charges.size(); i++) {
      Json::Value charge = charges[i];
      if (charge["aid"].isInt() && charge["value"].isInt()) {
        OBAtom* patom = pmol->GetAtomById(charge["aid"].asInt());
        if (patom) {
          patom->SetFormalCharge(charge["value"].asInt());
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom charge", obWarning);
        }
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom charge", obWarning);
      }
    }
    
    // Atom radicals
    Json::Value radicals = molRoot["atoms"]["radical"];
    for(Json::ArrayIndex i = 0; i < radicals.size(); i++) {
      Json::Value radical = radicals[i];
      if (radical["aid"].isInt() && radical["type"].isString()) {
        OBAtom* patom = pmol->GetAtomById(radical["aid"].asInt());
        if (patom) {
          string radicalstring = radical["type"].asString();
          if (radicalstring == "singlet") {
            patom->SetSpinMultiplicity(1);
          } else if (radicalstring == "doublet") {
            patom->SetSpinMultiplicity(2); 
          } else if (radicalstring == "triplet") {
            patom->SetSpinMultiplicity(3);
          } else if (radicalstring == "quartet") {
            patom->SetSpinMultiplicity(4);
          } else if (radicalstring == "quintet") {
            patom->SetSpinMultiplicity(5);
          } else if (radicalstring == "hextet") {
            patom->SetSpinMultiplicity(6);
          } else if (radicalstring == "heptet") {
            patom->SetSpinMultiplicity(7);
          } else if (radicalstring == "octet") {
            patom->SetSpinMultiplicity(8);
          } else if (radicalstring == "none") {
            patom->SetSpinMultiplicity(0);
          } else {
            obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom radical", obWarning);
          }
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom radical", obWarning);
        }
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom radical", obWarning);
      }
    }
    
    // Atom isotopes
    Json::Value isotopes = molRoot["atoms"]["isotope"];
    for(Json::ArrayIndex i = 0; i < isotopes.size(); i++) {
      Json::Value isotope = isotopes[i];
      if (isotope["aid"].isInt() && isotope["value"].isInt()) {
        OBAtom* patom = pmol->GetAtomById(isotope["aid"].asInt());
        if (patom) {
          patom->SetIsotope(isotope["value"].asInt());
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom isotope", obWarning);
        }
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom isotope", obWarning);
      }
    }
    
    // TODO: atom label, atom comment
    // array of (aid, value<string>)
    
    // Bond orders
    Json::Value oAid1s = molRoot["bonds"]["aid1"];
    Json::Value oAid2s = molRoot["bonds"]["aid2"];
    Json::Value orders = molRoot["bonds"]["order"];
    for(Json::ArrayIndex i = 0; i < oAid1s.size(); i++) {
      if (oAid1s[i].isInt() && oAid2s[i].isInt() && orders[i].isString()) {
        int order = 0; // Use zero bond order for other bond types (complex, ionic, dative, unknown)
        string orderstring = orders[i].asString();
        if (orderstring == "single") {
          order = 1;
        } else if (orderstring == "double") {
          order = 2;
        } else if (orderstring == "triple") {
          order = 3;
        } else if (orderstring == "quadruple") {
          order = 4;
        }
        OBAtom* beginAtom = pmol->GetAtomById(oAid1s[i].asInt());
        OBAtom* endAtom = pmol->GetAtomById(oAid2s[i].asInt());
        if (beginAtom && endAtom) {
          OBBond* pbond = pmol->NewBond();
          pbond->SetBegin(beginAtom);
          pbond->SetEnd(endAtom);
          pbond->SetBondOrder(order);
          beginAtom->AddBond(pbond);
          endAtom->AddBond(pbond);
          // Save type string as generic data on bond (useful for non-standard bonds)
          OBPairData *bondType = new OBPairData;
          bondType->SetAttribute("type");
          bondType->SetValue(orderstring);
          bondType->SetOrigin(fileformatInput);
          pbond->SetData(bondType);
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid bond", obWarning);
        }
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid bond", obWarning);
      }
    }
    
    int dim = 0;    // Set dimension to 0 unless coordinates are found
    
    // Atom coordinates
    Json::Value conf = molRoot["coords"][0]["conformers"][0];
    Json::Value cAids = molRoot["coords"][0]["aid"];
    if (!conf["x"].empty() && !conf["y"].empty()) {
      dim = 2;
    }
    if (!conf["z"].empty()) {
      dim = 3;
    }
    for(Json::ArrayIndex i = 0; i < cAids.size(); i++) {
      double x,y,z = 0;
      x = conf["x"].get(i, 0).asDouble();
      y = conf["y"].get(i, 0).asDouble();
      z = conf["z"].get(i, 0).asDouble();
      OBAtom* patom = pmol->GetAtomById(cAids[i].asInt());
      if (patom) {
        patom->SetVector(x,y,z);
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid coordinates", obWarning);
      }
    }
    
    // TODO: Coordinates type
    //Json::Value type = molRoot["coords"][0]["type"];
    // Generic data? An array of strings
    // twod, threed, submitted, experimental, computed, standardized, augmented, aligned, compact,
    // units-unknown, units-angstroms, units-nanometers, units-pixel, units-points, units-stdbonds
    
    // Bond styles and annotations
    Json::Value aid1s = conf["style"]["aid1"];
    Json::Value aid2s = conf["style"]["aid2"];
    Json::Value styles = conf["style"]["annotation"];
    for(Json::ArrayIndex i = 0; i < aid1s.size(); i++) {
      if (aid1s[i].isInt() && aid2s[i].isInt() && styles[i].isString()) {
        OBAtom* beginAtom = pmol->GetAtomById(aid1s[i].asInt());
        OBAtom* endAtom = pmol->GetAtomById(aid2s[i].asInt());
        if (beginAtom && endAtom) {
          OBBond* pbond = pmol->GetBond(beginAtom, endAtom);
          if (!pbond) {
            // Create zero order bond if none exists
            OBBond* pbond = pmol->NewBond();
            pbond->SetBondOrder(0);
            beginAtom->AddBond(pbond);
            endAtom->AddBond(pbond);
          }
          // Use annotations to add stereo information
          unsigned int flags = pbond->GetFlags();
          string stylestring = styles[i].asString();
          if (stylestring == "aromatic") {
            flags |= OBBond::Aromatic;
          } else if (stylestring == "wedge-up") {
            flags |= OBBond::Wedge;
          } else if (stylestring == "wedge-down") {
            flags |= OBBond::Hash;
          } else if (stylestring == "crossed") {
            flags |= OBBond::CisOrTrans;
          } else if (stylestring == "wavy") {
            flags |= OBBond::WedgeOrHash;
          } else if (stylestring == "dashed" || stylestring == "dotted" || 
                     stylestring == "arrow" || stylestring == "resonance" || 
                     stylestring == "bold" || stylestring == "fischer" || 
                     stylestring == "closeContact" || stylestring == "unknown") {
            // Save non-standard annotations as generic data on bond (multiple possible)
            vector<string> val;
            if (pbond->HasData("style")) {
              AnnotationData *data = dynamic_cast<AnnotationData*>(pbond->GetData("style"));
              val = data->GetGenericValue();
              pbond->DeleteData("style");
            }
            AnnotationData *data = new AnnotationData;
            data->SetAttribute("style");
            data->SetOrigin(fileformatInput);
            val.push_back(stylestring);
            data->SetValue(val);
            pbond->SetData(data);
          }
          pbond->Set(pbond->GetIdx(), beginAtom, endAtom, pbond->GetBondOrder(), flags);    
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid bond style", obWarning);
        }
      }
    }
    pmol->SetDimension(dim);
    
    // Total molecular charge (redundant due to atom charges?)
    if (!molRoot["charge"].empty() && molRoot["charge"].isInt()) {
      pmol->SetTotalCharge(molRoot["charge"].asInt());
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
      if (flag & OBBond::CisOrTrans && pbond->GetBondOrder()==2)
        bd = OBStereo::UnknownDir;
      if (bd != OBStereo::NotStereo)
        updown[&*pbond] = bd;
    }
    
    pmol->EndModify();
    
    if (pConv->IsOption("s",OBConversion::INOPTIONS)) {
      // Use the stereo information in the input file
      pmol->DeleteData(OBGenericDataType::StereoData);
      Json::Value stereos = molRoot["stereo"];
      for(Json::ArrayIndex i = 0; i < stereos.size(); i++) {
        Json::Value stereo = stereos[i];
        if (stereo.isMember("tetrahedral")) {
          Json::Value tet = stereo["tetrahedral"];
          OBTetrahedralStereo::Config config;
          config.center = tet["center"].asInt();
          config.from = (tet["top"].asInt() == -1) ? OBStereo::ImplicitRef : tet["top"].asInt();
          config.refs.push_back((tet["below"].asInt() == -1) ? OBStereo::ImplicitRef : tet["below"].asInt());
          if (tet["parity"].asString() == "clockwise") {
            config.specified = true;
            config.winding = OBStereo::Clockwise;
            config.refs.push_back((tet["bottom"].asInt() == -1) ? OBStereo::ImplicitRef : tet["bottom"].asInt());
            config.refs.push_back((tet["above"].asInt() == -1) ? OBStereo::ImplicitRef : tet["above"].asInt());
          } else if (tet["parity"].asString() == "counterclockwise") {
            config.specified = true;
            config.winding = OBStereo::AntiClockwise;
            config.refs.push_back((tet["above"].asInt() == -1) ? OBStereo::ImplicitRef : tet["above"].asInt());
            config.refs.push_back((tet["bottom"].asInt() == -1) ? OBStereo::ImplicitRef : tet["bottom"].asInt());
          } else {
            config.specified = false;
            config.winding = OBStereo::UnknownWinding;
            config.refs.push_back((tet["bottom"].asInt() == -1) ? OBStereo::ImplicitRef : tet["bottom"].asInt());
            config.refs.push_back((tet["above"].asInt() == -1) ? OBStereo::ImplicitRef : tet["above"].asInt());
          }
          OBTetrahedralStereo *ts = new OBTetrahedralStereo(pmol);
          ts->SetConfig(config);
          pmol->SetData(ts);
        } else if (stereo.isMember("planar")) {
          Json::Value pl = stereo["planar"];
          OBCisTransStereo::Config config;
          config.begin = pl["left"].asInt();
          config.end = pl["right"].asInt();
          config.refs.push_back((pl["ltop"].asInt() == -1) ? OBStereo::ImplicitRef : pl["ltop"].asInt());
          config.refs.push_back((pl["rtop"].asInt() == -1) ? OBStereo::ImplicitRef : pl["rtop"].asInt());
          config.refs.push_back((pl["rbottom"].asInt() == -1) ? OBStereo::ImplicitRef : pl["rbottom"].asInt());
          config.refs.push_back((pl["lbottom"].asInt() == -1) ? OBStereo::ImplicitRef : pl["lbottom"].asInt());
          if (pl["parity"].asString() == "any" || pl["parity"].asString() == "unknown") {
            config.specified = false;
          } else {
            config.specified = true;
            config.shape = OBStereo::ShapeU;
          }
          OBCisTransStereo *ct = new OBCisTransStereo(pmol);
          ct->SetConfig(config);
          pmol->SetData(ct);
        } else if (stereo.isMember("squareplanar")) {
          Json::Value sq = stereo["squareplanar"];
          OBSquarePlanarStereo::Config config;
          config.center = sq["center"].asInt();
          config.refs.push_back((sq["lbelow"].asInt() == -1) ? OBStereo::ImplicitRef : sq["lbelow"].asInt());
          config.refs.push_back((sq["rbelow"].asInt() == -1) ? OBStereo::ImplicitRef : sq["rbelow"].asInt());
          config.refs.push_back((sq["rabove"].asInt()) ? OBStereo::ImplicitRef : sq["rabove"].asInt());
          config.refs.push_back((sq["labove"].asInt() == -1) ? OBStereo::ImplicitRef : sq["labove"].asInt());
          if (sq["parity"].asString() == "any" || sq["parity"].asString() == "unknown") {
            config.specified = false;
          } else {
            config.specified = true;
            config.shape = OBStereo::ShapeU;
          }
          OBSquarePlanarStereo *ss = new OBSquarePlanarStereo(pmol);
          ss->SetConfig(config);
          pmol->SetData(ss);
        } else if (stereo.isMember("octahedral")) {
          obErrorLog.ThrowError("PubChemJSONFormat", "Octahedral stereochemistry not implemented", obWarning);
          // aids: center, top, bottom, lbelow, rbelow, labove, rabove
        } else if (stereo.isMember("bipyramid")) {
          obErrorLog.ThrowError("PubChemJSONFormat", "Bipyramidal stereochemistry not implemented", obWarning);
          // aids: above, below, bottom, center, top, right
        } else if (stereo.isMember("tshape")) {
          obErrorLog.ThrowError("PubChemJSONFormat", "T shape stereochemistry not implemented", obWarning);
          // aids: center, top, bottom, above
        } else if (stereo.isMember("pentagonal")) {
          obErrorLog.ThrowError("PubChemJSONFormat", "Pentagonal stereochemistry not implemented", obWarning);
          // aids: center, top, bottom, left, lbelow, rbelow, labove, rabove
        }
      }
      pmol->SetChiralityPerceived();
    } else {
      // Use OB stereo perception to get stereo from coordinates and bond styles
      if (pmol->Has3D()) {
        // Use 3D coordinates to determine stereochemistry
        StereoFrom3D(pmol);
        // For unspecified cis/trans stereos, set their Configs to unspecified
        map<OBBond*, OBStereo::BondDirection>::const_iterator bd_it;
        OpenBabel::OBStereoFacade facade(pmol);
        for(bd_it=updown.begin(); bd_it!=updown.end(); ++bd_it) {
          OBBond* bond = bd_it->first;
          if (bond->GetBondOrder()!=2 || bd_it->second != OBStereo::UnknownDir)
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
    }
    
    // TODO: Properties
    
    return true;
  }
  
  bool PubChemJSONFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol == NULL)
      return false;
    ostream& ofs = *pConv->GetOutStream();
    
    if (pmol->GetDimension()==0) {
      obErrorLog.ThrowError("PubChemJSONFormat", "No 2D or 3D coordinates exist. "
                 "To generate 2D or 3D coordinates use --gen2D or --gen3D.", obError);
      return false;
    }
    
    // Make all hydrogens explicit by default, unless -d option set
    if (pConv->IsOption("d", OBConversion::GENOPTIONS)) {
      obErrorLog.ThrowError("PubChemJSONFormat", "Stereo output may be invalid due to implicit hydrogens", obWarning);
    } else {
      pmol->AddHydrogens();
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
    
    Json::Value doc;  // Root of molecule JSON
    
    // CID
    if (pmol->HasData("cid")) {
      OBPairInteger *cid = dynamic_cast<OBPairInteger*>(pmol->GetData("cid"));
      doc["id"]["id"]["cid"] = cid->GetGenericValue();
    }
    
    // Atoms
    int id = 1;
    FOR_ATOMS_OF_MOL(patom, pmol) {
      // Id (Overwrite Id from input to ensure consecutive integers 1+)
      patom->SetId(id);
      doc["atoms"]["aid"].append(id);
      // Element
      if (patom->GetAtomicNum()) {
        string el = etab.GetSymbol(patom->GetAtomicNum());
        std::transform(el.begin(), el.end(), el.begin(), ::tolower);
        doc["atoms"]["element"].append(el);
      } else {
        doc["atoms"]["element"].append("a");
      }
      // Charge
      int c = patom->GetFormalCharge();
      if (c != 0) {
        Json::Value charge;
        charge["aid"] = id;
        charge["value"] = c;
        doc["atoms"]["charge"].append(charge);
      }
      // Isotope
      int m = patom->GetIsotope();
      if (m != 0) {
        Json::Value isotope;
        isotope["aid"] = id;
        isotope["value"] = m;
        doc["atoms"]["isotope"].append(isotope);
      }
      // Radical
      int sm = patom->GetSpinMultiplicity();
      if (sm > 0 && sm < 9) {
        Json::Value radical;
        radical["aid"] = id;
        if (sm == 1) {
          radical["type"] = "singlet";
        } else if (sm == 2) {
          radical["type"] = "doublet";
        } else if (sm == 3) {
          radical["type"] = "triplet";
        } else if (sm == 4) {
          radical["type"] = "quartet";
        } else if (sm == 5) {
          radical["type"] = "quintet";
        } else if (sm == 6) {
          radical["type"] = "hextet";
        } else if (sm == 7) {
          radical["type"] = "heptet";
        } else if (sm == 8) {
          radical["type"] = "octet";
        }
        doc["atoms"]["radical"].append(radical);
      }
      // Coordinates
      // TODO: An option to round coordinates to n decimal places?
      doc["coords"][0]["aid"].append(id);
      doc["coords"][0]["conformers"][0]["x"].append(patom->GetX());
      doc["coords"][0]["conformers"][0]["y"].append(patom->GetY());
      if (pmol->GetDimension() == 3) {
        doc["coords"][0]["conformers"][0]["z"].append(patom->GetY());
      }
      
      id++;
    }

    // Bonds
    FOR_BONDS_OF_MOL(pbond, pmol) {
    
      // Order
      int aid1 = (int)pbond->GetBeginAtom()->GetId();
      int aid2 = (int)pbond->GetEndAtom()->GetId();
      doc["bonds"]["aid1"].append(aid1);
      doc["bonds"]["aid2"].append(aid2);
      int order = pbond->GetBondOrder();
      if (order == 0) {
        if (pbond->HasData("type")) {
          // Check to see if a "type" string exists
          OBPairData *id = dynamic_cast<OBPairData*>(pbond->GetData("type"));
          doc["bonds"]["order"].append(id->GetValue());
        } else {
          doc["bonds"]["order"].append("unknown");
        }
      } else if (order == 1) {
        doc["bonds"]["order"].append("single");
      } else if (order == 2) {
        doc["bonds"]["order"].append("double");
      } else if (order == 3) {
        doc["bonds"]["order"].append("triple");
      } else if (order == 4) {
        doc["bonds"]["order"].append("quadruple");
      }
      
      // Styles and annotations
      vector<string> annotations;
      if (pConv->IsOption("w", pConv->OUTOPTIONS)) {
        // option w means just use input bond stereo annotations
        if (pbond->IsWedge()) {
          annotations.push_back("wedge-up");
        } else if (pbond->IsHash()) {
          annotations.push_back("wedge-down");
        } else if (pbond->IsWedgeOrHash()) {
          annotations.push_back("wavy");
        } else if (pbond->IsCisOrTrans()) {
          annotations.push_back("crossed");
        }
      } else {
        // No option w means use stereochemistry information
        from_cit = from.find(&*pbond);
        if (from_cit != from.end() && from_cit->second == pbond->GetEndAtom()->GetId()) {
          swap(aid1, aid2);  // Swap start and end atom if necessary
        }
        if (unspec_ctstereo.find(&*pbond) != unspec_ctstereo.end()) {
          annotations.push_back("crossed");
        }  
        if (updown.find(&*pbond) != updown.end()) {
          if (updown[&*pbond] == 1) {
            annotations.push_back("wedge-up");
          } else if (updown[&*pbond] == 4) {
            annotations.push_back("wavy");
          } else if (updown[&*pbond] == 6) {
            annotations.push_back("wedge-down");
          }
        }
      }
      if (pbond->IsAromatic()) {
        annotations.push_back("aromatic");
      }
      if (pbond->HasData("style")) {
        AnnotationData *data = dynamic_cast<AnnotationData*>(pbond->GetData("style"));
        vector<string> styles = data->GetGenericValue();
        for(vector<string>::const_iterator i = styles.begin(); i != styles.end(); ++i) {
          annotations.push_back(*i);
        }
      }
      annotations.erase(unique(annotations.begin(), annotations.end()), annotations.end());
      for(vector<string>::const_iterator i = annotations.begin(); i != annotations.end(); ++i) {
        doc["coords"][0]["conformers"][0]["style"]["aid1"].append(aid1);
        doc["coords"][0]["conformers"][0]["style"]["aid2"].append(aid2);
        doc["coords"][0]["conformers"][0]["style"]["annotation"].append(*i);
      }
    }
    
    // Stereochemistry
    OBStereoFacade facade(pmol);
    FOR_ATOMS_OF_MOL(patom, pmol) {
      if (facade.HasTetrahedralStereo(patom->GetId())) {
        OBTetrahedralStereo::Config config = facade.GetTetrahedralStereo(patom->GetId())->GetConfig();
        Json::Value tet;
        tet["tetrahedral"]["type"] = "tetrahedral";
        tet["tetrahedral"]["center"] = (int)config.center;
        tet["tetrahedral"]["top"] = (int)config.from;
        tet["tetrahedral"]["below"] = (int)config.refs[0];
        tet["tetrahedral"]["bottom"] = (int)config.refs[1];
        tet["tetrahedral"]["above"] = (int)config.refs[2];
        if (config.winding == OBStereo::UnknownWinding || !config.specified) {
          tet["tetrahedral"]["parity"] = "any";
        } else if (config.winding == OBStereo::Clockwise){
          tet["tetrahedral"]["parity"] = "clockwise";
        } else if (config.winding == OBStereo::AntiClockwise){
          tet["tetrahedral"]["parity"] = "counterclockwise";
          tet["tetrahedral"]["bottom"] = (int)config.refs[2];
          tet["tetrahedral"]["above"] = (int)config.refs[1];
        }
        doc["stereo"].append(tet);
      }
      if (facade.HasSquarePlanarStereo(patom->GetId())) {
        OBSquarePlanarStereo *sqs = facade.GetSquarePlanarStereo(patom->GetId());
        OBSquarePlanarStereo::Config config = sqs->GetConfig();
        Json::Value sq;
        sq["squareplanar"]["center"] = (int)config.center;
        if (config.specified) {
          if (config.shape == OBStereo::ShapeU) {
            sq["squareplanar"]["parity"] = "u-shape";
            sq["squareplanar"]["lbelow"] = (int)config.refs[0];
            sq["squareplanar"]["rbelow"] = (int)config.refs[1];
            sq["squareplanar"]["rabove"] = (int)config.refs[2];
            sq["squareplanar"]["labove"] = (int)config.refs[3];
          } else if (config.shape == OBStereo::ShapeZ) {
            sq["squareplanar"]["parity"] = "z-shape";
            sq["squareplanar"]["lbelow"] = (int)config.refs[0];
            sq["squareplanar"]["rbelow"] = (int)config.refs[1];
            sq["squareplanar"]["labove"] = (int)config.refs[2];
            sq["squareplanar"]["rabove"] = (int)config.refs[3];
          } else if (config.shape == OBStereo::Shape4) {
            sq["squareplanar"]["parity"] = "x-shape";
            sq["squareplanar"]["lbelow"] = (int)config.refs[0];
            sq["squareplanar"]["rabove"] = (int)config.refs[1];
            sq["squareplanar"]["rbelow"] = (int)config.refs[2];
            sq["squareplanar"]["labove"] = (int)config.refs[3];
          }
        } else {
          sq["squareplanar"]["parity"] = "any";
          sq["squareplanar"]["lbelow"] = (int)config.refs[0];
          sq["squareplanar"]["rbelow"] = (int)config.refs[1];
          sq["squareplanar"]["rabove"] = (int)config.refs[2];
          sq["squareplanar"]["labove"] = (int)config.refs[3];
        }
      }
    }
    FOR_BONDS_OF_MOL(pbond, pmol) {
      if (facade.HasCisTransStereo(pbond->GetId())) {
        OBCisTransStereo *cts = facade.GetCisTransStereo(pbond->GetId());
        OBCisTransStereo::Config config = cts->GetConfig();
        Json::Value ct;
        ct["planar"]["type"] = "planar";
        ct["planar"]["ltop"] = (int)config.refs[0];
        OBAtom *begin = pmol->GetAtomById(config.begin);
        OBAtom *ltop = pmol->GetAtomById(config.refs[0]);
        if (begin && ltop) {
           if (ltop->IsConnected(begin)) {
            ct["planar"]["left"] = (int)config.begin;
            ct["planar"]["right"] = (int)config.end;
          } else {
            ct["planar"]["left"] = (int)config.end;
            ct["planar"]["right"] = (int)config.begin;
          }
          ct["planar"]["rbottom"] = (int)cts->GetTransRef(ct["planar"]["ltop"].asInt());
          ct["planar"]["lbottom"] = (int)cts->GetCisRef(ct["planar"]["rbottom"].asInt());
          ct["planar"]["rtop"] = (int)cts->GetTransRef(ct["planar"]["lbottom"].asInt());
          if (config.specified) {
            // Open babel is not capable of determining parity? (need CIP rules?)
            ct["planar"]["parity"] = "unknown";
          } else {
            ct["planar"]["parity"] = "any";
          }        
          doc["stereo"].append(ct);
        }
      }
    }
    
    doc["charge"] = pmol->GetTotalCharge();
    
    outRoot["PC_Compounds"].append(doc);
    
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
