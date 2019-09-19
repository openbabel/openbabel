/**********************************************************************
Copyright (C) 2013, 2017, 2018 by Matt Swain <m.swain@me.com>

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
#include <openbabel/stereo/squareplanar.h>

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
      "The JSON format returned by the PubChem PUG REST service\n\n"

      "The data contained in this format closely resembles PubChem's internal data structure.\n\n"

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

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    rapidjson::Document inRoot;
    rapidjson::Document outRoot;
    int currentMolIndex;

  };
  
  PubChemJSONFormat thePubChemJSONFormat;
  
  bool PubChemJSONFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
        obErrorLog.ThrowError("PubChemJSONFormat", msg.str(), obError);
        return false;
      }
      // Clear ifs flags so it is "good" and any subsequent mols are read, but leave it at EOF position.
      // Therefore, when not at EOF position we know to parse next file and reset currentMolIndex
      ifs.clear();
      currentMolIndex = 0;
    }

    if (!inRoot.IsObject()) {
      obErrorLog.ThrowError("PubChemJSONFormat", "JSON file should be a single object", obError);
      return false;
    }

    // Get the root level of the molecule
    rapidjson::Value molRoot;
    rapidjson::Value subsRoot;
    if (inRoot.HasMember("PC_Compounds") && inRoot["PC_Compounds"].IsArray()) {
      // File is a PC_Compounds array
      if (inRoot["PC_Compounds"].Size() > currentMolIndex) {
        molRoot = inRoot["PC_Compounds"][currentMolIndex];
        currentMolIndex++;
      } else {
        // Finished the last molecule
        ifs.seekg(0, ios::end);
        return false;
      }
    } else if (inRoot.HasMember("PC_Substances") && inRoot["PC_Substances"].IsArray()) {
      // File is a PC_Substances array
      if (inRoot["PC_Substances"].Size() > currentMolIndex) {
        subsRoot = inRoot["PC_Substances"][currentMolIndex];
        // We are assuming the first item in compound array is the deposited entry
        molRoot = subsRoot["compound"][0];
        currentMolIndex++;
      } else {
        // Finished the last molecule
        ifs.seekg(0, ios::end);
        return false;
      }
    } else {
      obErrorLog.ThrowError("PubChemJSONFormat", "JSON file must contain a PC_Compounds or PC_Substances array",
                            obError);
    }

    // CID or SID
    if (molRoot.HasMember("id") && molRoot["id"].HasMember("id") &&
        molRoot["id"]["id"].HasMember("cid") && molRoot["id"]["id"]["cid"].IsInt()) {
      // PC_Compound
      ostringstream s;
      s << molRoot["id"]["id"]["cid"].GetInt();
      std::string title(s.str());
      pmol->SetTitle(title);
      OBPairData *cid = new OBPairData;
      cid->SetAttribute("cid");
      cid->SetValue(title);
      cid->SetOrigin(fileformatInput);
      pmol->SetData(cid);
    } else if (!subsRoot.IsNull() && subsRoot.HasMember("sid") &&
               subsRoot["sid"].HasMember("id") && subsRoot["sid"]["id"].IsInt()) {
      // PC_Substance
      ostringstream s;
      s << subsRoot["sid"]["id"].GetInt();
      std::string title(s.str());
      pmol->SetTitle(title);
      OBPairData *sid = new OBPairData;
      sid->SetAttribute("sid");
      sid->SetValue(title);
      sid->SetOrigin(fileformatInput);
      pmol->SetData(sid);
    } else {
      obErrorLog.ThrowError("PubChemJSONFormat", "Did not find CID or SID in PubChem JSON", obWarning);
    }

    // Atom elements
    if (!molRoot.HasMember("atoms") || !molRoot["atoms"].HasMember("aid") || !molRoot["atoms"].HasMember("element")) {
      obErrorLog.ThrowError("PubChemJSONFormat", "Could not read atom data", obError);
    }
    const rapidjson::Value &eAids = molRoot["atoms"]["aid"];
    const rapidjson::Value &elements = molRoot["atoms"]["element"];
    pmol->ReserveAtoms(eAids.Size());
    for (rapidjson::SizeType i = 0; i < eAids.Size(); i++) {
      if (eAids[i].IsInt() && elements[i].IsInt()) {
        // Element provided as integer atomic number
        int atomicNum = elements[i].GetInt();
        OBAtom *patom = pmol->NewAtom(eAids[i].GetUint());
        if (atomicNum == 255 || atomicNum == 254 || atomicNum == 253 || atomicNum == 252) {
          patom->SetAtomicNum(0);
        } else {
          patom->SetAtomicNum(atomicNum);
        }
      } else if (eAids[i].IsInt() && elements[i].IsString()) {
        // Element provided as string (old format)
        string elementstring = elements[i].GetString();
        OBAtom *patom = pmol->NewAtom(eAids[i].GetUint());
        if (elementstring == "a" || elementstring == "d" || elementstring == "r" || elementstring == "lp") {
          patom->SetAtomicNum(0);
        } else {
          // Ensure first letter is uppercase
          elementstring[0] = toupper(elementstring[0]);
          patom->SetAtomicNum(OBElements::GetAtomicNum(elementstring.c_str()));
        }
      } else {
        obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom", obWarning);
      }
    }

    // Atom charges
    if (molRoot["atoms"].HasMember("charge") && molRoot["atoms"]["charge"].IsArray()) {
      const rapidjson::Value &charges = molRoot["atoms"]["charge"];
      for (rapidjson::SizeType i = 0; i < charges.Size(); i++) {
        const rapidjson::Value &charge = charges[i];
        if (charge["aid"].IsInt() && charge["value"].IsInt()) {
          OBAtom *patom = pmol->GetAtomById(charge["aid"].GetUint());
          if (patom) {
            patom->SetFormalCharge(charge["value"].GetInt());
          } else {
            obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom charge", obWarning);
          }
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom charge", obWarning);
        }
      }
    }

    // Atom radicals
    if (molRoot["atoms"].HasMember("radical") && molRoot["atoms"]["radical"].IsArray()) {
      rapidjson::Value &radicals = molRoot["atoms"]["radical"];
      for (rapidjson::SizeType i = 0; i < radicals.Size(); i++) {
        rapidjson::Value &radical = radicals[i];
        if (radical["aid"].IsInt() && radical["type"].IsInt()) {
          // Radical provided as integer
          OBAtom *patom = pmol->GetAtomById(radical["aid"].GetUint());
          if (patom) {
            int sm = radical["type"].GetInt();
            if (sm == 255) {
              sm = 0;
            }
            patom->SetSpinMultiplicity(sm);
          } else {
            obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom radical", obWarning);
          }
        } else if (radical["aid"].IsInt() && radical["type"].IsString()) {
          // Radical provided as string (old format)
          OBAtom *patom = pmol->GetAtomById(radical["aid"].GetUint());
          if (patom) {
            string radicalstring = radical["type"].GetString();
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
    }

    // Atom isotopes
    if (molRoot["atoms"].HasMember("isotope") && molRoot["atoms"]["isotope"].IsArray()) {
      rapidjson::Value &isotopes = molRoot["atoms"]["isotope"];
      for (rapidjson::SizeType i = 0; i < isotopes.Size(); i++) {
        rapidjson::Value &isotope = isotopes[i];
        if (isotope["aid"].IsInt() && isotope["value"].IsInt()) {
          OBAtom *patom = pmol->GetAtomById(isotope["aid"].GetUint());
          if (patom) {
            patom->SetIsotope(isotope["value"].GetUint());
          } else {
            obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom isotope", obWarning);
          }
        } else {
          obErrorLog.ThrowError("PubChemJSONFormat", "Invalid atom isotope", obWarning);
        }
      }
    }

    // TODO: atom label, atom comment
    // array of (aid, value<string>)

    // Bond orders
    if (molRoot.HasMember("bonds") && molRoot["bonds"].HasMember("aid1") && molRoot["bonds"].HasMember("aid2") &&
        molRoot["bonds"].HasMember("order")) {
      rapidjson::Value &oAid1s = molRoot["bonds"]["aid1"];
      rapidjson::Value &oAid2s = molRoot["bonds"]["aid2"];
      rapidjson::Value &orders = molRoot["bonds"]["order"];
      for (rapidjson::SizeType i = 0; i < oAid1s.Size(); i++) {
        if (oAid1s[i].IsInt() && oAid2s[i].IsInt() && orders[i].IsInt()) {
          // Bond order provided as integer
          OBAtom *beginAtom = pmol->GetAtomById(oAid1s[i].GetUint());
          OBAtom *endAtom = pmol->GetAtomById(oAid2s[i].GetUint());
          if (beginAtom && endAtom) {
            OBBond *pbond = pmol->NewBond();
            pbond->SetBegin(beginAtom);
            pbond->SetEnd(endAtom);
            int order = orders[i].GetInt();
            // Other bond types: dative (5), complex (6), ionic (7), unknown (255)
            if (order > 4) {
              // Save type string as generic data on bond for non-standard bonds
              string orderstring = "unknown";
              if (order == 5) {
                orderstring = "dative";
              } else if (order == 6) {
                orderstring = "complex";
              } else if (order == 7) {
                orderstring = "ionic";
              }
              OBPairData *bondType = new OBPairData;
              bondType->SetAttribute("type");
              bondType->SetValue(orderstring);
              bondType->SetOrigin(fileformatInput);
              pbond->SetData(bondType);
              // Use zero bond order for non-standard bonds
              order = 0;
            }
            pbond->SetBondOrder(order);
            beginAtom->AddBond(pbond);
            endAtom->AddBond(pbond);
          } else {
            obErrorLog.ThrowError("PubChemJSONFormat", "Invalid bond", obWarning);
          }
        } else if (oAid1s[i].IsInt() && oAid2s[i].IsInt() && orders[i].IsString()) {
          // Bond order provided as string (old format)
          int order = 0; // Use zero bond order for other bond types (complex, ionic, dative, unknown)
          string orderstring = orders[i].GetString();
          if (orderstring == "single") {
            order = 1;
          } else if (orderstring == "double") {
            order = 2;
          } else if (orderstring == "triple") {
            order = 3;
          } else if (orderstring == "quadruple") {
            order = 4;
          }
          OBAtom *beginAtom = pmol->GetAtomById(oAid1s[i].GetUint());
          OBAtom *endAtom = pmol->GetAtomById(oAid2s[i].GetUint());
          if (beginAtom && endAtom) {
            OBBond *pbond = pmol->NewBond();
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
    }

    unsigned short dim = 0;    // Set dimension to 0 unless coordinates are found

    // Atom coordinates
    if (molRoot.HasMember("coords") && molRoot["coords"].IsArray() && molRoot["coords"].Size() > 0 &&
        molRoot["coords"][0].HasMember("aid")) {
      rapidjson::Value &coords = molRoot["coords"][0];
      rapidjson::Value &cAids = coords["aid"];
      if (coords.HasMember("conformers") && coords["conformers"].IsArray() && coords["conformers"].Size() > 0) {
        rapidjson::Value &conf = coords["conformers"][0];
        if (conf.HasMember("x") && conf["x"].IsArray() && conf["x"].Size() >= cAids.Size() && conf.HasMember("y") &&
            conf["y"].IsArray() && conf["y"].Size() >= cAids.Size()) {
          dim = 2;
          if (conf.HasMember("z") && conf["z"].IsArray() && conf["z"].Size() >= cAids.Size()) {
            dim = 3;
          }
          for (rapidjson::SizeType i = 0; i < cAids.Size(); i++) {
            double x, y, z = 0;
            x = conf["x"][i].GetDouble();
            y = conf["y"][i].GetDouble();
            if (dim == 3) {
              z = conf["z"][i].GetDouble();
            }
            OBAtom *patom = pmol->GetAtomById(cAids[i].GetUint());
            if (patom) {
              patom->SetVector(x, y, z);
            } else {
              obErrorLog.ThrowError("PubChemJSONFormat", "Invalid coordinates", obWarning);
            }
          }
        }

        // TODO: Coordinates type
        // rapidjson::Value &type = coords["type"];
        // Generic data? An array of strings
        // twod, threed, submitted, experimental, computed, standardized, augmented, aligned, compact,
        // units-unknown, units-angstroms, units-nanometers, units-pixel, units-points, units-stdbonds

        // Bond styles and annotations
        if (conf.HasMember("style") && conf["style"].IsObject() && conf["style"].HasMember("aid1") &&
            conf["style"]["aid1"].IsArray() &&
            conf["style"].HasMember("aid2") && conf["style"]["aid2"].IsArray() &&
            conf["style"].HasMember("annotation") && conf["style"]["annotation"].IsArray()) {
          rapidjson::Value &aid1s = conf["style"]["aid1"];
          rapidjson::Value &aid2s = conf["style"]["aid2"];
          rapidjson::Value &styles = conf["style"]["annotation"];
          for (rapidjson::SizeType i = 0; i < aid1s.Size(); i++) {
            if (aid1s[i].IsInt() && aid2s[i].IsInt()) {
              OBAtom *beginAtom = pmol->GetAtomById(aid1s[i].GetUint());
              OBAtom *endAtom = pmol->GetAtomById(aid2s[i].GetUint());
              if (beginAtom && endAtom) {
                OBBond *pbond = pmol->GetBond(beginAtom, endAtom);
                if (!pbond) {
                  // Create zero order bond if none exists
                  pbond = pmol->NewBond();
                  pbond->SetBondOrder(0);
                  beginAtom->AddBond(pbond);
                  endAtom->AddBond(pbond);
                }
                // Use annotations to add stereo information
                unsigned int flags = pbond->GetFlags();

                if (styles[i].IsInt()) {
                  // Bond style provided as integer
                  int style = styles[i].GetInt();
                  if (style == 8) {
                    flags |= OBBond::Aromatic;
                  } else if (style == 5) {
                    flags |= OBBond::Wedge;
                  } else if (style == 6) {
                    flags |= OBBond::Hash;
                  } else if (style == 1) {
                    flags |= OBBond::CisOrTrans;
                  } else if (style == 3) {
                    flags |= OBBond::WedgeOrHash;
                  } else {
                    // Save non-standard annotations as generic data on bond (multiple possible)
                    vector<string> val;
                    if (pbond->HasData("style")) {
                      AnnotationData *data = dynamic_cast<AnnotationData *>(pbond->GetData("style"));
                      val = data->GetGenericValue();
                      pbond->DeleteData("style");
                    }
                    AnnotationData *data = new AnnotationData;
                    data->SetAttribute("style");
                    data->SetOrigin(fileformatInput);
                    string stylestring = "unknown";
                    if (style == 2) {
                      stylestring = "dashed";
                    } else if (style == 4) {
                      stylestring = "dotted";
                    } else if (style == 7) {
                      stylestring = "arrow";
                    } else if (style == 9) {
                      stylestring = "resonance";
                    } else if (style == 10) {
                      stylestring = "bold";
                    } else if (style == 11) {
                      stylestring = "fischer";
                    } else if (style == 12) {
                      stylestring = "closeContact";
                    }
                    val.push_back(stylestring);
                    data->SetValue(val);
                    pbond->SetData(data);
                  }
                } else if (styles[i].IsString()) {
                  // Bond style provided as string (old format)
                  string stylestring = styles[i].GetString();
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
                      AnnotationData *data = dynamic_cast<AnnotationData *>(pbond->GetData("style"));
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
                }
                pbond->Set(pbond->GetIdx(), beginAtom, endAtom, pbond->GetBondOrder(), flags);
              } else {
                obErrorLog.ThrowError("PubChemJSONFormat", "Invalid bond style", obWarning);
              }
            }
          }
        }
        pmol->SetDimension(dim);
      }
    }

    // Total molecular charge (redundant due to atom charges?)
    if (molRoot.HasMember("charge") && molRoot["charge"].IsInt()) {
      pmol->SetTotalCharge(molRoot["charge"].GetInt());
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

    pmol->EndModify();

    if (pConv->IsOption("s", OBConversion::INOPTIONS)) {
      // Use the stereo information in the input file
      pmol->DeleteData(OBGenericDataType::StereoData);
      if (molRoot.HasMember("stereo") && molRoot["stereo"].IsArray()) {
        for (rapidjson::SizeType i = 0; i < molRoot["stereo"].Size(); i++) {
          rapidjson::Value &stereo = molRoot["stereo"][i];
          if (stereo.HasMember("tetrahedral")) {
            rapidjson::Value &tet = stereo["tetrahedral"];
            OBTetrahedralStereo::Config config;
            config.center = tet["center"].GetUint();
            config.from = (tet["top"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["top"].GetInt();
            config.refs.push_back((tet["below"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["below"].GetInt());
            if ((tet["parity"].IsInt() && tet["parity"].GetInt() == 1) ||
                (tet["parity"].IsString() && strcmp(tet["parity"].GetString(), "clockwise") == 0)) {
              config.specified = true;
              config.winding = OBStereo::Clockwise;
              config.refs.push_back((tet["bottom"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["bottom"].GetInt());
              config.refs.push_back((tet["above"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["above"].GetInt());
            } else if ((tet["parity"].IsInt() && tet["parity"].GetInt() == 2) ||
                       (tet["parity"].IsString() && strcmp(tet["parity"].GetString(), "counterclockwise") == 0)) {
              config.specified = true;
              config.winding = OBStereo::AntiClockwise;
              config.refs.push_back((tet["above"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["above"].GetInt());
              config.refs.push_back((tet["bottom"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["bottom"].GetInt());
            } else {
              config.specified = false;
              config.winding = OBStereo::UnknownWinding;
              config.refs.push_back((tet["bottom"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["bottom"].GetInt());
              config.refs.push_back((tet["above"].GetInt() == -1) ? OBStereo::ImplicitRef : tet["above"].GetInt());
            }
            OBTetrahedralStereo *ts = new OBTetrahedralStereo(pmol);
            ts->SetConfig(config);
            pmol->SetData(ts);
          } else if (stereo.HasMember("planar")) {
            rapidjson::Value &pl = stereo["planar"];
            OBCisTransStereo::Config config;
            config.begin = pl["left"].GetUint();
            config.end = pl["right"].GetUint();
            config.refs.push_back((pl["ltop"].GetInt() == -1) ? OBStereo::ImplicitRef : pl["ltop"].GetInt());
            config.refs.push_back((pl["rtop"].GetInt() == -1) ? OBStereo::ImplicitRef : pl["rtop"].GetInt());
            config.refs.push_back((pl["rbottom"].GetInt() == -1) ? OBStereo::ImplicitRef : pl["rbottom"].GetInt());
            config.refs.push_back((pl["lbottom"].GetInt() == -1) ? OBStereo::ImplicitRef : pl["lbottom"].GetInt());
            if ((pl["parity"].IsInt() && (pl["parity"].GetInt() == 3 || pl["parity"].GetInt() == 255)) ||
                (pl["parity"].IsString() &&
                 (strcmp(pl["parity"].GetString(), "any") == 0 || strcmp(pl["parity"].GetString(), "unknown") == 0))) {
              config.specified = false;
            } else {
              config.specified = true;
              config.shape = OBStereo::ShapeU;
            }
            OBCisTransStereo *ct = new OBCisTransStereo(pmol);
            ct->SetConfig(config);
            pmol->SetData(ct);
          } else if (stereo.HasMember("squareplanar")) {
            rapidjson::Value &sq = stereo["squareplanar"];
            OBSquarePlanarStereo::Config config;
            config.center = sq["center"].GetUint();
            config.refs.push_back((sq["lbelow"].GetInt() == -1) ? OBStereo::ImplicitRef : sq["lbelow"].GetInt());
            config.refs.push_back((sq["rbelow"].GetInt() == -1) ? OBStereo::ImplicitRef : sq["rbelow"].GetInt());
            config.refs.push_back((sq["rabove"].GetInt()) ? OBStereo::ImplicitRef : sq["rabove"].GetInt());
            config.refs.push_back((sq["labove"].GetInt() == -1) ? OBStereo::ImplicitRef : sq["labove"].GetInt());
            if ((sq["parity"].IsInt() && (sq["parity"].GetInt() == 4 || sq["parity"].GetInt() == 255)) ||
                (sq["parity"].IsString() &&
                 (strcmp(sq["parity"].GetString(), "any") == 0 || strcmp(sq["parity"].GetString(), "unknown") == 0))) {
              config.specified = false;
            } else {
              config.specified = true;
              config.shape = OBStereo::ShapeU;
            }
            OBSquarePlanarStereo *ss = new OBSquarePlanarStereo(pmol);
            ss->SetConfig(config);
            pmol->SetData(ss);
          } else if (stereo.HasMember("octahedral")) {
            obErrorLog.ThrowError("PubChemJSONFormat", "Octahedral stereochemistry not implemented", obWarning);
            // aids: center, top, bottom, lbelow, rbelow, labove, rabove
          } else if (stereo.HasMember("bipyramid")) {
            obErrorLog.ThrowError("PubChemJSONFormat", "Bipyramidal stereochemistry not implemented", obWarning);
            // aids: above, below, bottom, center, top, right
          } else if (stereo.HasMember("tshape")) {
            obErrorLog.ThrowError("PubChemJSONFormat", "T shape stereochemistry not implemented", obWarning);
            // aids: center, top, bottom, above
          } else if (stereo.HasMember("pentagonal")) {
            obErrorLog.ThrowError("PubChemJSONFormat", "Pentagonal stereochemistry not implemented", obWarning);
            // aids: center, top, bottom, left, lbelow, rbelow, labove, rabove
          }
        }
        pmol->SetChiralityPerceived();
      }
    } else {
      // Use OB stereo perception to get stereo from coordinates and bond styles
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
    }

    // TODO: Properties

    // Increment currentMolIndex for next run
    currentMolIndex++;

    return true;
  }

  bool PubChemJSONFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol *pmol = dynamic_cast<OBMol *>(pOb);
    if (pmol == NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    if (pmol->GetDimension() == 0) {
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

    // Set up all the stereochemistry information
    set<OBBond *> unspec_ctstereo = GetUnspecifiedCisTrans(*pmol);
    map<OBBond *, OBStereo::BondDirection> updown;
    map<OBBond *, OBStereo::Ref> from;
    map<OBBond *, OBStereo::Ref>::const_iterator from_cit;
    if (!pConv->IsOption("w", pConv->OUTOPTIONS))
      TetStereoToWedgeHash(*pmol, updown, from);

    // Must always pass an allocator when memory may need to be allocated
    rapidjson::Document::AllocatorType &al = outRoot.GetAllocator();

    rapidjson::Value doc(rapidjson::kObjectType);  // Root of molecule JSON

    // CID
    if (pmol->HasData("cid")) {
      OBPairData *cid = dynamic_cast<OBPairData *>(pmol->GetData("cid"));
      rapidjson::Value cidValue(rapidjson::kStringType);
      cidValue.SetString(cid->GetValue().c_str(), al);
      rapidjson::Value id2(rapidjson::kObjectType);
      id2.AddMember("cid", cidValue, al);
      rapidjson::Value id1(rapidjson::kObjectType);
      id1.AddMember("id", id2, al);
      doc.AddMember("id", id1, al);
    }

    // Atoms
    rapidjson::Value aids(rapidjson::kArrayType);
    rapidjson::Value element(rapidjson::kArrayType);
    rapidjson::Value charge(rapidjson::kArrayType);
    rapidjson::Value isotope(rapidjson::kArrayType);
    rapidjson::Value radical(rapidjson::kArrayType);
    rapidjson::Value xcoords(rapidjson::kArrayType);
    rapidjson::Value ycoords(rapidjson::kArrayType);
    rapidjson::Value zcoords(rapidjson::kArrayType);
    unsigned int id = 1;
    FOR_ATOMS_OF_MOL(patom, pmol) {
      // Id (Overwrite Id from input to ensure consecutive integers 1+)
      patom->SetId(id);
      aids.PushBack(rapidjson::Value(id).Move(), al);
      // Element
      if (patom->GetAtomicNum()) {
        element.PushBack(rapidjson::Value(patom->GetAtomicNum()).Move(), al);
      } else {
        element.PushBack(rapidjson::Value(255).Move(), al);
      }
      // Charge
      int c = patom->GetFormalCharge();
      if (c != 0) {
        rapidjson::Value chg(rapidjson::kObjectType);
        chg.AddMember("aid", rapidjson::Value(id).Move(), al);
        chg.AddMember("value", rapidjson::Value(c).Move(), al);
        charge.PushBack(chg, al);
      }
      // Isotope
      int m = patom->GetIsotope();
      if (m != 0) {
        rapidjson::Value iso(rapidjson::kObjectType);
        iso.AddMember("aid", rapidjson::Value(id).Move(), al);
        iso.AddMember("value", rapidjson::Value(m).Move(), al);
        isotope.PushBack(iso, al);
      }
      // Radical
      int sm = patom->GetSpinMultiplicity();
      if (sm > 0 && sm < 9) {
        rapidjson::Value rad(rapidjson::kObjectType);
        rad.AddMember("aid", rapidjson::Value(id).Move(), al);
        rad.AddMember("type", rapidjson::Value(sm).Move(), al);
        radical.PushBack(rad, al);
      }
      // Coordinates
      // TODO: An option to round coordinates to n decimal places?
      xcoords.PushBack(rapidjson::Value(patom->GetX()).Move(), al);
      ycoords.PushBack(rapidjson::Value(patom->GetX()).Move(), al);
      if (pmol->GetDimension() == 3) {
        zcoords.PushBack(rapidjson::Value(patom->GetZ()).Move(), al);
      }
      id++;
    }

    // Add atoms to doc
    if (aids.Size() > 0) {
      rapidjson::Value atoms(rapidjson::kObjectType);
      atoms.AddMember("aids", aids, al);
      atoms.AddMember("element", element, al);
      if (charge.Size() > 0) {
        atoms.AddMember("charge", charge, al);
      }
      if (isotope.Size() > 0) {
        atoms.AddMember("isotope", isotope, al);
      }
      if (radical.Size() > 0) {
        atoms.AddMember("radical", radical, al);
      }

      doc.AddMember("atoms", atoms, al);
    }

    rapidjson::Value aid1(rapidjson::kArrayType);
    rapidjson::Value aid2(rapidjson::kArrayType);
    rapidjson::Value order(rapidjson::kArrayType);
    rapidjson::Value annAid1(rapidjson::kArrayType);
    rapidjson::Value annAid2(rapidjson::kArrayType);
    rapidjson::Value annotation(rapidjson::kArrayType);

    // Bonds
    FOR_BONDS_OF_MOL(pbond, pmol) {

      // Order
      int ord = pbond->GetBondOrder();
      if (ord == 0) {
        if (pbond->HasData("type")) {
          // Check to see if a "type" string exists
          OBPairData *typeData = dynamic_cast<OBPairData *>(pbond->GetData("type"));
          const string &orderstring = typeData->GetValue();
          if (orderstring == "dative") {
            ord = 5;
          } else if (orderstring == "complex") {
            ord = 6;
          } else if (orderstring == "ionic") {
            ord = 7;
          }
        }
      }
      aid1.PushBack(rapidjson::Value((int) pbond->GetBeginAtom()->GetId()).Move(), al);
      aid2.PushBack(rapidjson::Value((int) pbond->GetEndAtom()->GetId()).Move(), al);
      order.PushBack(rapidjson::Value(ord).Move(), al);

      // Styles and annotations
      vector<int> annotations;
      if (pConv->IsOption("w", pConv->OUTOPTIONS)) {
        // option w means just use input bond stereo annotations
        if (pbond->IsWedge()) {
          annotations.push_back(5);
        } else if (pbond->IsHash()) {
          annotations.push_back(6);
        } else if (pbond->IsWedgeOrHash()) {
          annotations.push_back(3);
        } else if (pbond->IsCisOrTrans()) {
          annotations.push_back(1);
        }
      } else {
        // No option w means use stereochemistry information
        from_cit = from.find(&*pbond);
        if (from_cit != from.end() && from_cit->second == pbond->GetEndAtom()->GetId()) {
          swap(aid1, aid2);  // Swap start and end atom if necessary
        }
        if (unspec_ctstereo.find(&*pbond) != unspec_ctstereo.end()) {
          annotations.push_back(1);
        }
        if (updown.find(&*pbond) != updown.end()) {
          if (updown[&*pbond] == 1) {
            annotations.push_back(5);
          } else if (updown[&*pbond] == 4) {
            annotations.push_back(3);
          } else if (updown[&*pbond] == 6) {
            annotations.push_back(6);
          }
        }
      }
      if (pbond->IsAromatic()) {
        annotations.push_back(8);
      }
      if (pbond->HasData("style")) {
        AnnotationData *data = dynamic_cast<AnnotationData *>(pbond->GetData("style"));
        vector<string> styles = data->GetGenericValue();
        for (vector<string>::const_iterator i = styles.begin(); i != styles.end(); ++i) {
          string stylestring = *i;
          int style = 255;
          if (stylestring == "dashed") {
            style = 2;
          } else if (stylestring == "dotted") {
            style = 4;
          } else if (stylestring == "arrow") {
            style = 7;
          } else if (stylestring == "resonance") {
            style = 9;
          } else if (stylestring == "bold") {
            style = 10;
          } else if (stylestring == "fischer") {
            style = 11;
          } else if (stylestring == "closeContact") {
            style = 12;
          }
          annotations.push_back(style);
        }
      }
      annotations.erase(unique(annotations.begin(), annotations.end()), annotations.end());
      for (vector<int>::const_iterator i = annotations.begin(); i != annotations.end(); ++i) {
        annAid1.PushBack(rapidjson::Value((int) pbond->GetBeginAtom()->GetId()).Move(), al);
        annAid2.PushBack(rapidjson::Value((int) pbond->GetEndAtom()->GetId()).Move(), al);
        annotation.PushBack(rapidjson::Value(*i).Move(), al);
      }
    }

    // Add bonds to doc
    if (aid1.Size() > 0) {
      rapidjson::Value bonds(rapidjson::kObjectType);
      bonds.AddMember("aid1", aid1, al);
      bonds.AddMember("aid2", aid2, al);
      bonds.AddMember("order", order, al);
      doc.AddMember("bonds", bonds, al);
    }

    // Add coords to doc
    if (doc.HasMember("atoms")) {
      rapidjson::Value conf(rapidjson::kObjectType);
      conf.AddMember("x", xcoords, al);
      conf.AddMember("y", ycoords, al);
      if (pmol->GetDimension() == 3) {
        conf.AddMember("z", zcoords, al);
      }
      if (annotation.Size() > 0) {
        rapidjson::Value style(rapidjson::kObjectType);
        style.AddMember("annotation", annotation, al);
        style.AddMember("aid1", annAid1, al);
        style.AddMember("aid2", annAid2, al);
        conf.AddMember("style", style, al);
      }
      rapidjson::Value conformers(rapidjson::kArrayType);
      conformers.PushBack(conf, al);
      rapidjson::Value coordType(rapidjson::kArrayType);
      if (pmol->GetDimension() == 2) {
        coordType.PushBack(rapidjson::Value(1).Move(), al);
      } else if (pmol->GetDimension() == 3) {
        coordType.PushBack(rapidjson::Value(2).Move(), al);
      }
      rapidjson::Value coord(rapidjson::kObjectType);
      coord.AddMember("type", coordType, al);
      coord.AddMember("aids", rapidjson::Value(doc["atoms"]["aids"], al).Move(), al);  // Copy
      coord.AddMember("conformers", conformers, al);
      rapidjson::Value coords(rapidjson::kArrayType);
      coords.PushBack(coord, al);
      doc.AddMember("coords", coords, al);
    }

    // Stereochemistry
    rapidjson::Value stereo(rapidjson::kArrayType);
    OBStereoFacade facade(pmol);
    FOR_ATOMS_OF_MOL(patom, pmol) {
      if (facade.HasTetrahedralStereo(patom->GetId())) {
        OBTetrahedralStereo::Config config = facade.GetTetrahedralStereo(patom->GetId())->GetConfig();
        rapidjson::Value tet(rapidjson::kObjectType);
        tet.AddMember("type", rapidjson::Value(1).Move(), al);  // "tetrahedral"
        tet.AddMember("center", rapidjson::Value((int) config.center).Move(), al);
        tet.AddMember("top", rapidjson::Value((int) config.from).Move(), al);
        tet.AddMember("below", rapidjson::Value((int) config.refs[0]).Move(), al);
        tet.AddMember("bottom", rapidjson::Value((int) config.refs[1]).Move(), al);
        tet.AddMember("above", rapidjson::Value((int) config.refs[2]).Move(), al);
        if (config.winding == OBStereo::UnknownWinding || !config.specified) {
          tet.AddMember("parity", rapidjson::Value(3).Move(), al);  // "any"
        } else if (config.winding == OBStereo::Clockwise) {
          tet.AddMember("parity", rapidjson::Value(1).Move(), al);  // "clockwise"
        } else if (config.winding == OBStereo::AntiClockwise) {
          tet.AddMember("parity", rapidjson::Value(2).Move(), al);  // "counterclockwise"
          tet.AddMember("bottom", rapidjson::Value((int) config.refs[2]).Move(), al);
          tet.AddMember("above", rapidjson::Value((int) config.refs[1]).Move(), al);
        }
        rapidjson::Value stereoContainer(rapidjson::kObjectType);
        stereoContainer.AddMember("tetrahedral", tet, al);
        stereo.PushBack(stereoContainer, al);
      }
      if (facade.HasSquarePlanarStereo(patom->GetId())) {
        OBSquarePlanarStereo *sqs = facade.GetSquarePlanarStereo(patom->GetId());
        OBSquarePlanarStereo::Config config = sqs->GetConfig();
        rapidjson::Value sq(rapidjson::kObjectType);
        sq.AddMember("center", rapidjson::Value((int) config.center).Move(), al);
        if (config.specified) {
          if (config.shape == OBStereo::ShapeU) {
            sq.AddMember("parity", rapidjson::Value(1).Move(), al);  // "u-shape"
            sq.AddMember("lbelow", rapidjson::Value((int) config.refs[0]).Move(), al);
            sq.AddMember("rbelow", rapidjson::Value((int) config.refs[1]).Move(), al);
            sq.AddMember("rabove", rapidjson::Value((int) config.refs[2]).Move(), al);
            sq.AddMember("labove", rapidjson::Value((int) config.refs[3]).Move(), al);
          } else if (config.shape == OBStereo::ShapeZ) {
            sq.AddMember("parity", rapidjson::Value(2).Move(), al);  // "z-shape"
            sq.AddMember("lbelow", rapidjson::Value((int) config.refs[0]).Move(), al);
            sq.AddMember("rbelow", rapidjson::Value((int) config.refs[1]).Move(), al);
            sq.AddMember("labove", rapidjson::Value((int) config.refs[2]).Move(), al);
            sq.AddMember("rabove", rapidjson::Value((int) config.refs[3]).Move(), al);
          } else if (config.shape == OBStereo::Shape4) {
            sq.AddMember("parity", rapidjson::Value(3).Move(), al);  // "x-shape"
            sq.AddMember("lbelow", rapidjson::Value((int) config.refs[0]).Move(), al);
            sq.AddMember("rabove", rapidjson::Value((int) config.refs[1]).Move(), al);
            sq.AddMember("rbelow", rapidjson::Value((int) config.refs[2]).Move(), al);
            sq.AddMember("labove", rapidjson::Value((int) config.refs[3]).Move(), al);
          }
        } else {
          sq.AddMember("parity", rapidjson::Value(4).Move(), al);  // "any"
          sq.AddMember("lbelow", rapidjson::Value((int) config.refs[0]).Move(), al);
          sq.AddMember("rbelow", rapidjson::Value((int) config.refs[1]).Move(), al);
          sq.AddMember("rabove", rapidjson::Value((int) config.refs[2]).Move(), al);
          sq.AddMember("labove", rapidjson::Value((int) config.refs[3]).Move(), al);
        }
        rapidjson::Value stereoContainer(rapidjson::kObjectType);
        stereoContainer.AddMember("squareplanar", sq, al);
        stereo.PushBack(stereoContainer, al);
      }
    }
    FOR_BONDS_OF_MOL(pbond, pmol) {
      if (facade.HasCisTransStereo(pbond->GetId())) {
        OBCisTransStereo *cts = facade.GetCisTransStereo(pbond->GetId());
        OBCisTransStereo::Config config = cts->GetConfig();
        rapidjson::Value ct(rapidjson::kObjectType);
        ct.AddMember("type", rapidjson::Value(1).Move(), al);  // "planar"
        ct.AddMember("ltop", rapidjson::Value((int) config.refs[0]).Move(), al);
        OBAtom *begin = pmol->GetAtomById(config.begin);
        OBAtom *ltop = pmol->GetAtomById(config.refs[0]);
        if (begin && ltop) {
          if (ltop->IsConnected(begin)) {
            ct.AddMember("left", rapidjson::Value((int) config.begin).Move(), al);
            ct.AddMember("right", rapidjson::Value((int) config.end).Move(), al);
          } else {
            ct.AddMember("left", rapidjson::Value((int) config.end).Move(), al);
            ct.AddMember("right", rapidjson::Value((int) config.begin).Move(), al);
          }
          ct.AddMember("rbottom", rapidjson::Value((int) cts->GetTransRef(config.refs[0])).Move(), al);
          ct.AddMember("lbottom", rapidjson::Value((int) cts->GetCisRef(cts->GetTransRef(config.refs[0]))).Move(), al);
          ct.AddMember("rtop", rapidjson::Value(
              (int) cts->GetTransRef(cts->GetCisRef(cts->GetTransRef(config.refs[0])))).Move(), al);
          if (config.specified) {
            // Open babel is not capable of determining parity? (need CIP rules?)
            ct.AddMember("parity", rapidjson::Value(255).Move(), al);  // "unknown"
          } else {
            ct.AddMember("parity", rapidjson::Value(3).Move(), al);  // "any"
          }
          rapidjson::Value stereoContainer(rapidjson::kObjectType);
          stereoContainer.AddMember("planar", ct, al);
          stereo.PushBack(stereoContainer, al);
        }
      }
    }

    // Add stereo to doc
    if (stereo.Size() > 0) {
      doc.AddMember("stereo", stereo, al);
    }

    // Add charge to doc
    doc.AddMember("charge", rapidjson::Value(pmol->GetTotalCharge()).Move(), al);

    // Create root object and PC_Compounds array if this is the first molecule in the file
    if (!outRoot.IsObject() || !outRoot.HasMember("PC_Compounds")) {
      outRoot.SetObject();
      outRoot.AddMember("PC_Compounds", rapidjson::Value(rapidjson::kArrayType).Move(), al);
    }

    // Add molecule to PC_Compounds array
    outRoot["PC_Compounds"].PushBack(doc, al);

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
