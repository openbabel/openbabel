/**********************************************************************
Copyright (C) 2022 by Geoffrey R. Hutchison

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

#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/json.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/tetrahedral.h>

using namespace std;
namespace OpenBabel {

class ChemicalJSONFormat : public OBMoleculeFormat {
public:
  ChemicalJSONFormat() { OBConversion::RegisterFormat("cjson", this); }

  virtual const char *Description() {
    return "Chemical JSON\n"
           "The native file format for Avogadro2 and Open Chemistry\n\n"

           "Write Options, e.g. -xv\n"
           " m  minified output formatting, with no line breaks or indents\n"
           " v  verbose output (include default values)\n\n";
  };

  virtual const char *SpecificationURL() { return ""; };

  virtual bool ReadMolecule(OBBase *pOb, OBConversion *pConv);
  virtual bool WriteMolecule(OBBase *pOb, OBConversion *pConv);

private:
  rapidjson::Document inRoot;
  rapidjson::Document outRoot;
  int currentMolIndex;
};

ChemicalJSONFormat theChemicalJSONFormat;

bool ChemicalJSONFormat::ReadMolecule(OBBase *pOb, OBConversion *pConv) {
  OBMol *pmol = pOb->CastAndClear<OBMol>();
  if (pmol == nullptr)
    return false;
  istream &ifs = *pConv->GetInStream();

  if (!ifs.good())
    return false;

  // Parse entire file into memory once, then reuse inRoot for subsequent
  // molecules (It's really tricky to stream json)
  if (ifs.peek() != EOF) {
    rapidjson::IStreamWrapper isw(ifs);
    inRoot.ParseStream(isw);
    if (inRoot.HasParseError()) {
      stringstream msg;
      msg << "JSON parse error at offset " << inRoot.GetErrorOffset() << ": "
          << rapidjson::GetParseError_En(inRoot.GetParseError());
      obErrorLog.ThrowError("ChemicalJSONFormat", msg.str(), obError);
      return false;
    }

    // at this point, inRoot is the entire file
    // so we should be good to just parse and then we're at EOF
  }

  if (!inRoot.IsObject()) {
    obErrorLog.ThrowError("ChemicalJSONFormat",
                          "JSON file should be a single object", obError);
    return false;
  }

  // CJSON has a single object, but check that it's really CJSON
  if (!inRoot.HasMember("chemicalJson")) {
    obErrorLog.ThrowError("ChemicalJSONFormat",
                          "This does not have a chemicalJSON object", obError);
    return false;
  }

  // Atoms
  // should have several sections and we'll check for them
  // coords (usually 3D as an array)
  // elements / "number" as an array
  // formal charges (optional)
  // labels (optional)
  if (!inRoot.HasMember("atoms") || !inRoot["atoms"].IsObject()) {
    obErrorLog.ThrowError("ChemicalJSONFormat", "Atoms are not specified",
                          obError);
    return false;
  }

  // Sanity check for coordinates
  const rapidjson::Value &atoms = inRoot["atoms"];
  if (!atoms.HasMember("coords") || !atoms["coords"].IsObject()) {
    obErrorLog.ThrowError("ChemicalJSONFormat", "Coordinates are not specified",
                          obError);
    return false;
  }
  const rapidjson::Value &coords = atoms["coords"];
  // check for 3D coordinates
  if (!coords.HasMember("3d") || !coords["3d"].IsArray()) {
    obErrorLog.ThrowError("ChemicalJSONFormat",
                          "3D coordinates are not specified", obError);
    return false;
  }

  // Sanity check for elements / numbers
  if (!atoms.HasMember("elements") || !atoms["elements"].IsObject()) {
    obErrorLog.ThrowError("ChemicalJSONFormat", "Elements are not specified",
                          obError);
    return false;
  }
  const rapidjson::Value &elements = atoms["elements"];
  if (!elements.HasMember("number") || !elements["number"].IsArray()) {
    obErrorLog.ThrowError("ChemicalJSONFormat",
                          "Element numbers are not specified", obError);
    return false;
  }
  const rapidjson::Value &elementNumbers = elements["number"];

  // check that we have the right number of coordinates
  if (coords["3d"].Size() != 3 * elementNumbers.Size()) {
    obErrorLog.ThrowError(
        "ChemicalJSONFormat",
        "Number of coordinates does not match number of elements", obError);
    return false;
  }

  pmol->BeginModify();
  pmol->SetDimension(3);
  pmol->ReserveAtoms(elementNumbers.Size());

  // Add atoms
  double x, y, z;
  for (rapidjson::SizeType i = 0; i < elementNumbers.Size(); i++) {
    const rapidjson::Value &element = elementNumbers[i];
    OBAtom *patom = pmol->NewAtom();
    patom->SetAtomicNum(element.GetInt());

    // Set coordinates
    x = coords["3d"][3 * i].GetDouble();
    y = coords["3d"][3 * i + 1].GetDouble();
    z = coords["3d"][3 * i + 2].GetDouble();
    patom->SetVector(x, y, z);
  }

  // Bonds are optional (e.g., crystals)
  if (inRoot.HasMember("bonds") && inRoot["bonds"].IsObject()) {
    const rapidjson::Value &bonds = inRoot["bonds"];
    if (bonds.HasMember("connections") && bonds["connections"].IsObject()) {
      const rapidjson::Value &connections = bonds["connections"];
      if (connections.HasMember("index") && connections["index"].IsArray()) {
        const rapidjson::Value &index = connections["index"];
        // check the orders as well
        if (bonds.HasMember("order") && bonds["order"].IsArray()) {
          const rapidjson::Value &order = bonds["order"];

          // check that we have the right number of bonds
          if (index.Size() != 2 * order.Size()) {
            obErrorLog.ThrowError(
                "ChemicalJSONFormat",
                "Number of bonds does not match number of indices", obError);
            return false;
          }

          // Add bonds
          for (rapidjson::SizeType i = 0; i < order.Size(); i++) {
            auto start = index[2 * i].GetInt() + 1;
            auto end = index[2 * i + 1].GetInt() + 1;
            pmol->AddBond(start, end, order[i].GetInt());
          }
        }
      }
    }
    // todo
  }

  // Automatically determine spin multiplicity for atoms with hydrogens
  // specified
  pmol->AssignSpinMultiplicity();
  pmol->EndModify();

  if (pmol->Has3D()) {
    // Use 3D coordinates to determine stereochemistry
    StereoFrom3D(pmol);
  }

  return true;
}

bool ChemicalJSONFormat::WriteMolecule(OBBase *pOb, OBConversion *pConv) {
  OBMol *pmol = dynamic_cast<OBMol *>(pOb);
  if (pmol == nullptr)
    return false;
  ostream &ofs = *pConv->GetOutStream();

  if (pmol->GetDimension() != 3) {
    obErrorLog.ThrowError("ChemicalJSONFormat",
                          "No 3D coordinates exist. "
                          "To 3D coordinates use --gen3D.",
                          obError);
    return false;
  }

  // Must always pass an allocator when memory may need to be allocated
  rapidjson::Document::AllocatorType &al = outRoot.GetAllocator();

  rapidjson::Value doc(rapidjson::kObjectType); // Root of molecule JSON
  doc.AddMember("chemicalJson", 1, al);
  // doc.AddMember("name", pmol->GetTitle().c_str(), al);

  // Atoms
  rapidjson::Value atoms(rapidjson::kObjectType);
  rapidjson::Value coords3d(rapidjson::kArrayType);
  rapidjson::Value elementNumbers(rapidjson::kArrayType);
  rapidjson::Value formalCharges(rapidjson::kArrayType);
  rapidjson::Value partialCharges(rapidjson::kArrayType);
  rapidjson::Value nmrShifts(rapidjson::kArrayType);

  std::string chargeMethod = "Gasteiger"; // that's the default
  OBPairData *dp = (OBPairData *)pmol->GetData("PartialCharges");
  if (dp != nullptr)
    chargeMethod = dp->GetValue();

  FOR_ATOMS_OF_MOL(patom, pmol) {
    // Add coordinates
    coords3d.PushBack(patom->x(), al);
    coords3d.PushBack(patom->y(), al);
    coords3d.PushBack(patom->z(), al);

    // Add element number
    elementNumbers.PushBack(patom->GetAtomicNum(), al);

    // formal charges
    formalCharges.PushBack(patom->GetFormalCharge(), al);

    // partial charges
    partialCharges.PushBack(patom->GetPartialCharge(), al);

    // check for NMR shifts
    if (patom->HasData("NMR Isotropic Shift"))
      nmrShifts.PushBack(
          rapidjson::StringRef(patom->GetData("NMR Isotropic Shift")->GetValue().c_str()), al);
  }

  // conformers / multiple coordinates
  rapidjson::Value coords(rapidjson::kObjectType);
  coords.AddMember("3d", coords3d, al); // default coords
  if (pmol->NumConformers() > 1) {
    rapidjson::Value conformers(rapidjson::kArrayType);

    for (unsigned int i = 0; i < pmol->NumConformers(); i++) {
      pmol->SetConformer(i);
      rapidjson::Value conformer(rapidjson::kArrayType);
      FOR_ATOMS_OF_MOL(patom, pmol) {
        conformer.PushBack(patom->x(), al);
        conformer.PushBack(patom->y(), al);
        conformer.PushBack(patom->z(), al);
      }
      conformers.PushBack(conformer, al);
    }
    coords.AddMember("3dSets", conformers, al);
  }
  atoms.AddMember("coords", coords, al);

  rapidjson::Value elements(rapidjson::kObjectType);
  elements.AddMember("number", elementNumbers, al);
  atoms.AddMember("elements", elements, al);

  atoms.AddMember("formalCharges", formalCharges, al);
  doc.AddMember("atoms", atoms, al);

  rapidjson::Value charges(rapidjson::kObjectType);
  charges.AddMember(rapidjson::StringRef(chargeMethod.c_str()), partialCharges,
                    al);
  doc.AddMember("partialCharges", charges, al);

  // optionally add the NMR spectra
  // spectra: { "nmr": { "shifts": [1.123, 115.0, 3.75] } }
  rapidjson::Value spectra(rapidjson::kObjectType);
  if (nmrShifts.Size() > 0) {
    rapidjson::Value nmr(rapidjson::kObjectType);
    nmr.AddMember("shifts", nmrShifts, al);
    spectra.AddMember("nmr", nmr, al);
  }

  // Bonds
  if (pmol->NumBonds() > 0) {
    rapidjson::Value bonds(rapidjson::kObjectType);
    rapidjson::Value connectionIdx(rapidjson::kArrayType);
    rapidjson::Value bondOrder(rapidjson::kArrayType);

    FOR_BONDS_OF_MOL(pbond, pmol) {
      connectionIdx.PushBack(pbond->GetBeginAtomIdx() - 1, al);
      connectionIdx.PushBack(pbond->GetEndAtomIdx() - 1, al);
      bondOrder.PushBack(pbond->GetBondOrder(), al);
    }
    rapidjson::Value connections(rapidjson::kObjectType);
    connections.AddMember("index", connectionIdx, al);

    bonds.AddMember("connections", connections, al);
    bonds.AddMember("order", bondOrder, al);
    doc.AddMember("bonds", bonds, al);
  }

  // unit cells
  if (pmol->HasData(OBGenericDataType::UnitCell)) {
    OBUnitCell *uc = (OBUnitCell *)pmol->GetData(OBGenericDataType::UnitCell);
    rapidjson::Value unitCell(rapidjson::kObjectType);
    unitCell.AddMember("a", uc->GetA(), al);
    unitCell.AddMember("b", uc->GetB(), al);
    unitCell.AddMember("c", uc->GetC(), al);
    unitCell.AddMember("alpha", uc->GetAlpha(), al);
    unitCell.AddMember("beta", uc->GetBeta(), al);
    unitCell.AddMember("gamma", uc->GetGamma(), al);

    // also write the cell vectors
    rapidjson::Value cellVectors(rapidjson::kArrayType);
    vector<vector3> obVectors = uc->GetCellVectors();
    for (vector<vector3>::iterator i = obVectors.begin(); i != obVectors.end();
         ++i) {
      cellVectors.PushBack(i->x(), al);
      cellVectors.PushBack(i->y(), al);
      cellVectors.PushBack(i->z(), al);
    }

    doc.AddMember("unitCell", unitCell, al);
  }

  // vibrations
  if (pmol->HasData(OBGenericDataType::VibrationData)) {
    OBVibrationData *vib =
        (OBVibrationData *)pmol->GetData(OBGenericDataType::VibrationData);
    rapidjson::Value vibrations(rapidjson::kObjectType);

    rapidjson::Value frequencies(rapidjson::kArrayType);
    rapidjson::Value modes(rapidjson::kArrayType);
    vector<double> wavenumbers = vib->GetFrequencies();
    unsigned int mode = 1;
    unsigned int modeCount = vib->GetNumberOfFrequencies();
    for (unsigned int i = 0; i < modeCount; i++) {
      frequencies.PushBack(wavenumbers[i], al);
      modes.PushBack(mode++, al);
    }
    vibrations.AddMember("frequencies", frequencies, al);
    vibrations.AddMember("modes", modes, al);

    rapidjson::Value intensities(rapidjson::kArrayType);
    vector<double> intensitiesVec = vib->GetIntensities();
    for (unsigned int i = 0; i < modeCount; i++) {
      intensities.PushBack(intensitiesVec[i], al);
    }
    vibrations.AddMember("intensities", intensities, al);

    rapidjson::Value raman(rapidjson::kArrayType);
    vector<double> ramanVec = vib->GetRamanActivities();
    if (ramanVec.size() > 0) {
      for (unsigned int i = 0; i < modeCount; i++) {
        raman.PushBack(ramanVec[i], al);
      }
      vibrations.AddMember("ramanIntensities", raman, al);
    }

    rapidjson::Value displacements(rapidjson::kArrayType);
    auto lx = vib->GetLx();
    for (unsigned int i = 0; i < modeCount; i++) {
      rapidjson::Value displacement(rapidjson::kArrayType);
      auto obDisp = lx[i]; // this is a vector<vector3>
      for (auto j = obDisp.begin(); j != obDisp.end(); ++j) {
        displacement.PushBack(j->x(), al);
        displacement.PushBack(j->y(), al);
        displacement.PushBack(j->z(), al);
      }

      displacements.PushBack(displacement, al);
    }
    vibrations.AddMember("eigenVectors", displacements, al);

    doc.AddMember("vibrations", vibrations, al);
  }

  // check for electronic spectra (UV/Vis, CD)
  if (pmol->HasData(OBGenericDataType::ElectronicData)) {
    OBElectronicTransitionData *edata =
        (OBElectronicTransitionData *)pmol->GetData(
            OBGenericDataType::ElectronicTransitionData);
    rapidjson::Value electronic(rapidjson::kObjectType);
    rapidjson::Value energies(rapidjson::kArrayType);
    rapidjson::Value intensities(rapidjson::kArrayType);
    // get the energies and intensities
    std::vector<double> wavelengths = edata->GetWavelengths();
    std::vector<double> forces = edata->GetForces();

    // we need to convert the wavelengths to eV
    const double hc = 1239.841984332; // in eV*nm
    for (unsigned int i = 0; i < wavelengths.size(); i++) {
      energies.PushBack(hc / wavelengths[i], al);
      intensities.PushBack(forces[i], al);
    }
    electronic.AddMember("energies", energies, al);
    electronic.AddMember("intensities", intensities, al);

    std::vector<double> rotatoryStrengthsVec =
        edata->GetRotatoryStrengthsLength();
    if (rotatoryStrengthsVec.size() > 0) {
      rapidjson::Value rotatoryStrengths(rapidjson::kArrayType);
      for (unsigned int i = 0; i < rotatoryStrengthsVec.size(); i++) {
        rotatoryStrengths.PushBack(rotatoryStrengthsVec[i], al);
      }
      electronic.AddMember("rotation", rotatoryStrengths, al);
    }
    spectra.AddMember("electronic", electronic, al);
  }

  if (spectra.MemberCount() > 0) {
    doc.AddMember("spectra", spectra, al);
  }

  // @todo
  // residues / chains
  // other spectra
  // other properties

  // Write json to output stream if this is the last molecule in the file
  if (pConv->IsLast()) {
    rapidjson::OStreamWrapper osw(ofs);
    if (pConv->IsOption("m", pConv->OUTOPTIONS)) {
      rapidjson::Writer<rapidjson::OStreamWrapper> writer(osw);
      doc.Accept(writer);
    } else {
      rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(osw);
      writer.SetIndent(' ', 2);
      doc.Accept(writer);
    }
  }
  return true;
}

} // namespace OpenBabel
