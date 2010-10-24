#include <openbabel/query.h>
#include <openbabel/obconversion.h>

using namespace std;

namespace OpenBabel {

  OBQuery* CompileMoleculeQuery(OBMol *mol, const OBBitVec &mask)
  {
    // set all atoms to 1 if the mask is empty
    OBBitVec mask2 = mask;
    if (!mask2.CountBits())
      for (unsigned int i = 0; i < mol->NumAtoms(); ++i)
        mask2.SetBitOn(i + 1);

    OBQuery *query = new OBQuery;
    unsigned int offset = 0;
    std::vector<unsigned int> indexes;
    FOR_ATOMS_OF_MOL (obatom, mol) {
      indexes.push_back(obatom->GetIndex() - offset);
      if (!mask2.BitIsSet(obatom->GetIndex() + 1)) {
        offset++;
        continue;
      }
      query->AddAtom(new OBQueryAtom(obatom->GetAtomicNum(), obatom->IsInRing(), obatom->IsAromatic()));
    }
    FOR_BONDS_OF_MOL (obbond, mol) {
      unsigned int beginIndex = obbond->GetBeginAtom()->GetIndex();
      unsigned int endIndex = obbond->GetEndAtom()->GetIndex();
      if (!mask2.BitIsSet(beginIndex + 1) || !mask2.BitIsSet(endIndex + 1))
        continue;

      query->AddBond(new OBQueryBond(query->GetAtoms()[indexes[beginIndex]], query->GetAtoms()[indexes[endIndex]],
            obbond->GetBondOrder(), obbond->IsAromatic()));
    }

    return query;
  }

  OBQuery* CompileSmilesQuery(const std::string &smiles, const OBBitVec &mask)
  {
    OBConversion conv;
    conv.SetInFormat("smi");
    OBMol mol;
    conv.ReadString(&mol, smiles);
    return CompileMoleculeQuery(&mol, mask);
  }












}
