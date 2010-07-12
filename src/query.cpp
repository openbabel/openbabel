#include <openbabel/query.h>
#include <openbabel/obconversion.h>

using namespace std;

namespace OpenBabel {
 
  OBQuery* CompileMoleculeQuery(OBMol *mol)
  {
    OBQuery *query = new OBQuery;
    FOR_ATOMS_OF_MOL (obatom, mol) {
      query->AddAtom(new OBQueryAtom(obatom->GetAtomicNum()));
    }
    FOR_BONDS_OF_MOL (obbond, mol) {
      unsigned int beginIndex = obbond->GetBeginAtom()->GetIndex();
      unsigned int endIndex = obbond->GetEndAtom()->GetIndex();
      query->AddBond(new OBQueryBond(query->GetAtoms()[beginIndex], query->GetAtoms()[endIndex],
            obbond->GetBondOrder(), obbond->IsAromatic()));
    }

    return query;  
  }

  OBQuery* CompileSmilesQuery(const std::string &smiles)
  {
    OBConversion conv;
    conv.SetInFormat("smi");
    OBMol mol;
    conv.ReadString(&mol, smiles);
    return CompileMoleculeQuery(&mol);  
  }












}
