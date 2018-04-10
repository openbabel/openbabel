namespace OpenBabel
{
  class OBReactionFacade
  {
  public:
    OBReactionFacade(OBMol *mol): _mol(mol), roles((char*)0)
    {
    };
    ~OBReactionFacade()
    {
      free(roles);
    }
    unsigned int GetRole(OBAtom *atom)
    {
      if (roles==(char*)0)
        AssignRoles(); // lazily assign them
      return roles[atom->GetIdx()];
    }
  private:
    OBMol* _mol;
    char* roles;
    void AssignRoles()
    {
      unsigned int size = _mol->NumAtoms() + 1;
      // Note to self: consider alloca
      roles = (char*)malloc(size * sizeof(char));
      FOR_ATOMS_OF_MOL(atom, _mol) {
        unsigned int rxnrole = 0; // reactant/unspecified
        OBGenericData *data = atom->GetData("rxnrole");
        if (data) {
          OBPairInteger *pi = (OBPairInteger*)data;
          rxnrole = pi->GetGenericValue();
        }
        roles[atom->GetIdx()] = rxnrole;
      }
    }
  };
}
