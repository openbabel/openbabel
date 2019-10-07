/**********************************************************************
Copyright (C) 2019 by Pat Lorton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/* Populate an OBMol from a Schrödinger Maestro file using MaeParser.
   More information can be found at: https://github.com/schrodinger/maeparser
*/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <iostream>
#include <map>

#include <MaeConstants.hpp>
#include <Reader.hpp>
#include <Writer.hpp>

using namespace std;
using namespace schrodinger::mae;
using boost::dynamic_bitset;
namespace OpenBabel
{

class MAEFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  MAEFormat()
	{
		OBConversion::RegisterFormat("mae", this);
		OBConversion::RegisterFormat("maegz", this);
	}

	virtual const char* Description() override //required
	{
		return
		"Maestro format\n"
		"File format of Schrödinger Software\n";
    };

    //URL where the file format is specified
    virtual const char* SpecificationURL() override
    {
        return "https://github.com/schrodinger/maeparser";
    };

    virtual int SkipObjects(int n, OBConversion* pConv) override;

    ////////////////////////////////////////////////////
    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv) override;
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv) override;

private:
    // TODO:  This will be implemented in maeparser more completely,
    // we should migrate to that when available.
    const map<int, int> atomic_num_to_color = {
        {1, 21}, {6, 2}, {7, 43}, {8, 70}, {9, 8}, {16, 13}, {17, 9} };

    shared_ptr<Block> m_next_mae;
    shared_ptr<Reader> m_reader;
    shared_ptr<Writer> m_writer;

    shared_ptr<IndexedBlock> TranslateAtomBlock(OBMol* pmol);
    shared_ptr<IndexedBlock> TranslateBondBlock(OBMol* pmol);

    void setupReader(OBConversion* pConv);
    void checkEOF(OBConversion* pConv);

    string m_in_filename = "";
    int m_in_location = -1;


};
	////////////////////////////////////////////////////

//Make an instance of the format class
MAEFormat theMAEFormat;

/////////////////////////////////////////////////////////////////

int MAEFormat::SkipObjects(int n, OBConversion* pConv)
{
    // Required for the MaeParser interface, create a shared_ptr w/o
    // memory management
    setupReader(pConv);
    for(int i=0; i<n; i++) {
        m_next_mae = m_reader->next(CT_BLOCK);
        checkEOF(pConv);
        if(m_next_mae==nullptr) {
            return 0;
        }
    }
    return 0;
};


void MAEFormat::setupReader(OBConversion* pConv)
{
    if(m_in_filename == pConv->GetInFilename() && 
            pConv->GetInStream()->tellg() == m_in_location) return;
    m_in_filename = pConv->GetInFilename();

    // Required for the MaeParser interface, create a shared_ptr w/o
    // memory management
    shared_ptr<istream> ifs(shared_ptr<istream>(), pConv->GetInStream());
    m_reader = make_shared<Reader>(ifs);

    m_next_mae = m_reader->next(CT_BLOCK);
}


/* Guilt disclaimer:  This is an ugly hack, but I think required.
 *
 * Since maeparser buffers what it reads, we have to do some cheating
 * around the input stream in order to keep obconversion reading even when
 * the file pointer has reached the end of the file.
 */
void MAEFormat::checkEOF(OBConversion* pConv)
{
    if(m_next_mae == nullptr) {
        // At the end of the data, set the stream there so obconversion
        // stops iterating
        pConv->GetInStream()->setstate(ios::eofbit);
    } else if(pConv->GetInStream()->eof()) {
        // maeparser is done reading/buffering, but has data left to process
        // (additional molecules in its buffer), so move the input stream away
        // from the end and reset its flags
        pConv->GetInStream()->putback(1);
        pConv->GetInStream()->clear();
    }

    // Keep track of the last position we're at
    m_in_location = pConv->GetInStream()->tellg();

    return;
}


bool MAEFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
        return false;

    setupReader(pConv);

    pmol->BeginModify();
    pmol->SetDimension(3);
    pmol->SetTitle(m_next_mae->getStringProperty(CT_TITLE).c_str());

    const auto atom_data = m_next_mae->getIndexedBlock(ATOM_BLOCK);
    // All atoms are gauranteed to have these three field names:
    const auto atomic_numbers = atom_data->getIntProperty(ATOM_ATOMIC_NUM);
    const auto formal_charges = atom_data->getIntProperty(ATOM_FORMAL_CHARGE);
    const auto xs = atom_data->getRealProperty(ATOM_X_COORD);
    const auto ys = atom_data->getRealProperty(ATOM_Y_COORD);
    const auto zs = atom_data->getRealProperty(ATOM_Z_COORD);
    const auto natoms = atomic_numbers->size();

    pmol->ReserveAtoms(natoms);
    // atomic numbers, and x, y, and z coordinates
    for (size_t i = 0; i < natoms; ++i) {
        OBAtom* patom = pmol->NewAtom();
        patom->SetVector(xs->at(i), ys->at(i), zs->at(i));
        patom->SetAtomicNum(atomic_numbers->at(i));
        patom->SetFormalCharge(formal_charges->at(i));
    }

    const auto bond_data = m_next_mae->getIndexedBlock(BOND_BLOCK);
    // All bonds are gauranteed to have these three field names:
    auto bond_atom_1s = bond_data->getIntProperty(BOND_ATOM_1);
    auto bond_atom_2s = bond_data->getIntProperty(BOND_ATOM_2);
    auto orders = bond_data->getIntProperty(BOND_ORDER);
    const auto bond_count = bond_atom_1s->size();

    for (size_t i = 0; i < bond_count; ++i) {
        // Atom indices in the bond data structure are 1 indexed
        const auto bond_atom_1 = bond_atom_1s->at(i);
        const auto bond_atom_2 = bond_atom_2s->at(i);
        // Bonds may be duplicated in MAE format
        if(bond_atom_1 > bond_atom_2) continue;
        const auto order = orders->at(i);
        const unsigned int flag = 0; // Need to do work here around stereo/kekule
        if (!pmol->AddBond(bond_atom_1, bond_atom_2, order, flag)) {
            return false;
        }
    }

    pmol->EndModify();

    m_next_mae = m_reader->next(CT_BLOCK);
    checkEOF(pConv);

    return true;
}

static void addIntProp(string name, vector<int> values,
        shared_ptr<IndexedBlock>& block)
{
    auto prop = make_shared<IndexedProperty<int> >(values);
    block->setIntProperty(name, prop);
}

static void addRealProp(string name, vector<double> values,
        shared_ptr<IndexedBlock>& block)
{
    auto prop = make_shared<IndexedProperty<double> >(values);
    block->setRealProperty(name, prop);
}


////////////////////////////////////////////////////////////////
shared_ptr<IndexedBlock> MAEFormat::TranslateAtomBlock(OBMol* pmol)
{
    auto atom_block = make_shared<IndexedBlock>(ATOM_BLOCK);

    const auto num_atoms = pmol->NumAtoms();
    // Set up a real property
    vector<double> x, y, z;
    x.resize(num_atoms);
    y.resize(num_atoms);
    z.resize(num_atoms);

    vector<int> atomic_num, formal_charge, mmod_type, color;
    atomic_num.resize(num_atoms);
    formal_charge.resize(num_atoms);
    mmod_type.resize(num_atoms);
    color.resize(num_atoms);

    for (unsigned int i=0; i<num_atoms; i++) {
        const auto atom = pmol->GetAtom(i+1); // 1-based index
        x[i] = atom->x();
        y[i] = atom->y();
        z[i] = atom->z();

        atomic_num[i] = atom->GetAtomicNum();
        formal_charge[i] = atom->GetFormalCharge();
        mmod_type[i] = 62;

        auto pos = atomic_num_to_color.find(atomic_num[i]);
        if (pos == atomic_num_to_color.end()) {
            color[i] = 2;
        } else {
            color[i] = pos->second;
        }
    }

    addRealProp(ATOM_X_COORD, x, atom_block);
    addRealProp(ATOM_Y_COORD, y, atom_block);
    addRealProp(ATOM_Z_COORD, z, atom_block);

    addIntProp(ATOM_ATOMIC_NUM, atomic_num, atom_block);
    addIntProp(ATOM_FORMAL_CHARGE, formal_charge, atom_block);
    // Future versions of maeparaser will have const definitions for these
    addIntProp("i_m_mmod_type", mmod_type, atom_block);
    addIntProp("i_m_color", color, atom_block);

    return atom_block;
}

shared_ptr<IndexedBlock> MAEFormat::TranslateBondBlock(OBMol* pmol)
{
    auto bond_block = make_shared<IndexedBlock>(BOND_BLOCK);

    vector<int> from, to, order;

    OBAtom *nbr;
    OBBond *bond;
    vector<OBBond*>::iterator j;
    int bondline = 0;
    vector<int> zbos;

    const auto num_atoms = pmol->NumAtoms();
    for (unsigned int i=0; i<num_atoms; i++) {
        const auto atom = pmol->GetAtom(i+1); // 1-based index
        for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j)) {
            bond = (OBBond*) *j;
            from.push_back(atom->GetIdx());
            to.push_back(nbr->GetIdx());
            order.push_back(bond->GetBondOrder());
        }
    }
    addIntProp(BOND_ATOM_1, from, bond_block);
    addIntProp(BOND_ATOM_2, to, bond_block);
    addIntProp(BOND_ORDER, order, bond_block);

    return bond_block;
}

bool MAEFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    // The Writer automatically writes the format block at instantiation, so
    // must use a single writer for all writing
    if(pConv->GetOutputIndex()<=1) {
        // Required for the MaeParser interface, create a shared_ptr w/o
        // memory management
        shared_ptr<ostream> ofs(shared_ptr<ostream>(), pConv->GetOutStream());
        m_writer = make_shared<Writer>(ofs);
    }

    /** Write the representation of the OBMol molecule to the output stream **/
    auto mae_block = make_shared<Block>(CT_BLOCK);
    mae_block->setStringProperty(CT_TITLE, pmol->GetTitle());

    auto atom_block = TranslateAtomBlock(pmol);
    auto bond_block = TranslateBondBlock(pmol);

    auto ibm = make_shared<IndexedBlockMap>();
    ibm->addIndexedBlock(atom_block->getName(), atom_block);
    ibm->addIndexedBlock(bond_block->getName(), bond_block);
    mae_block->setIndexedBlockMap(ibm);

    m_writer->write(mae_block);

    return true; //or false to stop converting
}

} //namespace OpenBabel
