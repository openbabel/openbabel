/**********************************************************************
Copyright (C) 2005 Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "mol.h"
#include <set>
#include <vector>
#include "bitvec.h"
#include "fingerprint.h"
//#ifdef _DEBUG
// #include "stdafx.h"
// #undef AddAtom
//#endif
using namespace std;
namespace OpenBabel
{

/*
Similar to Fabien Fontain's fingerprint class, with a slightly improved
algorithm, but re-written using STL which makes it shorter.
 
A molecule structure is analysed to identify linear fragments of length
from one to MAX_FRAGMENT_SIZE(=7) atoms but single atom fragments of C,N,and O 
are ignored. A fragment is terminated when the atoms form a ring. 

For each of these fragments the atoms, bonding and whether
they constitute a complete ring is recorded and saved in a std::set
so that there is only one of each fragment type. The number of fragments is
further reduced by eliminating chemically identical versions: ones with 
the atoms listed in reverse order and rings listed starting at different atoms.
Each remaining fragment is assigned a hash number from 0 to 1020 which is used
to set a bit in a 1024 bit vector  
*/

void fingerprint2::GetFingerprint(OBMol& mol,OBBitVec& fp) 
{
	Fset fragset;

	//identify fragments starting at every atom
	OBAtom *patom;
	vector<OBNodeBase*>::iterator i;
	for (patom = mol.BeginAtom(i);patom;patom = mol.NextAtom(i))
	{
		if(patom->IsHydrogen()) continue;

		vector<int> levels(mol.NumAtoms());
		vector<int> curfrag;
	
		getFragments(fragset, levels, curfrag, 1, patom, NULL);
	}
	
//	TRACE("%s %d frags before; ",mol.GetTitle(),fragset.size());
	Fset::iterator itr;
	for(itr=fragset.begin();itr!=fragset.end();++itr)
	{
		//Remove some chemically identical fragments
		RemoveDuplicates(fragset,*itr);
		
		//Use hash of fragment to set a bit in the fingerprint
		int hash = CalcHash(*itr);
		fp.SetBitOn(hash);
//		PrintFpt(*itr,hash);
	}
//	TRACE("%d after\n",fragset.size());
}

//////////////////////////////////////////////////////////
void fingerprint2::getFragments(Fset& fragset, vector<int>& levels, vector<int> curfrag, 
					int level, OBAtom* patom, OBBond* pbond)
{
	//Recursive routine to populate fragset
	int bo=0;
	if(pbond)
		bo = pbond->IsAromatic() ? 5 : pbond->GetBO();
	curfrag.push_back(bo);
	curfrag.push_back(patom->GetAtomicNum());
	levels[patom->GetIdx()-1] = level;

	vector<OBEdgeBase*>::iterator itr;
	OBBond *pnewbond;
//	PrintFpt(curfrag,(int)patom);
	for (pnewbond = patom->BeginBond(itr);pnewbond;pnewbond = patom->NextBond(itr))
	{		
		if(pnewbond==pbond) continue; //don't retrace steps
		OBAtom* pnxtat = pnewbond->GetNbrAtom(patom);
		if(pnxtat->IsHydrogen()) continue;

		int atlevel = levels[pnxtat->GetIdx()-1];
		if(atlevel) //ring
		{
			//If complete ring (last bond is back to starting atom) add bond at front
			if(atlevel==1)
				curfrag[0] = bo;
		}
		else //no ring
		{
			if(level<=MAX_FRAGMENT_SIZE)
			{
				//Do the next atom; curfrag is passed by value and hence copied
//				TRACE("level=%d size=%d %p frag[0]=%p\n",level, curfrag.size(),&curfrag, &(curfrag[0])); 
				getFragments(fragset, levels, curfrag, level+1, pnxtat, pnewbond);
			}
		}
	}

	//do not save C,N,O single atom fragments
	if(level>1 || patom->GetAtomicNum()>8  || patom->GetAtomicNum()<6)
		fragset.insert(curfrag); //curfrag ignored if an identical fragment already present
}

///////////////////////////////////////////////////
void fingerprint2::RemoveDuplicates(Fset& fragset, vector<int>& frag)
{
	vector<int> t1(frag); //temporary copy
	if(t1[0])
	{
		//is a complete ring
		int i;
		for(i=0;i<(t1.size()-1)/2;++i)
		{
			//rotate atoms in ring and erase from set if present
			rotate(t1.begin(),t1.begin()+2,t1.end());
			if(t1!=frag)
				fragset.erase(t1);
			
			//reverse the direction around ring and erase from set if present
			vector<int> t2(t1);
			reverse(t2.begin()+1, t2.end());
			if(t2!=frag)
				fragset.erase(t2);
		}
	}
	else
	{
		//not a complete ring
		//reverse the order of the atoms and remove from set if present
		reverse(t1.begin()+1, t1.end()); //(leave 0 at front alone)
		if(t1!=frag)
			fragset.erase(t1);
	}
}
//////////////////////////////////////////////////////////
int fingerprint2::CalcHash(vector<int>& frag)
{
	//Something like... whole of fragment treated as a binary number modulus 1021
	const int MODINT = 108; //2^32 % 1021 
	int i, hash=0;
	for(i=0;i<frag.size();++i)
		hash= (hash*MODINT + (frag[i] % 1021)) % 1021;
	return hash;
}

void fingerprint2::PrintFpt(vector<int>& f, int hash)
{
	int i;
	for(i=0;i<f.size();++i)
//		TRACE("%d ",f[i]);
//	TRACE("<%d>\n",hash);
		cout << f[i] << " ";
	cout << "<" << hash << ">" << endl;
}

/* Structure of a fragment (vector<int>)
   For a complete ring: last atom bonded to first atom
    bo(0)(n), atno(1), bo(1)(2), atno(2), bo(2)(3),...atno(n)

 For the rest, even when stopped by encountering atoms already visited
       0    , atno(1), bo(1)(2), atno(2), bo(2)(3),...atno(n) 
*/

//***************************************************
bool GetFingerprint(OBMol& mol, vector<unsigned int>& vecwords, int nwords, int type)
{
	//Handle different fingerprint types in this simple way until the need is demonstrated
	//for a more sophisticated dynamic method similar to that for formats.
	if(type==0) type=2; //fingerprint2 is the default
	OBBitVec f(1024);
	switch(type)
	{
		case 1:
		{
				//Fabien Fontain's fingerprint
			fingerprint fpt(mol.GetTitle());
			fpt.HashMol(mol);	
			f = fpt.GetFpt();
			break;
		}
		case 2:
		{
				//fingerprint2
			fingerprint2 fing2;
			fing2.GetFingerprint(mol,f);
			break;
		}
		default:
			return false;
	}
	
	if(nwords)
		f.Fold(nwords*SETWORD);
	f.GetWords(vecwords);
	return true;
}

double Tanimoto(const vector<unsigned int> vec1, const vector<unsigned int> vec2) 
{
	//Independent of sizeof(unsigned int)
	if(vec1.size()!=vec2.size())
		return -1; //different number of bits
	double andbits=0, orbits=0;
	int i;
  for (i=0;i<vec1.size();++i)
	{
		int andfp = vec1[i] & vec2[i];
		int orfp = vec1[i] | vec2[i];
		//Count bits
		for(;andfp;andfp=andfp<<1)
			if(andfp<0) ++andbits;
		for(;orfp;orfp=orfp<<1)
			if(orfp<0) ++orbits;
	}
    return((double)andbits/(double)orbits);
}

} //namespace OpenBabel
