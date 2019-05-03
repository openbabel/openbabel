/**********************************************************************
obmolecformat.cpp - Implementation of subclass of OBFormat for conversion of OBMol.

Copyright (C) 2005 Chris Morley

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
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#ifdef HAVE_SHARED_POINTER
  #include <openbabel/reaction.h>
#endif

#include <algorithm>

using namespace std;
namespace OpenBabel
{
  bool OBMoleculeFormat::OptionsRegistered=false;
  std::map<std::string, OBMol*> OBMoleculeFormat::IMols;
  OBMol* OBMoleculeFormat::_jmol;
  std::vector<OBMol> OBMoleculeFormat::MolArray;
  bool OBMoleculeFormat::StoredMolsReady=false;

  bool OBMoleculeFormat::ReadChemObjectImpl(OBConversion* pConv, OBFormat* pFormat)
  {
    std::istream *ifs = pConv->GetInStream();
    if (!ifs || !ifs->good())
      return false;

    OBMol* pmol = new OBMol;

    if(pConv->IsOption("C",OBConversion::GENOPTIONS))
      return DeferMolOutput(pmol, pConv, pFormat);

    bool ret=true;
   if(pConv->IsOption("separate",OBConversion::GENOPTIONS))
   {
     //On first call, separate molecule and put fragments in MolArray.
     //On subsequent calls, remove a fragment from MolArray and send it for writing
     //Done this way so that each fragment can be written to its own file (with -m option)
     if(!StoredMolsReady)
     {
       while(ret) //do all the molecules in the file
       {
         ret = pFormat->ReadMolecule(pmol,pConv);

         if(ret && (pmol->NumAtoms() > 0 || (pFormat->Flags()&ZEROATOMSOK)))
         {
           vector<OBMol> SepArray = pmol->Separate(); //use un-transformed molecule
           //Add an appropriate title to each fragment
           if(SepArray.size()>1)
             for (unsigned int i=0; i<SepArray.size(); ++i)
             {
               stringstream ss;
               ss << pmol->GetTitle() << '#' << i+1;
               string title = ss.str();
               SepArray[i].SetTitle(title);
             }
           else
              SepArray[0].SetTitle(pmol->GetTitle());

           copy(SepArray.begin(),SepArray.end(),back_inserter(MolArray));
         }
       }
       reverse(MolArray.begin(),MolArray.end());
       StoredMolsReady = true;
       //Clear the flags of the input stream(which may have found eof) to ensure will
       //try to read anothe molecule and allow the stored ones to be sent for output.
       pConv->GetInStream()->clear();
     }

     if(MolArray.empty()) //normal end of fragments
       ret =false;
     else
     {
       // Copying is needed because the OBMol passed to AddChemObject will be deleted.
       // The OBMol in the vector is deleted here.
       OBMol* pMolCopy = new OBMol( MolArray.back());
       MolArray.pop_back();
       ret = pConv->AddChemObject(
           pMolCopy->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS), pConv))!=0;
     }
     if(!ret)
       StoredMolsReady = false;

     delete pmol;
     return ret;
   }

    ret=pFormat->ReadMolecule(pmol,pConv);

    OBMol* ptmol = NULL;
    //Molecule is valid if it has some atoms
    //or it represents a reaction
    //or the format allows zero-atom molecules and it has a title or properties
    if(ret && (pmol->NumAtoms() > 0 
      || pmol->IsReaction()
      || (pFormat->Flags()&ZEROATOMSOK && (*pmol->GetTitle() || pmol->HasData(1)))))
    {
      ptmol = static_cast<OBMol*>(pmol->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS),pConv));
      if(ptmol && (pConv->IsOption("j",OBConversion::GENOPTIONS)
                || pConv->IsOption("join",OBConversion::GENOPTIONS)))
      {
        //With j option, accumulate all mols in one stored in this class
        if(pConv->IsFirstInput())
          _jmol = new OBMol;
        pConv->AddChemObject(_jmol);
        //will be discarded in WriteChemObjectImpl until the last input mol. This complication
        //is needed to allow joined molecules to be from different files. pOb1 in AddChem Object
        //is zeroed at the end of a file and _jmol is in danger of not being output.
        *_jmol += *ptmol;
        delete ptmol;
        return true;
      }
    }
    else
      delete pmol;

    // Normal operation - send molecule to be written
    ret = ret && (pConv->AddChemObject(ptmol)!=0); //success of both writing and reading
    return ret;
  }

  bool OBMoleculeFormat::WriteChemObjectImpl(OBConversion* pConv, OBFormat* pFormat)
  {
    if(pConv->IsOption("C",OBConversion::GENOPTIONS))
      return OutputDeferredMols(pConv);
    if(pConv->IsOption("j",OBConversion::GENOPTIONS)
        || pConv->IsOption("join",OBConversion::GENOPTIONS))
      {
        //arrives here at the end of a file
        if(!pConv->IsLast())
          return true;
        bool ret=pFormat->WriteMolecule(_jmol,pConv);
        pConv->SetOutputIndex(1);
        delete _jmol;
        return ret;
      }


    //Retrieve the target OBMol
    OBBase* pOb = pConv->GetChemObject();

    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    bool ret=false;
    if(pmol)
      {
        if(pmol->NumAtoms()==0)
          {
            std::string auditMsg = "OpenBabel::Molecule ";
            auditMsg += pmol->GetTitle();
            auditMsg += " has 0 atoms";
            obErrorLog.ThrowError(__FUNCTION__,
                                  auditMsg,
                                  obInfo);
          }
        ret=true;

        ret = DoOutputOptions(pOb, pConv);

        if(ret)
          ret = pFormat->WriteMolecule(pmol,pConv);
    }

#ifdef HAVE_SHARED_POINTER
    //If sent a OBReaction* (rather than a OBMol*) output the consituent molecules
    OBReaction* pReact = dynamic_cast<OBReaction*> (pOb);
    if(pReact)
      ret = OutputMolsFromReaction(pReact, pConv, pFormat);
#endif
    delete pOb;
    return ret;
  }

  bool OBMoleculeFormat::DoOutputOptions(OBBase* pOb, OBConversion* pConv)
  {
    if(pConv->IsOption("addoutindex", OBConversion::GENOPTIONS)) {
      stringstream ss;
      ss << pOb->GetTitle() << " " << pConv->GetOutputIndex();
      pOb->SetTitle(ss.str().c_str());
    }

    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    if(pmol) {
      if(pConv->IsOption("writeconformers", OBConversion::GENOPTIONS)) {
        //The last conformer is written in the calling function
        int c = 0;
        for (; c < pmol->NumConformers()-1; ++c) {
          pmol->SetConformer(c);
          if(!pConv->GetOutFormat()->WriteMolecule(pmol, pConv))
            break;
        }
        pmol->SetConformer(c);
      }
    }
    return true;
  }

  /*! Instead of sending molecules for output via AddChemObject(), they are
    saved in here in OBMoleculeFormat or discarded. By default they are
    saved only if they are in the first input file. Parts of subsequent
    molecules, such as chemical structure, coordinates and OBGenericData
    can replace the parts in molecules with the same title that have already
    been stored, subject to a set of rules. After all input files have been
    read, the stored molecules (possibly now having augmented properties) are
    sent to the output format.

    Is a static function with *this as parameter so that it can be called from other
    format classes like XMLMoleculeFormat which are not derived from OBMoleculeFormat.
  */
  bool OBMoleculeFormat::DeferMolOutput(OBMol* pmol, OBConversion* pConv, OBFormat* pF )
  {
    static bool IsFirstFile;
    bool OnlyMolsInFirstFile=true;

    if(pConv->IsFirstInput())
      {
        IsFirstFile=true;
        IMols.clear();
        pConv->AddOption("OutputAtEnd", OBConversion::GENOPTIONS);
      }
    else
      {
        if((std::streamoff)pConv->GetInStream()->tellg()<=0)
          IsFirstFile=false;//File has changed
      }

    if (!pF->ReadMolecule(pmol,pConv))
      {
        delete pmol;
        return false;
      }
    const char* ptitle = pmol->GetTitle();
    if(*ptitle==0)
      obErrorLog.ThrowError(__FUNCTION__, "Molecule with no title ignored", obWarning);
    else
      {
        string title(ptitle);
        string::size_type pos = title.find_first_of("\t\r\n"); //some title have other data appended
        if(pos!=string::npos)
          title.erase(pos);

        map<std::string, OBMol*>::iterator itr;
        itr = IMols.find(title);
        if(itr!=IMols.end())
          {
            //Molecule with the same title has been input previously: update it
            OBMol* pNewMol = MakeCombinedMolecule(itr->second, pmol);
            if(pNewMol)
              {
                delete itr->second;
                IMols[title] = pNewMol;
              }
            else
              {
                //error: cleanup and return false
                delete pmol;
                return DeleteDeferredMols();
              }
          }
        else
          {
            //Molecule not already saved in IMols: save it if in first file
            if(!OnlyMolsInFirstFile || IsFirstFile)
              {
                IMols[title] = pmol;
                return true; //don't delete pmol
              }
          }
      }
    delete pmol;
    return true;
  }

  /*! Makes a new OBMol on the heap by combining two molecules according to the rule below.
    If both have OBGenericData of the same type, or OBPairData with the
    same attribute,  the version from pFirst is used.
    Returns a pointer to a new OBMol which will need deleting by the calling program
    (probably by being sent to an output format).
    If the molecules cannot be regarded as being the same structure a NULL
    pointer is returned and an error message logged.

    pFirst and pSecond and the objects they point to are not changed. (const
    modifiers difficult because class OBMol not designed appropriately)

    Combining molecules: rules for each of the three parts
    Title:
    Use the title of pFirst unless it has none, when use that of pSecond.
    Warning if neither molecule has a title.

    Structure
    - a structure with atoms replaces one with no atoms
    - a structure with bonds replaces one with no bonds,
    provided the formula is the same, else an error.
    - structures with atoms and bonds are compared by InChI; error if not the same.
    - a structure with 3D coordinates replaces one with 2D coordinates
    - a structure with 2D coordinates replace one with 0D coordinates

    OBGenericData
    OBPairData
  */
  OBMol* OBMoleculeFormat::MakeCombinedMolecule(OBMol* pFirst, OBMol* pSecond)
  {
    //Decide on which OBMol provides the new title
    string title("No title");
    if(*pFirst->GetTitle()!=0)
      title = pFirst->GetTitle();
    else
      {
        if(*pSecond->GetTitle()!=0)
          title = pSecond->GetTitle();
        else
          obErrorLog.ThrowError(__FUNCTION__,"Combined molecule has no title", obWarning);
      }

    //Decide on which OBMol provides the new structure
    bool swap=false;
    if(pFirst->NumAtoms()==0 && pSecond->NumAtoms()!=0)
      swap=true;
    else if(pSecond->NumAtoms()!=0)
      {
        if(pFirst->GetSpacedFormula()!=pSecond->GetSpacedFormula())
          {
            obErrorLog.ThrowError(__FUNCTION__,
                                  "Molecules with name = " + title + " have different formula",obError);
            return NULL;
          }
        else
          {
            if(pSecond->NumBonds()!=0 && pFirst->NumBonds()==0)
              swap=true;
            else
              {
                //Compare by inchi; error if different NOT YET IMPLEMENTED
                //Use the one with the higher dimension
                if(pSecond->GetDimension() > pFirst->GetDimension())
                  swap=true;
              }
          }
      }

    OBMol* pNewMol = new OBMol;
    pNewMol->SetTitle(title);

    OBMol* pMain = swap ? pSecond : pFirst;
    OBMol* pOther = swap ? pFirst : pSecond;

    *pNewMol = *pMain; //Now copies all data

    //Copy some OBGenericData from the OBMol which did not provide the structure
    vector<OBGenericData*>::iterator igd;
    for(igd=pOther->BeginData();igd!=pOther->EndData();++igd)
      {
        //copy only if not already data of the same type from molecule already copied to pNewMol
        unsigned datatype = (*igd)->GetDataType();
        OBGenericData* pData = pNewMol->GetData(datatype);
        if(datatype==OBGenericDataType::PairData)
          {
            if(pData->GetAttribute() == (*igd)->GetAttribute())
              continue;
          }
        else if(pNewMol->GetData(datatype)!=NULL)
          continue;

        OBGenericData* pCopiedData = (*igd)->Clone(pNewMol);
        pNewMol->SetData(pCopiedData);
      }
    return pNewMol;
  }

  bool OBMoleculeFormat::OutputDeferredMols(OBConversion* pConv)
  {
    std::map<std::string, OBMol*>::iterator itr, lastitr;
    bool ret=false;
    int i=1;
    lastitr = IMols.end();
    --lastitr;
    pConv->SetOneObjectOnly(false);
    for(itr=IMols.begin();itr!=IMols.end();++itr,++i)
      {
        if(!(itr->second)->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS), pConv))
          continue;
        pConv->SetOutputIndex(i);
        if(itr==lastitr)
          pConv->SetOneObjectOnly(); //to set IsLast

        ret = pConv->GetOutFormat()->WriteMolecule(itr->second, pConv);

        delete itr->second; //always delete OBMol object
        itr->second = NULL; // so can be deleted in DeleteDeferredMols()
        if (!ret) break;
      }
    DeleteDeferredMols();//cleans up in case there have been errors
    return ret;
  }

  bool OBMoleculeFormat::DeleteDeferredMols()
  {
    //Empties IMols, deteting the OBMol objects whose pointers are stored there
    std::map<std::string, OBMol*>::iterator itr;
    for(itr=IMols.begin();itr!=IMols.end();++itr)
      {
        delete itr->second; //usually NULL
      }
    IMols.clear();
    return false;
  }

  ///////////////////////////////////////////////////////////////////
#ifdef HAVE_SHARED_POINTER
  bool OBMoleculeFormat::OutputMolsFromReaction
    (OBReaction* pReact, OBConversion* pConv, OBFormat* pFormat)
  {
    //Output all the constituent molecules of the reaction

    //Collect the molecules first, just for convenience
    vector<obsharedptr<OBMol> > mols;
    for(int i=0;i<pReact->NumReactants();i++)
      mols.push_back(pReact->GetReactant(i));
    for(int i=0;i<pReact->NumProducts();i++)
      mols.push_back(pReact->GetProduct(i));
    for (int i = 0; i<pReact->NumAgents(); i++)
      mols.push_back(pReact->GetAgent(i));

    if(pReact->GetTransitionState())
      mols.push_back(pReact->GetTransitionState());

    pConv->SetOutputIndex(pConv->GetOutputIndex() - 1); // The OBReaction object is not output
    if((pFormat->Flags() & WRITEONEONLY) && mols.size()>1)
    {
      stringstream ss;
      ss << "There are " << mols.size() << " molecules to be output,"
         << "but this format is for single molecules only";
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
      mols.resize(1);
    }
    bool ok = true;
    for(unsigned int i=0;i<mols.size() && ok;++i)
    {
      if(mols[i])
      {
        //Have to do set these manually because not using "Convert" interface
        pConv->SetLast(i==mols.size()-1);
        pConv->SetOutputIndex(pConv->GetOutputIndex()+1);
        ok = pFormat->WriteMolecule(
          mols[i]->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS), pConv),pConv);
      }
    }
    return ok;
  }
#endif
  //////////////////////////////////////////////////////////////////
  /** Attempts to read the index file datafilename.obindx successively
      from the following directories:
      - the current directory
      - that in the environment variable BABEL_DATADIR or in the macro BABEL_DATADIR
      if the environment variable is not set
      - in a subdirectory of the BABEL_DATADIR directory with the version of OpenBabel as its name
      An index of type NameIndexType is then constructed. NameIndexType is defined
      in obmolecformat.h and may be a std::tr1::unordered_map (a hash_map) or std::map.
      In any case it is searched by
      @code
      NameIndexType::iterator itr = index.find(molecule_name);
      if(itr!=index.end())
      unsigned pos_in_datafile = itr->second;
      @endcode
      pos_in_datafile is used as a parameter in seekg() to read from the datafile

      If no index is found, it is constructed from the datafile by reading all of
      it using the format pInFormat, and written to the directory containing the datafile.
      This means that this function can be used without worrying whether there is an index.
      It will be slow to execute the first time, but subsequent uses get the speed benefit
      of indexed access to the datafile.

      The serialization and de-serialization of the NameIndexType is entirely in
      this routine and could possibly be improved. Currently re-hashing is done
      every time the index is read.
  **/

  bool OBMoleculeFormat::ReadNameIndex(NameIndexType& index,
                                       const string& datafilename, OBFormat* pInFormat)
  {
    struct headertype
    {
      char filename[256];
      size_t size;
    } header;

    NameIndexType::iterator itr;

    ifstream indexstream;
    OpenDatafile(indexstream, datafilename + ".obindx");
    if(!indexstream)
      {
        //Need to prepare the index
        ifstream datastream;
        string datafilepath = OpenDatafile(datastream, datafilename);
        if(!datastream)
          {
            obErrorLog.ThrowError(__FUNCTION__,
                                  datafilename + " was not found or could not be opened",  obError);
            return false;
          }

        OBConversion Conv(&datastream,NULL);
        Conv.SetInFormat(pInFormat);
        OBMol mol;
        streampos pos;
        while(Conv.Read(&mol))
          {
            string name = mol.GetTitle();
            if(!name.empty())
              index.insert(make_pair(name, pos));
            mol.Clear();
            pos = datastream.tellg();
          }
        obErrorLog.ThrowError(__FUNCTION__,
                              "Prepared an index for " + datafilepath, obAuditMsg);
        //Save index to file
        ofstream dofs((datafilepath + ".obindx").c_str(), ios_base::out|ios_base::binary);
        if(!dofs) return false;

        strncpy(header.filename,datafilename.c_str(), sizeof(header.filename));
        header.filename[sizeof(header.filename) - 1] = '\0';
        header.size = index.size();
        dofs.write((const char*)&header, sizeof(headertype));

        for(itr=index.begin();itr!=index.end();++itr)
          {
            //#chars; chars;  ofset(4bytes).
            const char n = static_cast<char> (itr->first.size());
            dofs.put(n);
            dofs.write(itr->first.c_str(),n);
            dofs.write((const char*)&itr->second,sizeof(unsigned));
          }
      }
    else
      {
        //Read index data from file and put into hash_map
        indexstream.read((char*)&header,sizeof(headertype));
        itr=index.begin(); // for hint
        for(unsigned int i=0;i<header.size;++i)
          {
            char len;
            indexstream.get(len);
            string title(len, 0);
            unsigned pos;
            indexstream.read(&title[0],len);
            indexstream.read((char*)&pos,sizeof(unsigned));
            index.insert(itr, make_pair(title,pos));
          }
      }
    return true;
  }

} //namespace OpenBabel

//! \file obmolecformat.cpp
//! \brief Subclass of OBFormat for conversion of OBMol.
