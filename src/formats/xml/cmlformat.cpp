/**********************************************************************
Copyright (C) 2005 by Chris Morley
Some portions Copyright (C) 2006 by Geoffrey R. Hutchison

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

#include <openbabel/math/matrix3x3.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/kinetics.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/obfunctions.h>
#include <openbabel/xml.h>
#include <float.h>
#ifdef HAVE_SHARED_POINTER
  #include <openbabel/reaction.h>
#endif


#ifdef WIN32
#pragma warning (disable : 4800)
#endif
#include <sstream>

using namespace std;
namespace OpenBabel
{

  class CMLFormat : public XMLMoleculeFormat
  {
  private:
    const char* CML1NamespaceURI()const
      {return "http://cml.sourceforge.net/schema/cmlCore/HTMLDOCS/cmlCore.pdf";}
    const char* CML2NamespaceURI()const{return "http://www.xml-cml.org/schema/cml2/core";}

  public:
    //Constuctor used on startup which registers this format type ID
    CMLFormat()
    {
      OBConversion::RegisterFormat("cml", this, "chemical/x-cml");
      OBConversion::RegisterFormat("mrv", this); //subset of Marvin only
			OBConversion::RegisterOptionParam("1", this);
			OBConversion::RegisterOptionParam("a", this);
			OBConversion::RegisterOptionParam("N", this, 1);
			OBConversion::RegisterOptionParam("m", this);
			OBConversion::RegisterOptionParam("x", this);
			OBConversion::RegisterOptionParam("c", this);
			OBConversion::RegisterOptionParam("p", this);
			OBConversion::RegisterOptionParam("2", this, 0, OBConversion::INOPTIONS);

			XMLConversion::RegisterXMLFormat(this, true);	//this is the default XLMformat
			XMLConversion::RegisterXMLFormat(this, false,CML1NamespaceURI());//CML1 also
			XMLConversion::RegisterXMLFormat(this, false,CML2NamespaceURI());//Old CML2 also
    }
    virtual const char* NamespaceURI()const{return "http://www.xml-cml.org/schema";}

    virtual const char* Description()
    {
      return
        "Chemical Markup Language\n"
        "An XML format for interchange of chemical information.\n\n"

        "This format writes and reads CML XML files. To write CML1 format rather than\n"
        "the default CML2, use the ``-x1`` option. To write the array form use ``-xa``\n"
        "and to specify all hydrogens using the hydrogenCount attribute on atoms use\n"
        "``-xh``.\n\n"

        "Crystal structures are written using the <crystal>, <xfract> (,...etc.)\n"
        "elements if the OBMol has a OBGenericDataType::UnitCell data.\n\n"

        "All these forms are handled transparently during reading. Only a subset of\n"
        "CML elements and attributes are recognised, but these include most of those\n"
        "which define chemical structure, see below.\n\n"

        "The following are read:\n\n"

        "- Elements:\n\n"
        "  - molecule, atomArray, atom, bondArray, bond, atomParity, bondStereo\n"
        "  - name, formula, crystal, scalar (contains crystal data)\n"
        "  - string, stringArray, integer, integerArray, float floatArray, builtin\n\n"

        "- Attributes:\n\n"
        "  - On <molecule>: id, title, ref(in CMLReact)\n"
        "  - On <atom>: id, atomId, atomID, elementType, x2, y2, x3, y3, z3, xy2, xyz3,\n"
        "    xFract, yFract, zFract, xyzFract, hydrogenCount, formalCharge, isotope,\n"
        "    isotopeNumber, spinMultiplicity, radical(from Marvin),\n"
        "    atomRefs4 (for atomParity)\n"
        "  - On <bond>: atomRefs2, order, CML1: atomRef, atomRef1, atomRef2\n\n"

        "Atom classes are also read and written. This is done using a specially\n"
        "formed atom id. When reading, if the atom id is of the form aN_M (where\n"
        "N and M are positive integers), then M is interpreted as the atom class.\n"
        "Such atom ids are automatically generated when writing an atom with an\n"
        "atom class.\n\n"

        "Write Options for CML: -x[flags] (e.g. -x1ac)\n"
        "  1  write CML1 (rather than CML2)\n"
        "  a  write array format for atoms and bonds\n"
        "  A  write aromatic bonds as such, not Kekule form\n"
        "  m  write metadata\n"
        "  x  omit XML and namespace declarations\n"
        "  c  continuous output: no formatting\n"
        "  p  write properties\n"
        "  N<prefix> add namespace prefix to elements\n\n"

        "Read Options, e.g. -a2\n"
        "  2  read 2D rather than 3D coordinates if both provided\n\n"

        "In the absence of hydrogenCount and any explicit hydrogen on\n"
        "an atom, implicit hydrogen is assumed to be present appropriate\n"
        "to the radical or spinMultiplicity attributes on the atom or\n"
        "its normal valency if they are not present.\n\n"

        "The XML formats require the XML text to be well formed but\n"
        "generally interpret it fairly tolerantly. Unrecognised elements\n"
        "and attributes are ignored and there are rather few error messages\n"
        "when any required structures are not found. This laxity allows, for\n"
        "instance, the reactant and product molecules to be picked out of a CML\n"
        "React file using CML. Each format has an element which is regarded as\n"
        "defining the object that OpenBabel will convert. For CML this is\n"
        "<molecule>. Files can have multiple objects and these can be treated\n"
        "the same as with other multiple object formats like SMILES and MDL\n"
        "Molfile. So conversion can start at the nth object using the ``-fn`` option\n"
        "and finish before the end using the ``-ln`` option. Multiple object XML files\n"
        "also can be indexed and searched using FastSearch, although this has\n"
        "not yet been extensively tested.\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.xml-cml.org/";}

    virtual const char* GetMIMEType()
    { return "chemical/x-cml"; };

    virtual unsigned int Flags()
    {
      return READXML | ZEROATOMSOK;
    };

    virtual bool WriteChemObject(OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  protected:
    virtual bool DoElement(const string& name);
    virtual bool EndElement(const string& name);
    virtual const char* EndTag(){ return "/molecule>"; };
  private:
    typedef vector< vector< pair<string,string> > > cmlArray;
    bool TransferArray(cmlArray& arr);
    bool TransferElement(cmlArray& arr);
    bool DoAtoms();
    bool DoBonds();
    bool DoHCounts();
    bool DoMolWideData();
    bool ParseFormula(string& formula, OBMol* pmol);
    void ReadNasaThermo();

    void MakeAtomIds(OBMol& mol, vector<string>& atomIDs);
    void WriteFormula(OBMol mol); //passes copy of mol
    void WriteMetadataList(OBMol& mol);
    string getTimestr();
    void WriteBondStereo(OBBond* pbond, vector<string>& atomIDs);
    void WriteCrystal(OBMol& mol);
    void WriteProperties(OBMol& mol, bool& propertyListWritten);
    void WriteThermo(OBMol& mol, bool& propertyListWritten);
    string GetMolID();//for error mesaages
    bool WriteInChI(OBMol& mol);
    bool WriteScalarProperty(OBMol& mol, const char* title, double value,
      const char* dictref=NULL, const char* units=NULL, const char* convention=NULL);

    bool WriteVibrationData(OBMol& mol);
    bool WriteRotationData(OBMol& mol);

  private:
    map<string,int> AtomMap; //key=atom id, value= ob atom index
    cmlArray AtomArray;
    cmlArray BondArray;
    vector<int> HCounts; // for each atom, either -1 or the value of the hydrogenCount
    vector< pair<string,string> > cmlBondOrAtom; //for cml1 only
    vector< pair<string,string> > molWideData;
    bool inBondArray; //for cml1 only
    bool inFormula;
    string RawFormula;
    xmlChar* prefix;
    string CurrentAtomID;
    int CrystalScalarsNeeded, PropertyScalarsNeeded, TransformsNeeded;
    vector<double> CrystalVals;
    OBUnitCell* pUnitCell;
    SpaceGroup _SpaceGroup;
    string SpaceGroupName;
    string titleonproperty;
  };

  ////////////////////////////////////////////////////////////
  //Make an instance of the format class
  CMLFormat theCMLFormat;


  /*
    There are 4 CML styles: CML1, CML2, both with and without array forms.
    All styles are converted into the same internal structure in AtomArray
    and BondArray which contains pairs of (attribute)name/value pairs for
    each atom or bond. At the end of molecule this is analysed in DoAtoms()
    and DoBonds() to construct an OBMol.
  */

  //Callback routines
  ///////////////////////////////////////////////////////
  bool CMLFormat::DoElement(const string& name)
  {
    //A linear search is good enough for <20 element names; commonest at start.
    string value;
    if(name=="atom")
      {
        cmlBondOrAtom.clear();
        int IsEmpty = xmlTextReaderIsEmptyElement(reader());
        TransferElement(AtomArray);
        if(IsEmpty==1) //have to push here because end atom may not be called
          AtomArray.push_back(cmlBondOrAtom);
      }
    else if(name=="bond")
      {
        cmlBondOrAtom.clear();
        int IsEmpty = xmlTextReaderIsEmptyElement(reader());
        TransferElement(BondArray);
        if(IsEmpty==1)
          BondArray.push_back(cmlBondOrAtom);
      }
    else if(name=="molecule" || name=="jobstep") //hack for molpro
      {
        //Ignore atoms with "ref" attributes
        if(xmlTextReaderGetAttribute(reader(), BAD_CAST "ref"))
          return true;
        _pmol->Clear();
        AtomArray.clear();
        BondArray.clear();
        HCounts.clear();
        inBondArray = false;
        inFormula = false;
        RawFormula.erase();
        molWideData.clear();
        CrystalScalarsNeeded=0;
        CrystalVals.clear();
        pUnitCell = NULL;
        PropertyScalarsNeeded=0;

        if(++_embedlevel)
          return true; //ignore if already inside a molecule
        _pmol->BeginModify();
        AtomMap.clear();

        const xmlChar* ptitle  = xmlTextReaderGetAttribute(reader(), BAD_CAST "title");
        if(!ptitle)
          ptitle  = xmlTextReaderGetAttribute(reader(), BAD_CAST "id");
        if(!ptitle)
          ptitle  = xmlTextReaderGetAttribute(reader(), BAD_CAST "molID");//Marvin
        if(ptitle)
          _pmol->SetTitle((const char*)ptitle);

        ptitle = xmlTextReaderGetAttribute(reader(), BAD_CAST "spinMultiplicity");
        if(ptitle)
          _pmol->SetTotalSpinMultiplicity(atoi((const char*)ptitle));

        // free((void*)ptitle);//libxml2 doc says "The string must be deallocated by the caller."

      }
    else if(name=="atomArray")
      {
        if(!inFormula) //do nothing when a child of <formula>
        {
          inBondArray=false;
          TransferArray(AtomArray);
        }
      }
    else if(name=="bondArray")
      {
        inBondArray=true;
        TransferArray(BondArray);
      }
    else if(name=="atomParity" || name=="bondStereo")
      {
        //Save in molWideData:
        //the content,  the atomRefs4 attribute, and (for atomParity only) the centralAtom
        string atrefs4("atomRefs4");
        value = _pxmlConv->GetAttribute(atrefs4.c_str());
        pair<string,string> atomrefdata(atrefs4,value);

        xmlTextReaderRead(reader());
        const xmlChar* pvalue = xmlTextReaderConstValue(reader());
        if(pvalue)
          {
            value = (const char*)pvalue;
            Trim(value);
            pair<string,string> nameAndvalue(name,value);
            molWideData.push_back(nameAndvalue);
            molWideData.push_back(atomrefdata);

            stringstream ss;
            if(name=="atomParity")
              ss << AtomArray.size()+1; //index of current atom
            else
              ss << BondArray.size(); //index of current bond
            pair<string,string> atdata("centralAtomOrBond",ss.str());
            molWideData.push_back(atdata);
          }
      }
    else if(name=="name")
      {
        if(_pmol)
          _pmol->SetTitle(_pxmlConv->GetContent().c_str());
      }
    else if(name=="formula")
      {
        if(!xmlTextReaderIsEmptyElement(reader()))
          inFormula=true;
        //Only concise form is currently supported
        const xmlChar* pformula = xmlTextReaderGetAttribute(reader(), BAD_CAST "concise");
        if(pformula)
        {
          RawFormula = (const char*)pformula;
          // free((void*)pformula);
        }
      }
    else if(name=="crystal")
      {
        CrystalScalarsNeeded = 6;
      }
    else if(name=="scalar")
      {
        if(CrystalScalarsNeeded)
          {
            xmlTextReaderRead(reader());
            const xmlChar* pvalue = xmlTextReaderConstValue(reader());
            if(pvalue)
              {
                CrystalVals.push_back(atof((const char*)pvalue));
                if(--CrystalScalarsNeeded==0)
                  {
                    pUnitCell = new OBUnitCell;
                    pUnitCell->SetOrigin(fileformatInput);
                    pUnitCell->SetData(CrystalVals[0],CrystalVals[1],CrystalVals[2],
                                       CrystalVals[3],CrystalVals[4],CrystalVals[5]);
                    _pmol->SetData(pUnitCell);
                  }
              }
          }
        else if(PropertyScalarsNeeded)
          {
            //Reads OBPairData(like SDF properties). Name is in scalar title or id attribute
            const xmlChar* pattr  = xmlTextReaderGetAttribute(reader(), BAD_CAST "title");
            if(!pattr)
              pattr  = xmlTextReaderGetAttribute(reader(), BAD_CAST "id");

            string attr;
            if(pattr)
              attr = (const char*)pattr;
            else
              attr = titleonproperty;
            // free((void*)pattr);//"The string must be deallocated by the caller."

            xmlTextReaderRead(reader());
            const xmlChar* pvalue = xmlTextReaderConstValue(reader());

            if(titleonproperty.find("ZPE")!=string::npos)
            {
              double energy;
              stringstream ss((const char*)pvalue);
              ss >> energy; //units are kJ/mol
              const double CALSTOJOULES = 4.1816;
              _pmol->SetEnergy(energy/CALSTOJOULES);
            }
            else
            {
              if(pvalue && !attr.empty())
              {
                OBPairData *dp = new OBPairData;
                dp->SetAttribute(attr);
                string val((const char*)pvalue);
                dp->SetValue(Trim(val));
                dp->SetOrigin(fileformatInput);
                _pmol->SetData(dp);
              }
            }
            PropertyScalarsNeeded=0;
          }
      }
    else if(name=="array" && PropertyScalarsNeeded)
      {
        //Read vibrational frequencies and rotational constants from properties
        xmlTextReaderRead(reader());
        const xmlChar* pvalue = xmlTextReaderConstValue(reader());
        string value;
        if(pvalue)
          value = (const char*)pvalue;
        vector<string> items;
        tokenize(items,value);

        if(titleonproperty.find("vibFreqs")!=string::npos)
        {
          vector< vector< vector3 > > vLx;
          vector<double> vFrequencies, vIntensities;
          for(unsigned i=0;i<items.size();++i)
            vFrequencies.push_back(atof(items[i].c_str()));

          OBVibrationData* vd = new OBVibrationData;
          vd->SetData(vLx, vFrequencies, vIntensities);
          vd->SetOrigin(fileformatInput);
          _pmol->SetData(vd);
        }

        else if(titleonproperty.find("rotConsts")!=string::npos)
        {
          const double WAVENUM_TO_GHZ=30.0;
          vector<double> rotConsts;
          for(unsigned i=0;i<items.size();++i)
            rotConsts.push_back(atof(items[i].c_str()) * WAVENUM_TO_GHZ);

          OBRotationData* rd = new OBRotationData;
          rd->SetData(OBRotationData::UNKNOWN, rotConsts, 1);//rotor type and symmetry number unknown
          rd->SetOrigin(fileformatInput);
          _pmol->SetData(rd);
        }

        PropertyScalarsNeeded = 0;
      }
    else if(name=="symmetry")
      {
        const xmlChar* pname  = xmlTextReaderGetAttribute(reader(), BAD_CAST "spaceGroup");
        if (pname)
		  {
            SpaceGroupName = (const char*)pname;
            // free((void*)pname);
          }
      }
     else if(name=="transform3")
      {
        xmlTextReaderRead(reader());
        const xmlChar* ptransform = xmlTextReaderConstValue(reader());
        if (ptransform)
		  {
            string t = (const char*)ptransform;
            _SpaceGroup.AddTransform(t);
            // free((void*)ptransform);
         }
      }
   else if(name=="property")
      {
        //***pattr need to be deleted***
        const char* pattr  = (const char*)xmlTextReaderGetAttribute(reader(), BAD_CAST "dictRef");
        if(pattr && !strcmp(pattr,"Thermo_OldNasa"))
          ReadNasaThermo();
        else
          {
            if(!pattr) // no dictRef; look for title on scalar
              pattr  = (const char*)xmlTextReaderGetAttribute(reader(), BAD_CAST "title");
            if(pattr)
              titleonproperty = pattr;
            else
              titleonproperty.clear();
            PropertyScalarsNeeded = 1;
          }
      }

    // CML1 elements
    else	if(name=="string" || name=="float" || name=="integer"
             || name=="coordinate3"|| name=="coordinate2")
      {
        string name = _pxmlConv->GetAttribute("builtin");
        xmlTextReaderRead(reader());
        const xmlChar* pvalue = xmlTextReaderConstValue(reader());
        if(!pvalue)
          return false;
        string value = (const char*)pvalue;
        Trim(value);
        pair<string,string> nameAndvalue(name,value);
        cmlBondOrAtom.push_back(nameAndvalue);
      }
    else	if(name=="stringArray" || name=="floatArray" || name=="integerArray")
      {
        string name = _pxmlConv->GetAttribute("builtin");
        //		cmlArray& arr = (name=="atomRef1" || name=="atomRef2" || name=="order")
        //			? BondArray : AtomArray;
        cmlArray& arr = inBondArray ? BondArray : AtomArray;

        xmlTextReaderRead(reader());
        const xmlChar* pvalue = xmlTextReaderConstValue(reader());
        if(!pvalue)
          return false;
        string value = (const char*)pvalue;

        vector<string> items;
        tokenize(items,value);
        if(arr.size()<items.size())
          arr.resize(items.size());
        unsigned int i;
        for(i=0;i<items.size();++i)
          {
            pair<string,string> nameAndvalue(name,items[i]);
            arr[i].push_back(nameAndvalue);
          }
      }

    //The end element event would not be called for <element/>, so call it explicitly.
    if(xmlTextReaderIsEmptyElement(reader())==1)
      return EndElement(name);

    return true;
  }

  //////////////////////////////////////////////////////
  bool CMLFormat::EndElement(const string& name)
  {
    if(name=="atom")
      {
        //ok for cml1 but is not called at end of <atom.../>
        AtomArray.push_back(cmlBondOrAtom);
      }
    else if(name=="bond")
      {
        BondArray.push_back(cmlBondOrAtom);
      }
    else if(name=="formula")
      inFormula=false;
    else if(name=="molecule" || name=="jobstep") //hack for molpro
      {
        if(!DoAtoms() || !DoBonds() || !DoHCounts() || !DoMolWideData())
          return false;

        if (_pmol->GetDimension()==0)
          StereoFrom0D(_pmol); // Remove any spurious stereos (due to symmetry)

        //Use formula only if nothing else provided
        if(_pmol->NumAtoms()==0 && !RawFormula.empty())
          if(!ParseFormula(RawFormula, _pmol))
            obErrorLog.ThrowError(_pmol->GetTitle(),"Error in formula", obError);

        _pmol->AssignSpinMultiplicity();
        _pmol->EndModify();
        return (--_embedlevel>=0); //false to stop parsing if no further embedded mols
        //		return false;//means stop parsing
      }
     else if(name=="symmetry")
      {
        if(!SpaceGroupName.empty())
        {
          const SpaceGroup *group = SpaceGroup::GetSpaceGroup(SpaceGroupName);
          if ((!group || !(_SpaceGroup == *group)) && _SpaceGroup.IsValid())
            group = SpaceGroup::Find(&_SpaceGroup);
          if (group)
            pUnitCell->SetSpaceGroup(group);
          else
            pUnitCell->SetSpaceGroup(SpaceGroupName);
        }
      }
    return true;
  }

  /////////////////////////////////////////////////////////

  static unsigned int GetAtomicNumAndIsotope(const char* symbol, int *isotope)
  {
    const char* p = symbol;
    switch (p[0]) {
    case 'D':
      if (p[1] == '\0') {
        *isotope = 2;
        return 1;
      }
      break;
    case 'T':
      if (p[1] == '\0') {
        *isotope = 3;
        return 1;
      }
      break;
    }
    return OBElements::GetAtomicNum(symbol);
  }

  static const char* FindStartOfAtomClass(const char* atomid)
  {
    // Try to find a match to 'a' followed by a number followed by _ followed by at least one digit
    if (atomid[0] != 'a')
      return (const char*)0; // Needs to start with 'a'
    const char *p = atomid + 1;
    while(*p >= '0' && *p <= '9')
      p++;
    if (p == atomid + 1)
      return (const char*)0; // No digits
    if (*p != '_')
      return (const char*)0;
    p++;
    if (*p >= '0' && *p <= '9')
      return p;
    return (const char*)0;
  }

  ///Interprets atoms from AtomArray and writes then to an OBMol
  bool CMLFormat::DoAtoms()
  {
    int dim=0; //dimension of molecule
    bool use2d = _pxmlConv->IsOption("2", OBConversion::INOPTIONS);

    int nAtoms=_pmol->NumAtoms();//was 0
    cmlArray::iterator AtomIter;
    for(AtomIter=AtomArray.begin();AtomIter!=AtomArray.end();++AtomIter)
      {
        //		OBAtom obatom;
        OBAtom* pAtom = _pmol->NewAtom();
        nAtoms++;
        int nhvy = nAtoms;
        int hcount = -1; // default value which may be overridden by hydrogenCount below

        double x=0,y=0,z=0;
        bool using3=false, using2=false, usingFract=false;

        vector<pair<string,string> >::iterator AttributeIter;
        for(AttributeIter=AtomIter->begin();AttributeIter!=AtomIter->end();++AttributeIter)
          {
            string& attrname = AttributeIter->first;
            string& value    = AttributeIter->second;

            if(attrname=="id" || attrname=="atomId" || attrname=="atomID")//which one correct?
              {
                Trim(value);
                if(AtomMap.count(value)>0)
                  obErrorLog.ThrowError(GetMolID(),"The atom id " + value + " is not unique", obWarning);
                AtomMap[value] = nhvy;//nAtoms;

                //If the id ends with "_NUMBER", then NUMBER is taken as an atom class
                const char* atomclass = FindStartOfAtomClass(value.c_str());
                if (atomclass) {
                  OBPairInteger *pi = new OBPairInteger();
                  pi->SetAttribute("Atom Class");
                  pi->SetValue(atoi(atomclass));
                  pi->SetOrigin(fileformatInput);
                  pAtom->SetData(pi);
                }
                continue;
              }
            else if(attrname=="elementType")
              {
                int atno, iso=0;
                atno = GetAtomicNumAndIsotope(value.c_str(), &iso);
                pAtom->SetAtomicNum(atno);
                if(iso)
                  pAtom->SetIsotope(iso);
                continue;
              }

            //If more than one set of coordinates provided,
            //prefer 3D over 2D over 3Dfractional,
            //but if use2d is true, prefer 2D over 3D
            else if((attrname=="x3" || attrname=="y3" || attrname=="z3" || attrname=="xyz3") && !use2d)
            {
              using3 = true;
              usingFract = false;
            }
            else if((attrname=="x2" || attrname=="y2" || attrname=="z2" || attrname=="xy2") && !using3)
            {
              using2 = true;
              usingFract = false;
            }
            else if(pUnitCell && !using3 && !using2
              && (attrname=="xFract" || attrname=="yFract" || attrname=="zFract"))
              usingFract=true;

            if ((using3     && attrname=="x3") ||
                (using2     && attrname=="x2") ||
                (usingFract && attrname=="xFract")) {
              x=strtod(value.c_str(),NULL);
            }
            else if ((using3     && attrname=="y3") ||
                     (using2     && attrname=="y2") ||
                     (usingFract && attrname=="yFract")) {
              y=strtod(value.c_str(),NULL);
            }
            else if ((using3     && attrname=="z3") ||
                     (using2     && attrname=="z2") ||
                     (usingFract && attrname=="zFract")) {
              z=strtod(value.c_str(),NULL);
            }
            else if(using2 && attrname=="xy2")
              {
                vector<string> vals;
                tokenize(vals,value);
                if(vals.size()==2)
                  {
                    x=strtod(vals[0].c_str(),NULL);
                    y=strtod(vals[1].c_str(),NULL);
                  }
              }
            else if(using3 && attrname=="xyz3")
              {
                vector<string> vals;
                tokenize(vals,value);
                if(vals.size()==3)
                  {
                    x=strtod(vals[0].c_str(),NULL);
                    y=strtod(vals[1].c_str(),NULL);
                    z=strtod(vals[2].c_str(),NULL);
                  }
              }

            if(attrname=="hydrogenCount")
              {
                //Actually adding H atoms to the structure is deferred until the explicit
                //structure is complete, because hydrogenCount may include explicit H, bug#3014855
                hcount = atoi(value.c_str());

               /* int nhvy = nAtoms;
                  int i;
                  for(i=0;i<atoi(value.c_str());++i)
                    {
                      OBAtom* hatom = _pmol->NewAtom();
                      hatom->SetAtomicNum(1);
                      hatom->SetType("H");
                      _pmol->AddBond(nhvy,_pmol->NumAtoms(),1);
                      ++nAtoms;
                    }
               */
              }

            else if(attrname=="formalCharge")
              pAtom->SetFormalCharge(atoi(value.c_str()));

            else if(attrname=="label")
              {
                OBPairData *label = new OBPairData();
                label->SetAttribute("label");
                label->SetValue(value.c_str());
                pAtom->SetData(label);
              }

            else if(attrname=="color")
              {
                OBPairData *color = new OBPairData();
                color->SetAttribute("color");
                color->SetValue(value.c_str());
                pAtom->SetData(color);
              }

            else if(attrname=="radius")
              {
                OBPairData *radius = new OBPairData();
                radius->SetAttribute("radius");
                radius->SetValue(value.c_str());
                pAtom->SetData(radius);
              }

            else if(attrname=="spinMultiplicity")
              pAtom->SetSpinMultiplicity(atoi(value.c_str()));

            /*else if(attrname=="atomRefs4")//from atomParity element (but there is no such thing!)
              {
                vector<string> ids;
                tokenize(ids, value);

                const xmlChar* pvalue = xmlTextReaderConstValue(reader());
                int parity = 0;
                if (pvalue)
                  parity = atoi((const char*)pvalue);
                if (parity != 0) { // Should be +1 or -1
                  TetSym ts;
                  ts.atomrefs = ids;
                  ts.parity = parity;
                  tetsyms.push_back(ts);
                }
              }*/

            else if(attrname=="radical") //Marvin extension
              {
                int spin=0;
                if(value=="monovalent")
                  spin=2;
                else if(value=="divalent")
                  spin=3;
                else if(value=="divalent3")
                  spin=3;
                else if(value=="divalent1")
                  spin=1;
                pAtom->SetSpinMultiplicity(spin);
              }
            else if(attrname=="isotopeNumber" || attrname=="isotope")
              pAtom->SetIsotope(atoi(value.c_str()));

          } //each attribute

          //Save hydrogen count
          HCounts.push_back(hcount);

          //Save atom coordinates
          if(using3 || usingFract)
            dim=3;
          else if(using2)
          {
            dim=2;
            z=0.0;
          }
          else
            dim=0;
          if(usingFract)
            {
              //Coordinates are fractional
              vector3 v;
              v.Set(x, y, z);
              v = pUnitCell->FractionalToCartesian(v);
              pAtom->SetVector(v);
            }
          else
            pAtom->SetVector(x, y, z);
      }//each atom

    _pmol->SetDimension(dim);
    return true;
  }
  /////////////////////////////////////////////////////////////////////

  ///Interprets bonds from BondArray and writes then to an OBMol
  bool CMLFormat::DoBonds()
  {
    vector<pair<string,string> >::iterator AttributeIter;
    cmlArray::iterator BondIter;
    bool HaveWarned = false;
    for(BondIter=BondArray.begin();BondIter!=BondArray.end();++BondIter)
      {
        int indx1=0,indx2=0, ord=0;
        string bondstereo, BondStereoRefs;
        string colour;
        string label;
        bool PossibleBond = false;

        for(AttributeIter=BondIter->begin();AttributeIter!=BondIter->end();++AttributeIter)
          {
            string attrname = AttributeIter->first;
            string value    = AttributeIter->second;
            Trim(value);


            if(attrname.compare(0, 7, "atomRef")==0) //generic
              {
                PossibleBond = true;
                string::size_type pos = value.find(' ');

                if(!HaveWarned && (attrname=="atomRefs1"
                                   || (attrname=="atomRefs2" && pos==string::npos)))
                  {
                    obErrorLog.ThrowError(GetMolID(),
                                          attrname + " is not legal CML in this context, "
                                          "but OpenBabel will attempt to understand what was meant.", obWarning);
                    HaveWarned = true;
                  }

                if(indx1==0)
                  {
                    if(pos!=string::npos)
                      {
                        indx1 = AtomMap[value.substr(0,pos)];
                        string temp =value.substr(pos+1);
                        indx2 = AtomMap[Trim(temp)];
//C4239                        indx2 = AtomMap[Trim(value.substr(pos+1))];

                      }
                    else
                      indx1 = AtomMap[value];
                  }
                else
                  {
                    if(indx2==0)
                      indx2 = AtomMap[value];
                    else
                      indx1=-1; //forces error
                  }
              }
            else if(attrname=="order")
              {
                const char bo = value[0];
                if(bo=='S')
                  ord=1;
                else if(bo=='D')
                  ord=2;
                else if(bo=='T')
                  ord=3;
                else if(bo=='A')
                  ord=5;
                else {
                  char* endptr;
                  ord = strtol(value.c_str(), &endptr, 10);
                }
              }

            else if(attrname=="color")
              colour=value[0];

            else if(attrname=="label")
              label = value;
          }

        if(PossibleBond)
          {
            if(indx1<=0 || indx2<=0)
              {
                obErrorLog.ThrowError(GetMolID(),"Incorrect bond attributes", obError);
                return false;
              }
            if(ord==0) //Bonds are single if order is not specified
              {
                ord=1;
                //But unspecied bond order means cannot assign spinmultiplicity
                _pmol->SetIsPatternStructure();
              }
            _pmol->AddBond(indx1,indx2,ord,0);

            if(!colour.empty())
              {
                OBPairData *dp = new OBPairData();
                dp->SetAttribute("color");
                dp->SetValue(colour.c_str());
                _pmol->GetBond(_pmol->NumBonds()-1)->SetData(dp);
              }
            if(!label.empty())
              {
                OBPairData *dp = new OBPairData();
                dp->SetAttribute("label");
                dp->SetValue(label.c_str());
                _pmol->GetBond(_pmol->NumBonds()-1)->SetData(dp);
              }
          }
      }

    return true;
  }

  /////////////////////////////////////////////////////////////////
  bool CMLFormat::DoHCounts()
  {
    //Add extra H atoms so that each atom has the value of its attribute hydrogenCount
    FOR_ATOMS_OF_MOL(atom, _pmol)
    {
      int hcount = HCounts[atom->GetIdx() - 1];
      if (hcount == -1)
      {
        OBAtomAssignTypicalImplicitHydrogens(&*atom);
        continue;
      }

      int explH = atom->ExplicitHydrogenCount(); // includes H isotopes
      if(explH > hcount)
      {
        map<string,int>::iterator it;
        for(it=AtomMap.begin();it!=AtomMap.end();++it)
          if(it->second == atom->GetIdx())
            break;
        stringstream ss;
        ss << "In atom " << it->first << " the number of explicit hydrogens exceeds the hydrogenCount attribute.";
        obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
        return false;
      }
      atom->SetImplicitHCount(hcount - explH);
    }
    return true;
  }

  bool CMLFormat::DoMolWideData()
  {
    //Handle atomParity and bondStereo
    vector<pair<string,string> >::iterator AttributeIter;
    for(AttributeIter=molWideData.begin();AttributeIter!=molWideData.end();++AttributeIter)
      {
        string name  = AttributeIter->first;
        string value = AttributeIter->second;

        if(name=="atomParity" || name=="bondStereo")
          {
            vector<unsigned int> AtomRefIdx;

            string nextname = (++AttributeIter)->first;
            string atrefsvalue = AttributeIter->second;
            if(nextname=="atomRefs4" && !atrefsvalue.empty())
              {
                vector<string> ids;
                tokenize(ids, atrefsvalue);
                int i;
                for(i=0;i<4;++i)
                  AtomRefIdx.push_back(AtomMap[ids[i]]);
              }

            nextname = (++AttributeIter)->first;
            if(!(nextname=="centralAtomOrBond"))
              return false;

            int Idx = atoi(AttributeIter->second.c_str());
            if(name=="atomParity")
              {
                OBAtom* patom = _pmol->GetAtom(Idx);
                if(!patom)
                  return false;

                OBStereo::Ref center = patom->GetId();
                OBStereo::Ref from = _pmol->GetAtom(AtomRefIdx[0])->GetId();
                if (from == center)
                  from = OBStereo::ImplicitRef;

                OBStereo::Refs refs;
                vector<unsigned int>::const_iterator idx_cit=AtomRefIdx.begin();
                ++idx_cit;
                for (; idx_cit!=AtomRefIdx.end(); ++idx_cit) {
                  OBStereo::Ref id = _pmol->GetAtom(*idx_cit)->GetId();
                  if (id == center)
                    id = OBStereo::ImplicitRef;
                  refs.push_back(id);
                }

                int parity = atoi(value.c_str());
                OBStereo::Winding winding = OBStereo::Clockwise; // parity > 0
                if (parity < 0)
                  winding = OBStereo::AntiClockwise;
                else if (parity == 0) // What to do with parity of 0?
                  return false;

                OBTetrahedralStereo::Config cfg = OBTetrahedralStereo::Config(
                                                        center, from, refs, winding, OBStereo::ViewFrom);
                OBTetrahedralStereo *th = new OBTetrahedralStereo(_pmol);
                th->SetConfig(cfg);
                _pmol->SetData(th);
              }
            else //bondStereo
              {
                OBBond* pbond1=NULL;
                OBBond* pbond2=NULL;
                if(atrefsvalue.empty()) // ToDo
                  {
                    OBBond* pDBond = _pmol->GetBond(Idx);
                    //With no atomRefs4, the specification is either W, H,
                    if(value=="W")
                      {
                        pDBond->SetWedge();
                      }
                    else if(value=="H")
                      {
                        pDBond->SetHash();
                      }
                    // ... or ordinary cis/trans
                    if(value!="C" && value!="T")
                      continue;
                    //which is valid only with one substituent on each C

                    OBAtom* pAt1 = pDBond->GetBeginAtom();
                    OBAtom* pAt2 = pDBond->GetEndAtom();
                    FOR_NBORS_OF_ATOM(a1,pAt1)
                      {
                        if (a1->GetAtomicNum() != OBElements::Hydrogen && &*a1 != pAt2)
                          break;
                        pbond1 = _pmol->GetBond(pAt1->GetIdx(),a1->GetIdx());
                      }

                    FOR_NBORS_OF_ATOM(a2,pAt2)
                      {
                        if (a2->GetAtomicNum() != OBElements::Hydrogen && &*a2 != pAt1)
                          break;
                        pbond2 = _pmol->GetBond(pAt2->GetIdx(),a2->GetIdx());
                      }
                  }
                else
                  {
                    pbond1 = _pmol->GetBond(AtomRefIdx[0],AtomRefIdx[1]);
                    pbond2 = _pmol->GetBond(AtomRefIdx[2],AtomRefIdx[3]);
                  }

                if(!pbond1 || !pbond2)
                  continue;

                // Create the list of 4 atomrefs
                OBStereo::Ref begin, end;
                begin = _pmol->GetAtom(AtomRefIdx[1])->GetId();
                end =   _pmol->GetAtom(AtomRefIdx[2])->GetId();
                OBStereo::Refs refs(4);
                refs[0] = _pmol->GetAtom(AtomRefIdx[0])->GetId();
                refs[1] = OBStereo::ImplicitRef;
                FOR_NBORS_OF_ATOM(nbr, _pmol->GetAtomById(begin))
                  if (nbr->GetId()!=end && nbr->GetId()!=refs[0]) {
                    refs[1] = nbr->GetId();
                    break;
                  }
                OBStereo::Ref tmpref = _pmol->GetAtom(AtomRefIdx[3])->GetId();
                OBStereo::Ref finalref = OBStereo::ImplicitRef;
                FOR_NBORS_OF_ATOM(nbr, _pmol->GetAtomById(end))
                  if (nbr->GetId()!=begin && nbr->GetId()!=tmpref) {
                    finalref = nbr->GetId();
                    break;
                  }
                if (value=="C") { // for Cis (tmpref and refs[0] on same side)
                  refs[2] = finalref; refs[3] = tmpref;
                }
                else { // for Trans
                  refs[2] = tmpref; refs[3] = finalref;
                }

                // Create the new stereo object
                OBCisTransStereo::Config ct_cfg = OBCisTransStereo::Config(
                                                       begin, end, refs, OBStereo::ShapeU);
                OBCisTransStereo *ct = new OBCisTransStereo(_pmol);
                ct->SetConfig(ct_cfg);
                _pmol->SetData(ct);
              }
          }
      }
    //Clear here to aid embedded molecules
    AtomArray.clear();
    BondArray.clear();
    molWideData.clear();

    return true;
  }

  //////////////////////////////////////////////////////////
  bool CMLFormat::TransferArray(cmlArray& arr)
  {
    //Reads attributes of the current node, e.g. atomID="a1 a2 a3"
    //parses each of them into their separate items, e.g. a1, a2, a3
    //and pushes them as a pairs in each of the members of the array
    // e.g. ("atomID", "a1") in AtomArray[0], ("atomID", "a2") in AtomArray[1]

    if(xmlTextReaderHasAttributes(reader()))
      {
        int ret = xmlTextReaderMoveToFirstAttribute(reader());
        while(ret==1)
          {
            const xmlChar* pname = xmlTextReaderConstName(reader());
            string name((const char*)pname);
            const xmlChar* pvalue = xmlTextReaderConstValue(reader());
            string value;
            if(pvalue)
              value = (const char*)pvalue;
            vector<string> items;
            tokenize(items,value);
            if(arr.size()<items.size())
              arr.resize(items.size());
            unsigned int i;
            for(i=0;i<items.size();++i)
              {
                pair<string,string> nameAndvalue(name,items[i]);
                arr[i].push_back(nameAndvalue);
              }
            ret = xmlTextReaderMoveToNextAttribute(reader());
          }
      }
    return true;
  }

  bool CMLFormat::TransferElement(cmlArray& arr)
  {
    //Reads the attributes of the current node, e.g. <atom id="a1" elementType="C"/>
    //pushes each of them as a pairs into each of the members of the array
    // e.g. ("id", "a1") and (elementType", "C") will be put into AtomArray[n]
    //where n is the number of times this routine has been called before.

    if(xmlTextReaderHasAttributes(reader()))
      {
        int ret = xmlTextReaderMoveToFirstAttribute(reader());
        while(ret==1)
          {
            const xmlChar* pname = xmlTextReaderConstName(reader());
            string name((const char*)pname);
            const xmlChar* pvalue = xmlTextReaderConstValue(reader());
            string value;
            if(pvalue)
              {
                value = (const char*)pvalue;
                Trim(value);
              }
            pair<string,string> nameAndvalue(name,value);
            cmlBondOrAtom.push_back(nameAndvalue);
            ret = xmlTextReaderMoveToNextAttribute(reader());
          }
      }
    return true;
  }

  bool CMLFormat::ParseFormula(string& formula, OBMol* pmol)
  {
    vector<string> items;
    tokenize(items, formula);
    vector<string>::iterator iSymbol, iNumber;
    for(iSymbol=items.begin();iSymbol!=items.end();++iSymbol)
      {
        iNumber = iSymbol+1;
        if(iNumber==items.end())
          return false;
        int n=atoi(iNumber->c_str());
        int atno, iso=0;
        atno = GetAtomicNumAndIsotope(iSymbol++->c_str(), &iso);
        if(atno<=0 || n<=0)
          return false;
        int i;
        for(i=0;i<n;++i)
          {
            OBAtom* pAtom = pmol->NewAtom();
            pAtom->SetAtomicNum(atno);
            if(iso)
              pAtom->SetIsotope(iso);
          }
      }
    return true;
  }

  void CMLFormat::ReadNasaThermo()
  {
    //Do all NasaThermo data here
    OBNasaThermoData* pTD = new OBNasaThermoData;
    pTD->SetOrigin(fileformatInput);
    _pmol->SetData(pTD);
    for(;;)
      {
        xmlTextReaderRead(reader());
        int typ = xmlTextReaderNodeType(reader());
        if(typ==XML_READER_TYPE_SIGNIFICANT_WHITESPACE)
          continue;
        const char* pname = (const char*)xmlTextReaderConstLocalName(reader());
        if(typ==XML_READER_TYPE_END_ELEMENT)
          {
            if(!strcmp(pname,"property"))//end of element
              return;
            else
              continue;
          }
        const char * pattr  = (const char*)xmlTextReaderGetAttribute(reader(), BAD_CAST "dictRef");
        xmlTextReaderRead(reader());
        const char* pvalue = (const char*)xmlTextReaderConstValue(reader());
        if(pattr && pvalue)
          {
            if(!strcmp(pattr,"NasaLowT"))
              pTD->SetLoT(atof(pvalue));
            else if(!strcmp(pattr,"NasaHighT"))
              pTD->SetHiT(atof(pvalue));
            else if(!strcmp(pattr,"NasaMidT"))
              pTD->SetMidT(atof(pvalue));
            else if(!strcmp(pattr,"NasaCoeffs"))
              {
                vector<string> vals;
                tokenize(vals, pvalue);
                for(int i=0;i<14;++i)
                  pTD->SetCoeff(i, atof(vals[i].c_str()));
              }
          }
      }
  }


  void CMLFormat::WriteMetadataList(OBMol& mol)
  {
    static const xmlChar C_METADATALIST[] = "metadataList";
    static const xmlChar C_METADATA[]     = "metadata";
    static const xmlChar C_TITLE[]        = "title";
    static const xmlChar C_NAME[]         = "name";
    static const xmlChar C_CONTENT[]      = "content";

    xmlTextWriterStartElement(writer(), C_METADATALIST);

    if(mol.HasData(OBGenericDataType::CommentData))
    {
      OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
      xmlTextWriterStartElement(writer(), C_METADATA);
      xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:description");
      xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST cd->GetData().c_str());
      xmlTextWriterEndElement(writer());
    }

    xmlTextWriterStartElement(writer(), C_METADATA);
    xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:source");
    xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "unknown");
    xmlTextWriterEndElement(writer());

    xmlTextWriterStartElement(writer(), C_METADATA);
    xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:creator");
    string version("OpenBabel version ");
    version += BABEL_VERSION;
    xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST version.c_str());
    xmlTextWriterEndElement(writer());

    xmlTextWriterStartElement(writer(), C_METADATA);
    xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:contributor");
    xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "unknown");
    xmlTextWriterEndElement(writer());

    xmlTextWriterStartElement(writer(), C_METADATA);
    xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:date");
    xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST getTimestr().c_str());
    xmlTextWriterEndElement(writer());

/*    xmlTextWriterStartElement(writer(), C_METADATA);
    xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "cmlm:structure");
    xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "yes");
    xmlTextWriterEndElement(writer());
*/
    xmlTextWriterEndElement(writer());
  }

  string CMLFormat::getTimestr()
  {
    const int TIME_STR_SIZE = 64;
    time_t akttime;                              /* Systemtime                        */
    char timestr[TIME_STR_SIZE + 1] = "";        /* Timestring                        */
    size_t time_res;                             /* Result of strftime                */

    /* ---- Get the system-time ---- */
    akttime = time((time_t *) NULL);
    time_res = strftime(timestr,
                        TIME_STR_SIZE,
                        "%a %b %d %H:%M:%S %Z %Y",
                        localtime((time_t *) &akttime)
                        );
    return string(timestr);
  }

  /////////////////////////////////////////////////////////////

  bool CMLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    static const xmlChar C_MOLECULE[]   = "molecule";
    static const xmlChar C_CML[]        = "cml";
    static const xmlChar C_ATOMARRAY[]  = "atomArray";
    static const xmlChar C_BONDARRAY[]  = "bondArray";
    static const xmlChar C_ATOM[]       = "atom";
    static const xmlChar C_BOND[]       = "bond";
    static const xmlChar C_ID[]         = "id";
    // static const xmlChar C_TITLE[]      = "title";
    static const xmlChar C_NAME[]       = "name";
    static const xmlChar C_ATOMPARITY[] = "atomParity";
    static const xmlChar C_BONDSTEREO[] = "bondStereo";

    static const xmlChar C_X2[]               = "x2";
    static const xmlChar C_Y2[]               = "y2";
    static const xmlChar C_X3[]               = "x3";
    static const xmlChar C_Y3[]               = "y3";
    static const xmlChar C_Z3[]               = "z3";
    static const xmlChar C_XFRACT[]               = "xFract";
    static const xmlChar C_YFRACT[]               = "yFract";
    static const xmlChar C_ZFRACT[]               = "zFract";
    static const xmlChar C_ATOMID[]           = "atomID";
    static const xmlChar C_ELEMENTTYPE[]      = "elementType";
    static const xmlChar C_ISOTOPE[]          = "isotope";
    static const xmlChar C_SPINMULTIPLICITY[] = "spinMultiplicity";
    static const xmlChar C_HYDROGENCOUNT[]    = "hydrogenCount";
    static const xmlChar C_FORMALCHARGE[]     = "formalCharge";
    static const xmlChar C_ATOMREFS2[]        = "atomRefs2";
    static const xmlChar C_ATOMREF1[]         = "atomRef1";
    static const xmlChar C_ATOMREF2[]         = "atomRef2";
    static const xmlChar C_ORDER[]            = "order";
    static const xmlChar C_ATOMREFS4[]        = "atomRefs4";
    static const xmlChar C_DESCRIPTION[]      = "description";
    /* defined in other functions
       static const xmlChar C_FORMULA[] = "formula";
       static const xmlChar C_CONCISE[] = "concise";
       static const xmlChar C_PROPERTYLIST[] = "propertyList";
       static const xmlChar C_PROPERTY[] = "property";
       static const xmlChar C_SCALAR[] = "scalar";
    */
    static const xmlChar C_LABEL[] = "label";
    static const xmlChar C_COLOR[] = "color";
    static const xmlChar C_RADIUS[] = "radius";
    //CML1
    static const xmlChar C_STRING[]       = "string";
    static const xmlChar C_INTEGER[]      = "integer";
    static const xmlChar C_FLOAT[]        = "float";
    static const xmlChar C_BUILTIN[]      = "builtin";
    static const xmlChar C_STRINGARRAY[]  = "stringArray";
    static const xmlChar C_INTEGERARRAY[] = "integerArray";
    static const xmlChar C_FLOATARRAY[]   = "floatArray";
    /* used as ordinary text
       atomRef
    */

    const xmlChar* C_X3orFRACT = C_X3; //Non-fraction coordinates are the default
    const xmlChar* C_Y3orFRACT = C_Y3;
    const xmlChar* C_Z3orFRACT = C_Z3;

    _pxmlConv = XMLConversion::GetDerived(pConv,false);
    if(!_pxmlConv)
      return false;

    bool cml1 = _pxmlConv->IsOption("1");
    bool arrayform = _pxmlConv->IsOption("a");
    bool WriteAromaticBonds =  _pxmlConv->IsOption("A");
    prefix = BAD_CAST _pxmlConv->IsOption("N");
    xmlChar* uri=NULL;

    //Write the header on the first object (incl OBReaction)
    //unless x option set or if has been called from elsewhere (e.g. CMLReact)
    if(!_pxmlConv->IsOption("MolsNotStandalone") && _pxmlConv->GetOutputIndex()==1)
      {
        if(!_pxmlConv->IsOption("x"))
          {
            xmlTextWriterStartDocument(writer(), NULL, NULL, NULL);
            if(cml1)
              uri = BAD_CAST CML1NamespaceURI();
            else
              uri=BAD_CAST NamespaceURI();// not the old CML2NamespaceURI();
          }
        //If more than one molecule to be output, write <cml> at start and </cml> at end.
        if(!_pxmlConv->IsLast())
          {
            xmlTextWriterStartElementNS(writer(), prefix, C_CML, uri);
            uri=NULL;
          }
      }

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
    {
#ifdef HAVE_SHARED_POINTER
        OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
        if(!pReact)
          return false;
        //Use CMLReact to convert OBReaction object
        OBFormat* pCMLRFormat = pConv->FindFormat("cmlr");
        if(!pCMLRFormat)
        {
          obErrorLog.ThrowError(__FUNCTION__, "Cannot find CMLReact format", obError);
          return false;
        }
        //Disable list option and supress topping and tailing in CMLReactFormat.
        _pxmlConv->AddOption("l", OBConversion::OUTOPTIONS);
        _pxmlConv->AddOption("ReactionsNotStandalone", OBConversion::OUTOPTIONS);
        bool ret = pCMLRFormat->WriteMolecule(pOb,_pxmlConv);
        _pxmlConv->RemoveOption("ReactionsNotStandalone", OBConversion::OUTOPTIONS);
        return ret;
#else
        return false;
#endif
    }


    OBMol &mol = *pmol;

    int numbonds = mol.NumBonds();

    bool UseFormulaWithNoBonds=false; //before 2.3.1 was true;

    int dim = mol.GetDimension();


    xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, uri);

    const char* id = mol.GetTitle();
    if(*id)
      {
        string name(id);
        //If name is a filename with a path, remove path and extension
        string::size_type pos;
        pos = name.find_last_of("/\\:");
        if(pos!=string::npos)
        {
          name.erase(0, pos+1);
          pos = name.rfind('.');
          if(pos!=string::npos)
            name.erase(pos);
        }

        if(!isalpha(name[0])) //since ids have to start with a letter, add "id" to those that don't...
          name = "id" + name;
        xmlTextWriterWriteAttribute(writer(), C_ID, BAD_CAST name.c_str());
        if(!isalpha(name[0])) //...and write <name> orig title </name>
          {
            xmlTextWriterStartElementNS(writer(), prefix, C_NAME, NULL);
            xmlTextWriterWriteFormatString(writer(),"%s", id);
            xmlTextWriterEndElement(writer());//name
          }
      }

    int TotalCharge = mol.GetTotalCharge();
    if(TotalCharge!=0)
      xmlTextWriterWriteFormatAttribute(writer(), C_FORMALCHARGE, "%d", TotalCharge);

    int TotalSpin = mol.GetTotalSpinMultiplicity();
    if(TotalSpin!=1)
      xmlTextWriterWriteFormatAttribute(writer(), C_SPINMULTIPLICITY, "%d", TotalSpin);

    if(_pxmlConv->IsOption("m") && _pxmlConv->GetOutputIndex()==1) //only on first molecule
      WriteMetadataList(mol);

    pUnitCell = NULL;
    if (!cml1 && mol.HasData(OBGenericDataType::UnitCell))
      {
        WriteCrystal(mol);//Output will be in crystallographic form
        UseFormulaWithNoBonds = false;
      }

    WriteInChI(mol);

    // Create map (tetStereos) from atom Ids to tetstereos
    std::map<unsigned int, OBTetrahedralStereo::Config > tetStereos;
    std::map<unsigned int, OBTetrahedralStereo::Config >::const_iterator tetStereo_cit;
    if (mol.GetDimension()!=3) {
      std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
      for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
        if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
          OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
          // Always get the clockwise version (it's the default anyway) as this has
          // a positive signed volume (i.e. CML atomParity of 1)
          OBTetrahedralStereo::Config cfg = ts->GetConfig(OBStereo::Clockwise);
          if(cfg.specified)
            tetStereos[cfg.center] = cfg;
        }
    }

    vector<string> atomIds;
    if(mol.NumAtoms()>0)
      {
        //if molecule has no bonds and atoms doesn't have coordinates, just output formula
        if(numbonds==0 && UseFormulaWithNoBonds && !mol.Has2D())
          WriteFormula(mol);
        else
          {
            xmlTextWriterStartElementNS(writer(), prefix, C_ATOMARRAY, NULL);

            MakeAtomIds(mol, atomIds);//Pre-construct to take into account atom class data


            stringstream id, eltyp, iso, chg, spn, hct, x, y, z;
            bool anyChg=false, anySpin=false, anyIsotope=false;
            double X, Y, Z; //atom coordinates

            OBAtom *patom;
            vector<OBAtom*>::iterator i;
            for (patom = mol.BeginAtom(i);patom;patom = mol.NextAtom(i))
              {
               string el(OBElements::GetSymbol(patom->GetAtomicNum()));
                if(el=="Xx")
                  el="R";

                int charge = patom->GetFormalCharge();
                int spin = patom->GetSpinMultiplicity();
                int isotope =patom->GetIsotope();

                int hcount=patom->GetImplicitHCount() + patom->ExplicitHydrogenCount(); //includes H isotopes

                X = patom->GetX();
                Y = patom->GetY();
                Z = patom->GetZ();

                if(pUnitCell)
                  {
                    //Convert to fractional coordinates
                    vector3 v = patom->GetVector();
                    v = pUnitCell->CartesianToFractional(v);
                    X = v.x();
                    Y = v.y();
                    Z = v.z();
                    C_X3orFRACT = C_XFRACT;
                    C_Y3orFRACT = C_YFRACT;
                    C_Z3orFRACT = C_ZFRACT;
                    dim=3; //should already be, but make sure
                  }

                if(arrayform)
                  {
                    if(charge)
                      anyChg = true;
                    if(spin)
                      anySpin = true;
                    if(isotope)
                      anyIsotope = true;
                    id << " " << atomIds[patom->GetIdx()];
                    eltyp << " " << el;
                    iso << " " << isotope;
                    chg << " " << charge;
                    spn << " " << spin;
                    hct << " " << hcount;

                    x << " " << X;
                    y << " " << Y;
                    z << " " << Z;
                  }
                else
                  {
                    //Non-array form
                    xmlTextWriterStartElementNS(writer(), prefix, C_ATOM, NULL);
                      xmlTextWriterWriteFormatAttribute(writer(), C_ID,"%s", atomIds[patom->GetIdx()].c_str());

                    if(!cml1)
                      {
                        xmlTextWriterWriteFormatAttribute(writer(), C_ELEMENTTYPE,"%s", el.c_str());
                        if(isotope)
                          xmlTextWriterWriteFormatAttribute(writer(), C_ISOTOPE,"%d", isotope);

                        if(charge)
                          xmlTextWriterWriteFormatAttribute(writer(), C_FORMALCHARGE,"%d", charge);

                        if(spin)
                          xmlTextWriterWriteFormatAttribute(writer(), C_SPINMULTIPLICITY,"%d", spin);

                        xmlTextWriterWriteFormatAttribute(writer(), C_HYDROGENCOUNT,"%d", hcount);

                        if(patom->HasData("label"))
                            xmlTextWriterWriteFormatAttribute(writer(), C_LABEL,"%s",
                            patom->GetData("label")->GetValue().c_str());

                        if(patom->HasData("color"))
                          xmlTextWriterWriteFormatAttribute(writer(), C_COLOR,"%s",
                          patom->GetData("color")->GetValue().c_str());

                        if(patom->HasData("radius"))
                          xmlTextWriterWriteFormatAttribute(writer(), C_RADIUS,"%s",
                          patom->GetData("radius")->GetValue().c_str());

                        if(dim==2)
                          {
                            xmlTextWriterWriteFormatAttribute(writer(), C_X2,"%f", X);
                            xmlTextWriterWriteFormatAttribute(writer(), C_Y2,"%f", Y);
                          }
                        if(dim==3)
                          {
                            xmlTextWriterWriteFormatAttribute(writer(), C_X3orFRACT,"%f", X);
                            xmlTextWriterWriteFormatAttribute(writer(), C_Y3orFRACT,"%f", Y);
                            xmlTextWriterWriteFormatAttribute(writer(), C_Z3orFRACT,"%f", Z);
                          }

                        if( (tetStereo_cit=tetStereos.find(patom->GetId())) != tetStereos.end() )
                          {
                            OBTetrahedralStereo::Config cfg = tetStereo_cit->second;
                            OBStereo::Refs refs = cfg.refs;
                            vector<string> atomrefs;
                            // According to http://cml.sourceforge.net/schema/cmlCore/HTMLDOCS/cmlCore.pdf,
                            // "if there are only 3 ligands, the current atom should be included
                            //  in the 4 atomRefs.".
                            if (cfg.from == OBStereo::ImplicitRef) // e.g. for [S@@](Cl)(Br)I
                                atomrefs.push_back(atomIds[mol.GetAtomById(cfg.center)->GetIdx()]); // Add the central atom again
                              else
                                atomrefs.push_back(atomIds[mol.GetAtomById(cfg.from)->GetIdx()]);

                            for (OBStereo::RefIter ref = refs.begin(); ref!=refs.end(); ++ref) {
                              if ( (OBStereo::Ref)*ref == OBStereo::ImplicitRef) // e.g. for Cl[S@@](Br)I
                                atomrefs.push_back(atomIds[mol.GetAtomById(cfg.center)->GetIdx()]); // Add the central atom again
                              else
                                atomrefs.push_back(atomIds[mol.GetAtomById(*ref)->GetIdx()]);
                            }

                            xmlTextWriterStartElementNS(writer(), prefix, C_ATOMPARITY, NULL);
                            xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS4, "%s %s %s %s",
                              atomrefs[0].c_str(), atomrefs[1].c_str(),
                              atomrefs[2].c_str(), atomrefs[3].c_str());
                            // Set the atomParity - this is always 1 as the atomRefs are arranged
                            // to make this so
                            xmlTextWriterWriteFormatString(writer(), "%d", 1);
                            xmlTextWriterEndElement(writer());//atomParity
                          }
                      }
                    else
                      {
                        //CML1
                        xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "elementType");
                        xmlTextWriterWriteFormatString(writer(),"%s", el.c_str());
                        xmlTextWriterEndElement(writer());

                        if(charge)
                          {
                            xmlTextWriterStartElementNS(writer(), prefix, C_INTEGER, NULL);
                            xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "formalCharge");
                            xmlTextWriterWriteFormatString(writer(),"%d", charge);
                            xmlTextWriterEndElement(writer());
                          }

                        xmlTextWriterStartElementNS(writer(), prefix, C_INTEGER, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "hydrogenCount");
                        xmlTextWriterWriteFormatString(writer(),"%d", hcount);
                        xmlTextWriterEndElement(writer());

                        if(dim==2 || dim==3)
                          {
                            xmlTextWriterStartElementNS(writer(), prefix, C_FLOAT, NULL);
                            xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "x",dim);
                            xmlTextWriterWriteFormatString(writer(),"%f", X);
                            xmlTextWriterEndElement(writer());

                            xmlTextWriterStartElementNS(writer(), prefix, C_FLOAT, NULL);
                            xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "y",dim);
                            xmlTextWriterWriteFormatString(writer(),"%f", Y);
                            xmlTextWriterEndElement(writer());
                          }

                        if(dim==3)
                          {
                            xmlTextWriterStartElementNS(writer(), prefix, C_FLOAT, NULL);
                            xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "z",dim);
                            xmlTextWriterWriteFormatString(writer(),"%f", Z);
                            xmlTextWriterEndElement(writer());
                          }
                        //Stereochemistry currently not written for CML1
                      }
                    xmlTextWriterEndElement(writer());//atom
                  }
              }

            if(arrayform)
              {
                if(!cml1)
                  {
                    xmlTextWriterWriteFormatAttribute(writer(), C_ATOMID,"%s", id.str().c_str());
                    xmlTextWriterWriteFormatAttribute(writer(), C_ELEMENTTYPE,"%s", eltyp.str().c_str());

                    if(anyIsotope)
                      xmlTextWriterWriteFormatAttribute(writer(), C_ISOTOPE,"%s", iso.str().c_str());

                    if(anyChg)
                      xmlTextWriterWriteFormatAttribute(writer(), C_FORMALCHARGE,"%s", chg.str().c_str());

                    if(anySpin)
                      xmlTextWriterWriteFormatAttribute(writer(), C_SPINMULTIPLICITY,"%s", spn.str().c_str());

                    xmlTextWriterWriteFormatAttribute(writer(), C_HYDROGENCOUNT,"%s", hct.str().c_str());

                    if(dim==2)
                      {
                        xmlTextWriterWriteFormatAttribute(writer(), C_X2,"%s", x.str().c_str());
                        xmlTextWriterWriteFormatAttribute(writer(), C_Y2,"%s", y.str().c_str());
                      }
                    if(dim==3)
                      {
                        xmlTextWriterWriteFormatAttribute(writer(), C_X3orFRACT,"%s", x.str().c_str());
                        xmlTextWriterWriteFormatAttribute(writer(), C_Y3orFRACT,"%s", y.str().c_str());
                        xmlTextWriterWriteFormatAttribute(writer(), C_Z3orFRACT,"%s", z.str().c_str());
                      }
                  }
                else
                  {
                    //CML1
                    xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
                    xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomID");
                    xmlTextWriterWriteFormatString(writer(),"%s", id.str().c_str());
                    xmlTextWriterEndElement(writer());

                    xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
                    xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "elementType");
                    xmlTextWriterWriteFormatString(writer(),"%s", eltyp.str().c_str());
                    xmlTextWriterEndElement(writer());

                    if(anyChg)
                      {
                        xmlTextWriterStartElementNS(writer(), prefix, C_INTEGERARRAY, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "formalCharge");
                        xmlTextWriterWriteFormatString(writer(),"%s", chg.str().c_str());
                        xmlTextWriterEndElement(writer());
                      }

                    xmlTextWriterStartElementNS(writer(), prefix, C_INTEGERARRAY, NULL);
                    xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "hydrogenCount");
                    xmlTextWriterWriteFormatString(writer(),"%s", hct.str().c_str());
                    xmlTextWriterEndElement(writer());

                    if(dim==2 || dim==3)
                      {
                        xmlTextWriterStartElementNS(writer(), prefix, C_FLOATARRAY, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "x",dim);
                        xmlTextWriterWriteFormatString(writer(),"%s", x.str().c_str());
                        xmlTextWriterEndElement(writer());

                        xmlTextWriterStartElementNS(writer(), prefix, C_FLOATARRAY, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "y",dim);
                        xmlTextWriterWriteFormatString(writer(),"%s", y.str().c_str());
                        xmlTextWriterEndElement(writer());
                      }
                    if(dim==3)
                      {
                        xmlTextWriterStartElementNS(writer(), prefix, C_FLOATARRAY, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "z",dim);
                        xmlTextWriterWriteFormatString(writer(),"%s", z.str().c_str());
                        xmlTextWriterEndElement(writer());
                      }
                  }
              }
            xmlTextWriterEndElement(writer());//atomArray
          }
      }

    // Create map (ctStereos) from bond Idxs to cistrans stereos
    std::map<unsigned int, OBCisTransStereo* > ctStereos;
    std::map<unsigned int, OBCisTransStereo* >::const_iterator ctStereo_cit;
    if (mol.GetDimension()!=3) {
      std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
      for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
        if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
          OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
          if(ct->GetConfig().specified) {
            unsigned int dblbond = mol.GetBond(mol.GetAtomById(ct->GetConfig().begin),
                                               mol.GetAtomById(ct->GetConfig().end  ))->GetIdx();
            ctStereos[dblbond] = ct;
          }
        }
    }

    if(mol.NumBonds()>0)
      {
        xmlTextWriterStartElementNS(writer(), prefix, C_BONDARRAY, NULL);

        stringstream ord;
        string ref1, ref2;
        OBBond *pbond;
        vector<OBBond*>::iterator ib;
        for (pbond = mol.BeginBond(ib);pbond;pbond = mol.NextBond(ib))
          {
            int bo = pbond->GetBondOrder();

            if(!arrayform)
              {
                if(bo==5 || (WriteAromaticBonds && pbond->IsAromatic())) //aromatic
                  ord << 'A';
                else
                  ord << bo;

                ref1 = atomIds[pbond->GetBeginAtomIdx()];
                ref2 = atomIds[pbond->GetEndAtomIdx()];
                xmlTextWriterStartElementNS(writer(), prefix, C_BOND, NULL);
                //				xmlTextWriterWriteFormatAttribute(writer(), C_ID,"b%d", pbond->GetIdx()); remove bond id
                if(!cml1)
                  {
                    xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS2,"%s %s",
                          ref1.c_str(), ref2.c_str());
                    xmlTextWriterWriteFormatAttribute(writer(), C_ORDER,"%s", ord.str().c_str());

                    if(pbond->HasData("color"))
                      xmlTextWriterWriteFormatAttribute(writer(), C_COLOR,"%s",
                          pbond->GetData("color")->GetValue().c_str());

                    if(pbond->HasData("label"))
                      xmlTextWriterWriteFormatAttribute(writer(), C_LABEL,"%s",
                          pbond->GetData("label")->GetValue().c_str());

                    if( (ctStereo_cit=ctStereos.find(pbond->GetIdx())) != ctStereos.end() )
                      {
                        OBCisTransStereo *ct = ctStereo_cit->second;
                        OBCisTransStereo::Config ct_cfg = ct->GetConfig();

                        // Find a non-implicit ref at either end of the dbl bond
                        OBStereo::Ref beginref, endref;
                        beginref = (ct_cfg.refs[0] == OBStereo::ImplicitRef) ? ct_cfg.refs[1] : ct_cfg.refs[0];
                        endref =   (ct_cfg.refs[2] == OBStereo::ImplicitRef) ? ct_cfg.refs[3] : ct_cfg.refs[2];
                        char cis_or_trans = ct->IsCis(beginref, endref) ? 'C' : 'T';

                        // Prepare the "atomrefs4"
                        vector<string> atomrefs(4);
                        atomrefs[0] = atomIds[mol.GetAtomById(beginref)->GetIdx()];
                        atomrefs[1] = atomIds[mol.GetAtomById(ct_cfg.begin)->GetIdx()];
                        atomrefs[2] = atomIds[mol.GetAtomById(ct_cfg.end)->GetIdx()];
                        atomrefs[3] = atomIds[mol.GetAtomById(endref)->GetIdx()];

                        // Create the XML tags
                        xmlTextWriterStartElementNS(writer(), prefix, C_BONDSTEREO, NULL);
                        xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS4, "%s %s %s %s",
                              atomrefs[0].c_str(), atomrefs[1].c_str(),
                              atomrefs[2].c_str(), atomrefs[3].c_str());
                        xmlTextWriterWriteFormatString(writer(),"%c", cis_or_trans);
                        xmlTextWriterEndElement(writer());//bondStereo
                      }
                  }
                else
                  {
                    //CML1
                    xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
                    xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
                    xmlTextWriterWriteFormatString(writer(),"%s", ref1.c_str());
                    xmlTextWriterEndElement(writer());

                    xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
                    xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
                    xmlTextWriterWriteFormatString(writer(),"%s", ref2.c_str());
                    xmlTextWriterEndElement(writer());

                    xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
                    xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "order");
                    xmlTextWriterWriteFormatString(writer(),"%d", bo);
                    xmlTextWriterEndElement(writer());
                  }
                xmlTextWriterEndElement(writer());//bond
                ord.str(""); //clear (For array form it accumulates.)
              }
            else
              {
                if(bo==5 || (WriteAromaticBonds && pbond->IsAromatic())) //aromatic
                  ord << " " << 'A';
                else
                  ord << " " << bo;

                ref1 += ' ' + atomIds[pbond->GetBeginAtomIdx()];
                ref2 += ' ' + atomIds[pbond->GetEndAtomIdx()];
              }
          }
        if(arrayform)
          {
            if(!cml1)
              {
                xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREF1, "%s", ref1.c_str());
                xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREF2, "%s", ref2.c_str());
                xmlTextWriterWriteFormatAttribute(writer(), C_ORDER, "%s", ord.str().c_str());
              }
            else
              {
                //CML1
                xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
                xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
                xmlTextWriterWriteFormatString(writer(),"%s", ref1.c_str());
                xmlTextWriterEndElement(writer());

                xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
                xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
                xmlTextWriterWriteFormatString(writer(),"%s", ref2.c_str());
                xmlTextWriterEndElement(writer());

                xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
                xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "order");
                xmlTextWriterWriteFormatString(writer(),"%s", ord.str().c_str());
                xmlTextWriterEndElement(writer());
              }
          }

        xmlTextWriterEndElement(writer());//bondArray

        //When array form, write bondStereo here
        if(arrayform)
          {
            for (pbond = mol.BeginBond(ib);pbond;pbond = mol.NextBond(ib))
              {
                if(pbond->GetBondOrder()==2 || pbond->IsWedge() || pbond->IsHash())
                  WriteBondStereo(pbond, atomIds);
              }
          }
      }

    bool propertyListWritten=false;
    if(mol.HasData(ThermoData))
      WriteThermo(mol, propertyListWritten);

    if(_pxmlConv->IsOption("p"))
      WriteProperties(mol, propertyListWritten);

    if(propertyListWritten)
      xmlTextWriterEndElement(writer());//propertList

    xmlTextWriterEndElement(writer());//molecule

    //Note that nothing will be written unless the next block is executed
    //IsLast() MUST return true for the last molecule.
    if(!_pxmlConv->IsOption("MolsNotStandalone") && _pxmlConv->IsLast())
      {
        xmlTextWriterEndDocument(writer());
        OutputToStream();
      }
    return true;
  }

  ///Constructs a unique id for each atom.
  void CMLFormat::MakeAtomIds(OBMol& mol, vector<string>& atomIDs)
  {
    /* If there is no atom class data for the atom, the id is a followed by the atom index.
       If there is atom class data, an underscore is appended followed by the atom class.
     */

    stringstream ss;
    atomIDs.push_back("Error"); //atom idex stats at 1. atomIDs[0] is not used
    for (unsigned int idx=1; idx<=mol.NumAtoms(); ++idx)
    {
      ss.str("");
      ss << 'a' << idx;
      OBGenericData* pac = mol.GetAtom(idx)->GetData("Atom Class");
      if(pac)
      {
        OBPairInteger* acdata = dynamic_cast<OBPairInteger*>(pac);
        if (acdata) {
          int ac = acdata->GetGenericValue();
          if (ac >= 0) { // Allow 0, why not?
            ss << '_' << ac;
          }
        }
      }
      atomIDs.push_back(ss.str());
    }
  }

  void CMLFormat::WriteFormula(OBMol mol)
  {
    //mol is a copy
    static const xmlChar C_FORMULA[] = "formula";
    static const xmlChar C_CONCISE[] = "concise";
    if(mol.NumAtoms()==1)
      mol.AddHydrogens(false,false);
    xmlTextWriterStartElementNS(writer(), prefix, C_FORMULA, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_CONCISE,"%s", mol.GetSpacedFormula().c_str());
    xmlTextWriterEndElement(writer());//formula
  }

  void CMLFormat::WriteBondStereo(OBBond* pbond, vector<string>& atomIDs)
  {
    static const xmlChar C_ATOMREFS4[]  = "atomRefs4";
    static const xmlChar C_BONDSTEREO[] = "bondStereo";

    char ch=0;
    if(pbond->IsWedge())
      ch='W';
    else if(pbond->IsHash())
      ch='H';

    if(ch)
      //this line here because element may not be written with double bond
      xmlTextWriterStartElementNS(writer(), prefix, C_BONDSTEREO, NULL);
    else
    {
      return; // TODO: This code has bit-rotted
      //double bond stereo
      int ud1=0, ud2=0;
      int idx1=0, idx2=0;
      OBAtom* patomA = pbond->GetBeginAtom();
      //FOR_BONDS_OF_ATOM(b1,patomA)
      //  {
      //    if(b1->IsUp() || b1->IsDown() )
      //      {
      //        idx1=(b1->GetNbrAtom(patomA))->GetIdx();
      //        ud1 = b1->IsDown() ? -1 : 1;
      //        // Conjugated double bonds have to be treated differently, see comments
      //        // in OBMol2Smi::GetCisTransBondSymbol(). Reverse symbol for other than first double bond.
      //        if((b1->GetNbrAtom(patomA))->HasDoubleBond())
      //          ud1 = -ud1;
      //        break;
      //      }
      //  }
      OBAtom* patomB = pbond->GetEndAtom();
      //FOR_BONDS_OF_ATOM(b2,patomB)
      //  {
      //    if(b2->IsUp() || b2->IsDown() )
      //      {
      //        idx2=(b2->GetNbrAtom(patomB))->GetIdx();
      //        ud2 = b2->IsDown() ? -1 : 1;
      //        break;
      //      }
      //  }
      if(!ud1 || !ud2)
        return;

      xmlTextWriterStartElementNS(writer(), prefix, C_BONDSTEREO, NULL);
      xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS4, "%s %s %s %s",
//                        "a%d a%d a%d a%d", idx1, patomA->GetIdx(), patomB->GetIdx(), idx2);
            atomIDs[idx1].c_str(), atomIDs[patomA->GetIdx()].c_str(),
            atomIDs[patomB->GetIdx()].c_str(), atomIDs[idx2].c_str());
      ch = (ud1==ud2) ? 'C' : 'T';
    }

    xmlTextWriterWriteFormatString(writer(),"%c", ch);
    xmlTextWriterEndElement(writer());//bondStereo
  }

  void CMLFormat::WriteCrystal(OBMol& mol)
  {
    static const xmlChar C_CRYSTAL[]  = "crystal";
    static const xmlChar C_SCALAR[] = "scalar";
    // static const xmlChar C_Z[] = "z";
    static const xmlChar C_TITLE[] = "title";
    static const xmlChar C_UNITS[] = "units";
    static const xmlChar C_SYMMETRY[]  = "symmetry";
    static const xmlChar C_SPACEGROUP[]  = "spaceGroup";
    static const xmlChar C_TRANSFORM3[]  = "transform3";

    pUnitCell = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);

    xmlTextWriterStartElementNS(writer(), prefix, C_CRYSTAL, NULL);
    //	xmlTextWriterWriteFormatAttribute(writer(), C_z,"%d", number of molecules per cell);

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "a");
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:angstrom");
    xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetA());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "b");
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:angstrom");
    xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetB());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "c");
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:angstrom");
    xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetC());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "alpha");
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:degree");
    xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetAlpha());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "beta");
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:degree");
    xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetBeta());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "gamma");
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:degree");
    xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetGamma());
    xmlTextWriterEndElement(writer());//scalar

    const SpaceGroup *group = pUnitCell->GetSpaceGroup();
    string s;
    if (group)
	  {
        xmlTextWriterStartElementNS(writer(), prefix, C_SYMMETRY, NULL);
        xmlTextWriterWriteAttribute (writer(), C_SPACEGROUP, (const xmlChar*)group->GetHallName().c_str());
        transform3dIterator ti;
        const transform3d *t = group->BeginTransform(ti);
        string s;
        while(t)
          {
			s = t->DescribeAsValues() + " 0 0 0 1";
            xmlTextWriterWriteElement(writer(), C_TRANSFORM3, (const xmlChar*)s.c_str());
            t = group->NextTransform(ti);
          }
        xmlTextWriterEndElement(writer());//symmetry
	  }
    else
	  {
        //s = pUnitCell.GetSpaceGroupName();
        s = pUnitCell->GetSpaceGroupName();
        if (s.length())
	      {
            xmlTextWriterStartElementNS(writer(), prefix, C_SYMMETRY, NULL);
            xmlTextWriterWriteAttribute (writer(), C_SPACEGROUP, (const xmlChar*)s.c_str());
            xmlTextWriterEndElement(writer());//symmetry
	      }
	  }

    xmlTextWriterEndElement(writer());//crystal
  }

  void CMLFormat::WriteProperties(OBMol& mol, bool& propertyListWritten)
  {
    static const xmlChar C_DICTREF[]      = "dictRef";
    static const xmlChar C_PROPERTYLIST[] = "propertyList";
    static const xmlChar C_PROPERTY[]     = "property";
    static const xmlChar C_SCALAR[]       = "scalar";
    static const xmlChar C_TITLE[]        = "title";

    vector<OBGenericData*>::iterator k;
    vector<OBGenericData*> vdata = mol.GetData();
    for (k = vdata.begin();k != vdata.end();++k)
      {
        if  ((*k)->GetDataType() == OBGenericDataType::PairData
          && (*k)->GetOrigin()   != local //internal OBPairData is not written
          && (*k)->GetAttribute()!= "InChI" //InChI is output in <identifier>
          && (*k)->GetAttribute()!= "PartialCharges")//annotation not needed since partial charges are not output in this format
          {
            if(!propertyListWritten)
              {
                xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTYLIST, NULL);
                propertyListWritten=true;
              }

            xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTY, NULL);
            //Title is now on <property>. If the attribute name has a namespace, use dictRef instead.
            string att((*k)->GetAttribute());
            xmlTextWriterWriteFormatAttribute(writer(),
              (att.find(':')==string::npos) ? C_TITLE : C_DICTREF,
              "%s",att.c_str());
            xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);

            //Title used to be on <scalar>...
            //xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s",(*k)->GetAttribute().c_str());
            xmlTextWriterWriteFormatString(writer(),"%s", (static_cast<OBPairData*>(*k))->GetValue().c_str());
            xmlTextWriterEndElement(writer());//scalar
            xmlTextWriterEndElement(writer());//property
          }
      }

    static const double CALSTOJOULES = 4.184;
    //Energy is output when it is not zero
    //This is the molecular energy, probably originally in Hartrees,
    // stored in OB as kcal/mol, but output here in kJ/mol
    if(fabs(mol.GetEnergy()) > 1e-3)
      WriteScalarProperty(mol, "Energy", mol.GetEnergy() * CALSTOJOULES,
        "me:ZPE", "kJ/mol", "computational");

    //spinMultiplicity is written only when it is not 1
    int smult = mol.GetTotalSpinMultiplicity();
    if(smult!=1)
      WriteScalarProperty(mol, "SpinMultiplicity", smult, "me:spinMultiplicity");

    if(mol.HasData(OBGenericDataType::VibrationData))
      WriteVibrationData(mol);
    if(mol.HasData(OBGenericDataType::RotationData))
      WriteRotationData(mol);
  }

  void CMLFormat::WriteThermo(OBMol& mol, bool& propertyListWritten)
  {
    static const xmlChar C_PROPERTYLIST[] = "propertyList";
    static const xmlChar C_PROPERTY[]     = "property";
    static const xmlChar C_SCALAR[]       = "scalar";
    static const xmlChar C_ARRAY[]        = "array";
    static const xmlChar C_DICTREF[]      = "dictRef";
    static const xmlChar C_SIZE[]         = "size";

    OBNasaThermoData* pThermoData = static_cast<OBNasaThermoData*>(mol.GetData(ThermoData));

    if(!propertyListWritten)
      {
        xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTYLIST, NULL);
        propertyListWritten=true;
      }

    xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","Thermo_OldNasa");

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","NasaLowT");
    xmlTextWriterWriteFormatString(writer(),"%.1f", pThermoData->GetLoT());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","NasaHighT");
    xmlTextWriterWriteFormatString(writer(),"%.1f", pThermoData->GetHiT());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","NasaMidT");
    xmlTextWriterWriteFormatString(writer(),"%.1f", pThermoData->GetMidT());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","Phase");
    xmlTextWriterWriteFormatString(writer(),"%c", pThermoData->GetPhase());
    xmlTextWriterEndElement(writer());//scalar

    xmlTextWriterStartElementNS(writer(), prefix, C_ARRAY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","NasaCoeffs");
    xmlTextWriterWriteFormatAttribute(writer(), C_SIZE,"%d",14);
    for(int i=0;i<14;++i)
      xmlTextWriterWriteFormatString(writer()," %e", pThermoData->GetCoeff(i));
    xmlTextWriterEndElement(writer());//array

    xmlTextWriterEndElement(writer());//property
  }

  std::string getSeparator()
  {
#ifdef WIN32
    return "\\";
#else
    return "/";
#endif
  }

  ///Returns molecule title or molecule number if there is no title together with the file name
  string CMLFormat::GetMolID()
  {
    stringstream molID;
    if(strlen(_pmol->GetTitle())==0)
      molID << "Mol #" << _pxmlConv->GetOutputIndex()+1;
    else
      molID << _pmol->GetTitle();

    string fn(_pxmlConv->GetInFilename());
    //Get file name: remove path
    string::size_type pos = fn.rfind(getSeparator());
    if(pos!=string::npos)
      fn.erase(0,pos+1);
    molID << " (in " << fn << ')';
    return molID.str();
  }

  bool CMLFormat::WriteInChI(OBMol& mol)
  {
    //If OBPair data has an entry with attribute "inchi" it is not
    //output in the property list but as a separate element in the form:
    //<identifier convention="iupac:inchi" value="InChI=1/CH4/h1H4"/>
    static const xmlChar C_IDENTIFIER[] = "identifier";
    static const xmlChar C_CONVENTION[] = "convention";
    static const xmlChar C_VALUE[]      = "value";
    OBPairData* pData = dynamic_cast<OBPairData*>(mol.GetData("InChI"));
    if(pData)
    {
      xmlTextWriterStartElementNS(writer(), prefix, C_IDENTIFIER, NULL);
      xmlTextWriterWriteFormatAttribute(writer(), C_CONVENTION,"%s","iupac:inchi");
      xmlTextWriterWriteFormatAttribute(writer(), C_VALUE,"%s", pData->GetValue().c_str());
      xmlTextWriterEndElement(writer());//identifier
      return true;
    }
    return false; //not written
  }

  bool CMLFormat::WriteVibrationData(OBMol& mol)
  {
    static const xmlChar C_PROPERTY[]     = "property";
    static const xmlChar C_SCALAR[]       = "scalar";
    static const xmlChar C_ARRAY[]        = "array";
    static const xmlChar C_DICTREF[]      = "dictRef";
    static const xmlChar C_UNITS[]        = "units";
    static const xmlChar C_TITLE[]        = "title";

    OBVibrationData* vd = (OBVibrationData*)mol.GetData(OBGenericDataType::VibrationData);

    xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s","Vibrational Frequencies");
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","me:vibFreqs");

    xmlTextWriterStartElementNS(writer(), prefix, C_ARRAY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s","cm-1");

    double imaginaryFrequency = 0.0;
    //A negative frequency is output separately as an imaginary frequency (for transition states)
    for (unsigned int i=0; i<vd->GetNumberOfFrequencies(); ++i)
    {
      double freq = vd->GetFrequencies()[i];
      if(freq>0.0)
        xmlTextWriterWriteFormatString(writer(),"%.2lf ", freq);
      else
        imaginaryFrequency = -freq;
    }
    xmlTextWriterEndElement(writer());//array
    xmlTextWriterEndElement(writer());//property

    if(imaginaryFrequency>0.0)
      WriteScalarProperty(mol, "ImaginaryFrequency", imaginaryFrequency, "me:imFreqs", "cm-1");

    return true;
  }

  bool CMLFormat::WriteRotationData(OBMol& mol)
  {
    static const xmlChar C_PROPERTY[]     = "property";
    static const xmlChar C_SCALAR[]       = "scalar";
    static const xmlChar C_ARRAY[]        = "array";
    static const xmlChar C_DICTREF[]      = "dictRef";
    static const xmlChar C_UNITS[]        = "units";
    static const xmlChar C_TITLE[]        = "title";

    OBRotationData* rd = (OBRotationData*)mol.GetData(OBGenericDataType::RotationData);

    xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s","Rotational Constants");
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","me:rotConsts");

    xmlTextWriterStartElementNS(writer(), prefix, C_ARRAY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s","cm-1");
    const double WAVENUM_TO_GHZ=30.0;
    for (unsigned int i=0; i<rd->GetRotConsts().size(); ++i)
      if(rd->GetRotConsts()[i]!=0.0)
        xmlTextWriterWriteFormatString(writer(),"%.3lf ", rd->GetRotConsts()[i]/WAVENUM_TO_GHZ);
    xmlTextWriterEndElement(writer());//array
    xmlTextWriterEndElement(writer());//property
    xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s","Symmetry Number");
    xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s","me:symmetryNumber");

    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    xmlTextWriterWriteFormatString(writer(),"%d ", rd->GetSymmetryNumber());
    xmlTextWriterEndElement(writer());//scalar
    xmlTextWriterEndElement(writer());//property
    return true;
  }


  bool CMLFormat::WriteScalarProperty(OBMol& mol,
    const char* title, double value, const char* dictref, const char* units, const char* convention)
  {
    static const xmlChar C_PROPERTY[]     = "property";
    static const xmlChar C_SCALAR[]       = "scalar";
    static const xmlChar C_DICTREF[]      = "dictRef";
    static const xmlChar C_UNITS[]        = "units";
    static const xmlChar C_TITLE[]        = "title";
    static const xmlChar C_CONVENTION[]   = "convention";
    static const xmlChar C_ZPEADDED[]   = "zeroPointVibEnergyAdded";

    xmlTextWriterStartElementNS(writer(), prefix, C_PROPERTY, NULL);
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s",title);
    if(dictref)
      xmlTextWriterWriteFormatAttribute(writer(), C_DICTREF,"%s",dictref);
    xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
    if(units)
      xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s",units);
    if(convention)
    {
      xmlTextWriterWriteFormatAttribute(writer(), C_CONVENTION,"%s",convention);
      if(strcmp(convention, "computational")==0)
        xmlTextWriterWriteFormatAttribute(writer(), C_ZPEADDED,"%s","false");
    }
    xmlTextWriterWriteFormatString(writer(),"%.2lf ", value);
    xmlTextWriterEndElement(writer());//scalar
    xmlTextWriterEndElement(writer());//property
    return true;
  }

  bool CMLFormat::WriteChemObject(OBConversion* pConv)
  {
    int OIndex = pConv->GetOutputIndex();
    OBBase* pOb = pConv->GetChemObject();
    if(dynamic_cast<OBMol*> (pOb))
    {
      //With an OBMol object, do the same as if this function wasn't defined,
      //i.e.access the functionality in OBMoleculeFormat

      //restore output index which is (unhelpfully) incremented by GetChemObject
      pConv->SetOutputIndex(OIndex);
      return XMLMoleculeFormat::WriteChemObject(pConv);
    }

    //With OBReaction object, handle directly in CMLFormat::WriteMolecule
    bool ret = WriteMolecule(pOb,pConv);
    delete pOb;
    return ret;
  }

}//namespace
