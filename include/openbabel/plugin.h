/**********************************************************************
plugin.h - facilitates construction of plugin classes

Copyright (C) 2007 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#ifndef OB_PLUGIN_H
#define OB_PLUGIN_H

#include <openbabel/babelconfig.h>
#include <openbabel/dlhandler.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <cstring>

#ifndef OBERROR
 #define OBERROR
#endif

namespace OpenBabel
{
///@addtogroup plugins Plugins
///@{

/// @brief Case insensitive string comparison for PluginMapType key.
struct OBERROR CharPtrLess : public std::binary_function<const char*,const char*, bool>
{
  bool operator()(const char* p1,const char* p2) const
  { return strcasecmp(p1,p2)<0; }
};

/** \class OBPlugin plugin.h <openbabel/plugin.h>
    \brief Base class for all types of dynamic classes discovered at runtime
    \since version 2.2
    \sa plugin.h
 */
class OBERROR OBPlugin
{
public:

  //Maps of thistype are used to store
  // (a)a list of the plugin types in OBPlugin, and
  // (b)a list of the sub-types in each type class derived from OBPlugin.
  typedef std::map<const char*, OBPlugin*, CharPtrLess> PluginMapType;
  typedef PluginMapType::const_iterator PluginIterator;

  ///Virtual destructor necessary for classes with virtual functions
  virtual ~OBPlugin(){};

  ///Required description of a sub-type
  virtual const char* Description() { return nullptr;} ;

  ///Redefined by each plugin type: "formats", "fingerprints", etc.
  virtual const char* TypeID(){ return "plugins"; }

  ///Write information on a plugin class to the string txt.
  ///Return false if not written.
  ///The default implementation outputs:
  /// the ID, a tab character, and the first line of the Description.
  ///The param string can be used in derived types to provide different outputs.
  virtual bool Display(std::string&txt, const char* param, const char* ID=nullptr);

  ///Make a new instance of the class.
  ///See OpTransform, OBGroupContrib, SmartsDescriptor classes for derived versions.
  ///Usually, the first parameter is the classname, the next three are
  ///parameters(ID, filename, description) for a constructor, and the rest data.
  virtual OBPlugin* MakeInstance(const std::vector<std::string>&){return nullptr;}

  ///Initialize the plugin.
  ///The default version does nothing.
  virtual void Init(){};

  ///Get a pointer to a plugin from its type and ID. Return NULL if not found.
  ///If Type is NULL, search all types. Not cast to Type*
  static OBPlugin* GetPlugin(const char* Type, const char* ID);

  ///Return the ID of the sub-type instance.
  const char* GetID()const{return _id;};

  ///Output a list of sub-type classes of the the type PluginID,
  ///or, if PluginID is "plugins" or empty, a list of the base types.
  ///If PluginID is not recognized or is NULL, the base types are output and the return is false.
  static bool ListAsVector(const char* PluginID, const char* param, std::vector<std::string>& vlist);

  ///As ListAsVector but sent to an ostream with a default of cout if not specified
  static void List(const char* PluginID, const char* param=nullptr, std::ostream* os=&std::cout);

  ///As ListAsVector but returns a string containing the list
  static std::string ListAsString(const char* PluginID, const char* param=nullptr);

  ///Utility function to return only the first line of a string
  static std::string FirstLine(const char* txt);

  ///Return an iterator at the start of the map of the plugin types PluginID or,
  /// if there is no such map, the end of the top level plugin map.
  static PluginIterator Begin(const char* PluginID)
  {
    if( !strcmp(PluginID, "plugins") || GetTypeMap(PluginID)!=PluginMap())
      return GetTypeMap(PluginID).begin();
    else
      return PluginMap().end();
  }

  static PluginIterator End(const char* PluginID)
  {
    return GetTypeMap(PluginID).end();
  }

  ///Returns the map of the subtypes
  virtual PluginMapType& GetMap() const =0;

  ///Load all plugins (formats, fingerprints, forcefields etc.)
  static void LoadAllPlugins();

protected:
  ///\brief Returns a reference to the map of the plugin types.
  /// Is a function rather than a static member variable to avoid initialization problems.
  static PluginMapType& PluginMap()
  {
    static PluginMapType m;
    return m;
  }

  ///Keep a record if all plugins have been loaded
  static int AllPluginsLoaded;

  ///Returns the map of a particular plugin type, e.g. GetMapType("fingerprints")
  static PluginMapType& GetTypeMap(const char* PluginID);

  ///\brief Returns the type with the specified ID, or NULL if not found.
  ///Needs to be cast to the appropriate class in the calling routine.
  static OBPlugin* BaseFindType(PluginMapType& Map, const char* ID);

protected:
  const char* _id;
};

#if defined(__CYGWIN__) || defined(__MINGW32__)

//Macro to be added to definition of the base class
#define MAKE_PLUGIN(BaseClass)\
protected:\
  static PluginMapType& Map();\
  virtual PluginMapType& GetMap() const {\
    return Map();\
  }\
public:\
  static BaseClass*& Default() {\
    static BaseClass* d;\
    return d;\
  }\
  BaseClass(const char* ID, bool IsDefault=false) {\
    _id=ID;\
    if (ID&&*ID) {\
      if (IsDefault || Map().empty()) {\
        Default() = this;\
      }\
      if (Map().count(ID) == 0) {\
        Map()[ID] = this;\
        PluginMap()[TypeID()] = this;\
      }\
    }\
  }\
  static BaseClass* FindType(const char* ID) {\
    if (!ID || *ID==0 || *ID==' ') {\
      return Default();\
    }\
    return static_cast<BaseClass*>(BaseFindType(Map(),ID));\
  }

#define PLUGIN_CPP_FILE(BaseClass)\
OBPlugin::PluginMapType& BaseClass::Map() {\
  static OBPlugin::PluginMapType map;\
  return map;\
}

#else // __CYGWIN__ || __MINGW32__

//Macro to be added to definition of the base class
#define MAKE_PLUGIN(BaseClass)\
protected:\
  static PluginMapType& Map() {\
    static PluginMapType m;\
    return m;\
  }\
  virtual PluginMapType& GetMap() const {\
    return Map();\
  }\
public:\
  static BaseClass*& Default() {\
    static BaseClass* d;\
    return d;\
  }\
  BaseClass(const char* ID, bool IsDefault=false) {\
    _id=ID;\
    if (ID&&*ID) {\
      if (IsDefault || Map().empty()) {\
        Default() = this;\
      }\
      if (Map().count(ID) == 0) {\
        Map()[ID] = this;\
        PluginMap()[TypeID()] = this;\
      }\
    }\
  }\
  static BaseClass* FindType(const char* ID) {\
    if (!ID || *ID==0 || *ID==' ') {\
      return Default();\
    }\
    return static_cast<BaseClass*>(BaseFindType(Map(),ID));\
  }

#endif // __CYGWIN__ || __MINGW32__

/** \file plugin.h
   \brief Simplify 'plugin' classes to be discovered and/or loaded at runtime.

The code in this file makes it easy to make 'plugin' classes. These classes are
derived from a base class, like OBFingerprint. The derived classes
('sub-types' like fingerprint2) usually have a single instance. Plugin classes
are only discovered at runtime, so no existing code needs to be changed when
adding a new derived class. In some builds the new code can be added or removed
by just moving a DLL or so file. The plugin classes derived from any base class
(including new ones) type can be listed from the commandline.

<h2>Step-by-Step Instructions</h2>

1) In the header file for YourBaseClass (which handles whatsits).
Make sure to include the plugin.h header ,
derive the class from OBPlugin
and in its definition add
  the MAKE_PLUGIN macro
  and a function TypeID() containing a simple descriptor of the type
\code
#include <openbabel/plugin.h>
class YourBaseClass : public OBPlugin
{
  MAKE_PLUGIN(YourBaseClass)

  const char* TypeID()
  { return "whatsits"; };

  ...rest of implementation, probably involving virtual functions redefined
  in the sub-type classes
};
\endcode
See below for what the macro contains.

2) Declare each sub-type in a class derived from the base class
and give it a constructor which calls OBPlugin constructor as shown:
\code
class YourSubType1 : public YourBaseClass
{
public:
  YourSubtype1(const char* ID, bool IsDefault=false)
    : YourBaseClass(ID, IsDefault){}

  virtual string Description()
  { return "A description with one or more lines";};

  ...rest of implementation
};
\endcode
Only the first line of the description is used when the subclasses are listed.

3) Declare a global instance of the sub-type class which specifies its ID.
and, optionally, whether it is to be regarded as the default type of YourBaseClass.
\code
YourSubType1 theType1("whatsit2",true);
\endcode

4) The following functions are available:

YourBaseClass* YourBaseClass::FindType(const char* ID);
This returns the default type when ID is NULL or empty.

To list the sub-types of any plugin class use the List which sends to cout
by default (or any other ostream if specified).
\code
  OBPlugin::List("whatsits")
\endcode
The ListAsString and ListAsVector functions are alternatives, usable with scripting.

It is also possible to iterate through each sub-type by the following code:
\code
  OBPlugin::PluginIterator itr;
  for(itr=OBPlugin::Begin("whatsits");itr!=OBPlugin::End("whatsits");++itr)
  {
    itr is a std::map::const_iterator
    itr->first is the ID of the subtype;
    itr->second is The OBPlugin* which you will have to cast to your type
  }
\endcode
Since this is not the most beautiful code, it is probably better to use the
List methods if possible.

YourBaseClass* MakeNewInstance();



<h2>How it works</h2>

MAKE_PLUGIN(YourBaseClass) inserts the following code into YourBaseClass:
\code
protected:

  //The collection of sub-types is in a local static variable to avoid
  //any difficulties with the order of initialization of static objects.
  static PluginMapType& Map()
  {
    static PluginMapType m;
    return m;
  }

  //Making the map accessible to the base class (Cannot be used during construction)
  virtual PluginMapType& GetMap()const
  {
   return Map();
  }

public:
  static YourBaseClass*& Default()
  {
    static YourBaseClass* d;
    return d;
  }

  //Constructor registers the sub-type
  YourBaseClass(const char* ID, bool IsDefault=false)
  {
    _id = ID;
    if(ID && *ID && *ID!=' ') //do not register if ID is empty or starts with a space
    {
      if(IsDefault || Map().empty())
        Default() = this;
      Map()[ID]=this;
      //Ensure YourBaseClass is registered in OBPlugin so it can be accessed from the commandline
      PluginMap()[TypeID()] =this;
    }
  }

  ///Returns the sub-type associated with the ID, or the default subtype if ID NULL or empty.
  static YourBaseClass* FindType(const char* ID)
  {
    if(!ID || *ID==0)
      return Default();
    return static_cast<YourBaseClass*>(BaseFindType(Map(),ID));
  }
\endcode
*/
/*\@}*/

#ifndef SWIG // Skipped by SWIG (for the moment)

#ifndef USING_DYNAMIC_LIBS

#define OB_STATIC_PLUGIN(className,instanceName) \
  class className; \
  OBAPI OB_EXTERN className instanceName;

  // formats
  OB_STATIC_PLUGIN(ABINITFormat, theABINITFormat)
  OB_STATIC_PLUGIN(AcesOutputFormat, theAcesOutputFormat)
  OB_STATIC_PLUGIN(AcesInputFormat, theAcesInputFormat)
  OB_STATIC_PLUGIN(ACRFormat, theACRFormat)
  OB_STATIC_PLUGIN(ADFOutputFormat, theADFOutputFormat)
  OB_STATIC_PLUGIN(ADFInputFormat, theADFInputFormat)
  OB_STATIC_PLUGIN(AlchemyFormat, theAlchemyFormat)
  OB_STATIC_PLUGIN(AmberPrepFormat, theAmberPrepFormat)
  OB_STATIC_PLUGIN(AoforceFormat, theAoforceFormat)
  OB_STATIC_PLUGIN(OBAPIInterface, theOBAPIInterface)
  OB_STATIC_PLUGIN(BallStickFormat, theBallStickFormat)
  OB_STATIC_PLUGIN(BGFFormat, theBGFFormat)
  OB_STATIC_PLUGIN(BoxFormat, theBoxFormat)
  OB_STATIC_PLUGIN(CacaoFormat, theCacaoFormat)
  OB_STATIC_PLUGIN(CacheFormat, theCacheFormat)
  OB_STATIC_PLUGIN(CARFormat, theCARFormat)
  OB_STATIC_PLUGIN(CASTEPFormat, theCASTEPFormat)
  OB_STATIC_PLUGIN(CCCFormat, theCCCFormat)
  OB_STATIC_PLUGIN(CHEM3D1Format, theCHEM3D1Format)
  OB_STATIC_PLUGIN(CHEM3D2Format, theCHEM3D2Format)
  OB_STATIC_PLUGIN(ChemDrawBinaryXFormat, theChemDrawBinaryXFormat)
  OB_STATIC_PLUGIN(ChemDrawFormat, theChemDrawFormat)
  OB_STATIC_PLUGIN(ChemKinFormat, theChemKinFormat)
  OB_STATIC_PLUGIN(CHTFormat, theCHTFormat)
  OB_STATIC_PLUGIN(CIFFormat, theCIFFormat)
  OB_STATIC_PLUGIN(CopyFormat, theCopyFormat)
  OB_STATIC_PLUGIN(CRK2DFormat, theCRK2DFormat)
  OB_STATIC_PLUGIN(CRK3DFormat, theCRK3DFormat)
  OB_STATIC_PLUGIN(CSRFormat, theCSRFormat)
  OB_STATIC_PLUGIN(CSSRFormat, theCSSRFormat)
  OB_STATIC_PLUGIN(DlpolyConfigFormat, theDlpolyConfigFormat)
  OB_STATIC_PLUGIN(DlpolyHISTORYFormat, theDlpolyHISTORYFormat)
  OB_STATIC_PLUGIN(DMolFormat, theDMolFormat)
  OB_STATIC_PLUGIN(EXYZFormat, theEXYZFormat)
  OB_STATIC_PLUGIN(FASTAFormat, theFASTAFormat)
  OB_STATIC_PLUGIN(FastSearchFormat, theFastSearchFormat)
  OB_STATIC_PLUGIN(FCHKFormat, theFCHKFormat)
  OB_STATIC_PLUGIN(FEATFormat, theFEATFormat)
  OB_STATIC_PLUGIN(FenskeZmatFormat, theFenskeZmatFormat)
  OB_STATIC_PLUGIN(FHIaimsFormat,theFHIaimsFormat)
  OB_STATIC_PLUGIN(FingerprintFormat, theFingerprintFormat)
  OB_STATIC_PLUGIN(FreeFormFractionalFormat, theFreeFormFractionalFormat)
  OB_STATIC_PLUGIN(GAMESSOutputFormat, theGAMESSOutputFormat)
  OB_STATIC_PLUGIN(GAMESSInputFormat, theGAMESSInputFormat)
  OB_STATIC_PLUGIN(OBGaussianCubeFormat, theGaussianCubeFormat)
  OB_STATIC_PLUGIN(GaussianOutputFormat, theGaussianOutputFormat)
  OB_STATIC_PLUGIN(GaussianInputFormat, theGaussianInputFormat)
  OB_STATIC_PLUGIN(GaussianZMatrixInputFormat, theGaussianZMatrixInputFormat)
  OB_STATIC_PLUGIN(GenBankFormat, theGenBankFormat)
  OB_STATIC_PLUGIN(GhemicalFormat, theGhemicalFormat)
  OB_STATIC_PLUGIN(GROFormat, theGROFormat)
  OB_STATIC_PLUGIN(GROMOS96Format, theGROMOS96Format)
  OB_STATIC_PLUGIN(GULPFormat, theGULPFormat)
  OB_STATIC_PLUGIN(HINFormat, theHINFormat)
  OB_STATIC_PLUGIN(JaguarOutputFormat, theJaguarOutputFormat)
  OB_STATIC_PLUGIN(JaguarInputFormat, theJaguarInputFormat)
  OB_STATIC_PLUGIN(LMPDATFormat, theLMPDATFormat)
  OB_STATIC_PLUGIN(MCDLFormat, theMCDLFormat)
  OB_STATIC_PLUGIN(MOLFormat, theMOLFormat)
  OB_STATIC_PLUGIN(SDFormat, theSDFormat)
  OB_STATIC_PLUGIN(OBT41Format, t41Format__)
  OB_STATIC_PLUGIN(OBMoldenFormat, moldenFormat__)
  OB_STATIC_PLUGIN(mmCIFFormat, themmCIFFormat)
  OB_STATIC_PLUGIN(MacroModFormat, theMacroModFormat)
  OB_STATIC_PLUGIN(MNAFormat, theMNAFormat)
  OB_STATIC_PLUGIN(MOL2Format, theMOL2Format)
  OB_STATIC_PLUGIN(MolproOutputFormat, theMolproOutputFormat)
  OB_STATIC_PLUGIN(MolproInputFormat, theMolproInputFormat)
  OB_STATIC_PLUGIN(MolReportFormat, theMolReportFormat)
  OB_STATIC_PLUGIN(MOPACFormat, theMOPACFormat)
  OB_STATIC_PLUGIN(MOPACCARTFormat, theMOPACCARTFormat)
  OB_STATIC_PLUGIN(MOPACINTFormat, theMOPACINTFormat)
  OB_STATIC_PLUGIN(MPDFormat, theMPDFormat)
  OB_STATIC_PLUGIN(MPQCFormat, theMPQCFormat)
  OB_STATIC_PLUGIN(MPQCInputFormat, theMPQCInputFormat)
  OB_STATIC_PLUGIN(MSIFormat, theMSIFormat)
  OB_STATIC_PLUGIN(OBMSMSFormat, msmsFormat__)
  OB_STATIC_PLUGIN(NulFormat, theNulFormat)
  OB_STATIC_PLUGIN(NWChemOutputFormat, theNWChemOutputFormat)
  OB_STATIC_PLUGIN(NWChemInputFormat, theNWChemInputFormat)
  OB_STATIC_PLUGIN(OBOpenDXCubeFormat, theOpenDXCubeFormat)
  OB_STATIC_PLUGIN(OrcaOutputFormat, theOrcaOutputFormat)
  OB_STATIC_PLUGIN(OrcaInputFormat, theOrcaInputFormat)
  OB_STATIC_PLUGIN(OutputFormat, theOutputFormat)
  OB_STATIC_PLUGIN(PCModelFormat, thePCModelFormat)
  OB_STATIC_PLUGIN(PDBFormat, thePDBFormat)
  OB_STATIC_PLUGIN(PDBQTFormat, thePDBQTFormat)
#ifdef HAVE_LIBZ
  OB_STATIC_PLUGIN(PNGFormat, thePNGFormat)
#endif
  OB_STATIC_PLUGIN(PointCloudFormat, thePointCloudFormat)
  OB_STATIC_PLUGIN(PovrayFormat, thePovrayFormat)
  OB_STATIC_PLUGIN(PQRFormat, thePQRFormat)
  OB_STATIC_PLUGIN(PQSFormat, thePQSFormat)
  OB_STATIC_PLUGIN(PWscfFormat, thePWscfFormat)
  OB_STATIC_PLUGIN(QChemOutputFormat, theQChemOutputFormat)
  OB_STATIC_PLUGIN(QChemInputFormat, theQChemInputFormat)
  OB_STATIC_PLUGIN(ReportFormat, theReportFormat)
  OB_STATIC_PLUGIN(SmiReactFormat, theSmiReactFormat)
  OB_STATIC_PLUGIN(RXNFormat, theRXNFormat)
  OB_STATIC_PLUGIN(ShelXFormat, theShelXFormat)
  OB_STATIC_PLUGIN(SMIFormat, theSMIFormat)
  OB_STATIC_PLUGIN(STLFormat, theSTLFormat)
  OB_STATIC_PLUGIN(CANSMIFormat, theCANSMIFormat)
  OB_STATIC_PLUGIN(FIXFormat, theFIXFormat)
  OB_STATIC_PLUGIN(SVGFormat, theSVGFormat)
  OB_STATIC_PLUGIN(TextFormat, theTextFormat)
  OB_STATIC_PLUGIN(ThermoFormat, theThermoFormat)
  OB_STATIC_PLUGIN(TinkerFormat, theTinkerFormat)
  OB_STATIC_PLUGIN(TitleFormat, theTitleFormat)
  OB_STATIC_PLUGIN(TurbomoleFormat, theTurbomoleFormat)
  OB_STATIC_PLUGIN(UniChemFormat, theUniChemFormat)
  OB_STATIC_PLUGIN(VASPFormat, theVASPFormat)
  OB_STATIC_PLUGIN(ViewMolFormat, theViewMolFormat)
  OB_STATIC_PLUGIN(XEDFormat, theXEDFormat)
  OB_STATIC_PLUGIN(XSFFormat, theXSFFormat)
  OB_STATIC_PLUGIN(XYZFormat, theXYZFormat)
  OB_STATIC_PLUGIN(YOBFormat, theYOBFormat)
  OB_STATIC_PLUGIN(ZINDOFormat, theZINDOFormat)
#ifdef HAVE_STATIC_LIBXML
  OB_STATIC_PLUGIN(ChemDrawXMLFormat, theChemDrawXMLFormat)
  OB_STATIC_PLUGIN(CMLFormat, theCMLFormat)
  OB_STATIC_PLUGIN(CMLReactFormat, theCMLReactFormat)
  OB_STATIC_PLUGIN(PubChemFormat, thePubChemFormat)
  OB_STATIC_PLUGIN(XMLFormat, theXMLFormat)
#endif
#ifdef HAVE_STATIC_INCHI
  OB_STATIC_PLUGIN(InChIFormat, theInChIFormat)
#endif
#ifdef HAVE_REGEX_H
  OB_STATIC_PLUGIN(GAMESSUKInputFormat, theGAMESSUKInputFormat)
  OB_STATIC_PLUGIN(GAMESSUKOutputFormat, theGAMESSUKOutputFormat)
#endif
#ifdef HAVE_RPC_XDR_H
  OB_STATIC_PLUGIN(XTCFormat, theXTCFormat)
#endif

  // descriptors
  OB_STATIC_PLUGIN(CanSmiles, theCanSmiles)
  OB_STATIC_PLUGIN(CompoundFilter, dummyCmpFilter)
  OB_STATIC_PLUGIN(MWFilter, theMWFilter)
  OB_STATIC_PLUGIN(SmartsFilter, firstSmartsFilter)
  OB_STATIC_PLUGIN(SmartsFilter, secondSmartsFilter)
  OB_STATIC_PLUGIN(TitleFilter, theTitleFilter)
  OB_STATIC_PLUGIN(FormulaDescriptor, TheFormulaDescriptor)
  //OB_STATIC_PLUGIN(FPCount, theFPCount)
  OB_STATIC_PLUGIN(InChIFilter, theInChIFilter)
  // smarts descriptors
  OB_STATIC_PLUGIN(SmartsDescriptor, theHBD)
  OB_STATIC_PLUGIN(SmartsDescriptor, theHBA1)
  OB_STATIC_PLUGIN(SmartsDescriptor, theHBA2)
  OB_STATIC_PLUGIN(SmartsDescriptor, thenF)
  // group contribution descriptors
  OB_STATIC_PLUGIN(OBGroupContrib, thelogP)
  OB_STATIC_PLUGIN(OBGroupContrib, theTPSA)
  OB_STATIC_PLUGIN(OBGroupContrib, theMR)

  // fingerprints
  OB_STATIC_PLUGIN(fingerprint2, thefingerprint2)
  OB_STATIC_PLUGIN(PatternFP, FP3PatternFP)
  OB_STATIC_PLUGIN(PatternFP, FP4PatternFP)
  OB_STATIC_PLUGIN(fingerprintECFP, theECFP0)
  OB_STATIC_PLUGIN(fingerprintECFP, theECFP2)
  OB_STATIC_PLUGIN(fingerprintECFP, theECFP4)
  OB_STATIC_PLUGIN(fingerprintECFP, theECFP6)
  OB_STATIC_PLUGIN(fingerprintECFP, theECFP8)
  OB_STATIC_PLUGIN(fingerprintECFP, theECFP10)

  // forcefields
  OB_STATIC_PLUGIN(OBForceFieldGaff, theForceFieldGaff)
  OB_STATIC_PLUGIN(OBForceFieldGhemical, theForceFieldGhemical)
  OB_STATIC_PLUGIN(OBForceFieldMMFF94, theForceFieldMMFF94)
  OB_STATIC_PLUGIN(OBForceFieldMMFF94, theForceFieldMMFF94s)
  OB_STATIC_PLUGIN(OBForceFieldUFF, theForceFieldUFF)

  // operations
  OB_STATIC_PLUGIN(OpAddInIndex, theOpAddInIndex)
  OB_STATIC_PLUGIN(OpAddPolarH, theOpAddPolarH)
  OB_STATIC_PLUGIN(OpAddNonPolarH, theOpAddNonPolarH)
  OB_STATIC_PLUGIN(OpChangeCell, theOpChangeCell)
  OB_STATIC_PLUGIN(OpCanonical, theOpCanonical)
  OB_STATIC_PLUGIN(OpDelPolarH, theOpDelPolarH)
  OB_STATIC_PLUGIN(OpDelNonPolarH, theOpDelNonPolarH)
  OB_STATIC_PLUGIN(OpFillUC, theOpFillUC)
  OB_STATIC_PLUGIN(OpEnergy, theOpEnergy)
  OB_STATIC_PLUGIN(OpMinimize, theOpMinimize)
  OB_STATIC_PLUGIN(OpGen2D, theOpGen2D)
  OB_STATIC_PLUGIN(OpGen3D, theOpGen3D)
  OB_STATIC_PLUGIN(OpNewS, theOpNewS)
  OB_STATIC_PLUGIN(OpPartialCharge, theOpPartialCharge)
  OB_STATIC_PLUGIN(OpReadConformers, theOpReadConformers)
  OB_STATIC_PLUGIN(OpSort, theOpSort)
  OB_STATIC_PLUGIN(OpExtraOut, theOpExtraOut)
#ifdef HAVE_STATIC_INCHI
  OB_STATIC_PLUGIN(OpUnique, theOpUnique)
#endif
#ifdef HAVE_EIGEN
  OB_STATIC_PLUGIN(OpConformer, theOpConformer)
#endif

  // charges
  OB_STATIC_PLUGIN(GasteigerCharges, theGasteigerCharges)
  OB_STATIC_PLUGIN(MMFF94Charges, theMMFF94Charges)
  OB_STATIC_PLUGIN(NoCharges, theNoCharges)
  OB_STATIC_PLUGIN(FromFileCharges, theFromFileCharges)
#ifdef HAVE_EIGEN
  OB_STATIC_PLUGIN(QEqCharges, theQEqCharges)
  OB_STATIC_PLUGIN(QTPIECharges, theQTPIECharges)
#endif
#ifdef HAVE_EIGEN3
  OB_STATIC_PLUGIN(EQEqCharges, theEQEqCharges)
#endif
  OBAPI std::vector<std::string> EnableStaticPlugins();

#endif // USING_DYNAMIC_LIBS

#endif // SWIG

} // end namespce

#endif
