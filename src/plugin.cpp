/**********************************************************************
plugin.cpp - facilitates construction of plugin classes

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

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>
#include <openbabel/oberror.h>

#include <iterator>

using namespace std;
namespace OpenBabel
{

OBPlugin::PluginMapType& OBPlugin::GetTypeMap(const char* PluginID)
{
  PluginMapType::iterator itr;

  // Make sure the plugins are loaded
  if (AllPluginsLoaded == 0) {
    OBPlugin::LoadAllPlugins();
  }

  itr = PluginMap().find(PluginID);
  if(itr!=PluginMap().end())
    return itr->second->GetMap();
  return PluginMap();//error: type not found; return plugins map
}

int OBPlugin::AllPluginsLoaded = 0;

void OBPlugin::LoadAllPlugins()
{
  int count = 0;
#if  defined(USING_DYNAMIC_LIBS)
  // Depending on availability, look successively in
  // FORMATFILE_DIR, executable directory or current directory
  string TargetDir;

#ifdef FORMATFILE_DIR
  TargetDir="FORMATFILE_DIR";
#endif

  DLHandler::getConvDirectory(TargetDir);

  vector<string> files;
  if(!DLHandler::findFiles(files,DLHandler::getFormatFilePattern(),TargetDir)) {
    obErrorLog.ThrowError(__FUNCTION__, "Unable to find OpenBabel plugins. Try setting the BABEL_LIBDIR environment variable.", obError);
    return;
  }

  vector<string>::iterator itr;
  for(itr=files.begin();itr!=files.end();++itr) {
    if(DLHandler::openLib(*itr))
      count++;
  }
  if(!count) {
    string error = "No valid OpenBabel plugs found in "+TargetDir;
    obErrorLog.ThrowError(__FUNCTION__, error, obError);
    return;
  }
#else
  count = 1; // Avoid calling this function several times
#endif //USING_DYNAMIC_LIBS

  // Status have to be updated now
  AllPluginsLoaded = count;

  // Make instances for plugin classes defined in the data file.
  // This is hook for OBDefine, but does nothing if it is not loaded
  // or if plugindefines.txt is not found.
  OBPlugin* pdef = OBPlugin::GetPlugin("loaders","define");
  if(pdef) {
    static vector<string> vec(3);
    vec[1] = string("define");
    vec[2] = string("plugindefines.txt");
    pdef->MakeInstance(vec);
  }

  return;
}

OBPlugin* OBPlugin::BaseFindType(PluginMapType& Map, const char* ID)
{
  // Make sure the plugins are loaded
  if (AllPluginsLoaded == 0) {
    OBPlugin::LoadAllPlugins();
  }

  if(!ID || !*ID)
    return NULL;
  PluginMapType::iterator itr = Map.find(ID);
  if(itr==Map.end())
    return NULL;
  else
    return itr->second;
}

OBPlugin* OBPlugin::GetPlugin(const char* Type, const char* ID)
{
  if(Type!=NULL)
    return BaseFindType(GetTypeMap(Type), ID);

  // Make sure the plugins are loaded
  if (AllPluginsLoaded == 0) {
    OBPlugin::LoadAllPlugins();
  }

  //When Type==NULL, search all types for matching ID and stop when found
  PluginMapType::iterator itr;
  for(itr=PluginMap().begin();itr!= PluginMap().end();++itr)
  {
    OBPlugin* result = BaseFindType(itr->second->GetMap(), ID);
    if(result)
      return result;
  }
  return NULL; //not found
}

bool OBPlugin::ListAsVector(const char* PluginID, const char* param, vector<string>& vlist)
{
  PluginMapType::iterator itr;
  bool ret=true;

  // Make sure the plugins are loaded
  if (AllPluginsLoaded == 0) {
    LoadAllPlugins();
  }

  if(PluginID)
  {
    if(*PluginID!=0 && strcmp(PluginID, "plugins"))
    {
      //List the sub classes of the specified type
      itr = PluginMap().find(PluginID);
      if(itr!=PluginMap().end())
      {
        bool onlyIDs = param!=NULL && strstr(param,"ids")!=NULL;
        //Get map of plugin type (like OBFingerprint) and output its contents
        PluginMapType Map = itr->second->GetMap();
        for(itr=Map.begin(); itr!=Map.end(); ++itr)
        {
          if(*(itr->first)=='_')//no listing when ID starts with '_'
            continue;
          if(onlyIDs)
            vlist.push_back(itr->first);
          else
          {
            string txt;
            if((itr->second)->Display(txt, param, itr->first))
              vlist.push_back(txt);
          }
        }
        return true;
      }
      ret=false; //asked for a type not available; provide plugin types instead
    }
  }
  //List the plugin types
  for(itr=PluginMap().begin();itr!= PluginMap().end();++itr)
    vlist.push_back(itr->first);
  return ret;
}

void OBPlugin::List(const char* PluginID, const char* param, ostream* os)
{
  vector<string> vlist;
  if(!ListAsVector(PluginID,param, vlist))
    *os << PluginID << " is not a recognized plugin type. Those with instances of sub-types loaded are:" << endl;
  copy(vlist.begin(), vlist.end(), std::ostream_iterator<string>(*os, "\n"));
}

string OBPlugin::ListAsString(const char* PluginID, const char* param)
{
  stringstream ss;
  List(PluginID, param, &ss);
  return ss.str();
}

string OBPlugin::FirstLine(const char* txt)
{
  string stxt(txt);
  string::size_type pos = stxt.find('\n');
  if(pos==string::npos)
    return stxt;
  else
    return stxt.substr(0,pos);
}

//Default version
bool OBPlugin::Display(string& txt, const char* param, const char* ID)
{
  //Use the provided ID if possible.
  if(ID)
    txt = ID;
  else
    txt = GetID();
  txt += "    ";// was '\t'; but caused problems in GUI menu
  if(param && !strcasecmp(param, "verbose"))
  {
    txt += Description();
    txt += '\n';
  }
  else
    txt += FirstLine(Description());
  return true;
}

#ifndef USING_DYNAMIC_LIBS

std::vector<std::string> EnableStaticPlugins()
{
  // Return list of all plugin ids. This also ensures the code is not removed
  // by compiler optimization.
  std::vector<std::string> plugin_ids;
  // formats
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theABINITFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theAcesOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theAcesInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theACRFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theADFOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theADFInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&t41Format__)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theAlchemyFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theAmberPrepFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theAoforceFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOBAPIInterface)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theBallStickFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theBGFFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theBoxFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCacaoFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCacheFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCARFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCASTEPFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCCCFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCHEM3D1Format)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCHEM3D2Format)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theChemDrawBinaryXFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theChemDrawFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theChemKinFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCHTFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCIFFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCopyFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCRK2DFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCRK3DFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCSRFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCSSRFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theDlpolyConfigFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theDlpolyHISTORYFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theDMolFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theEXYZFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFASTAFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFastSearchFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFCHKFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFEATFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFenskeZmatFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFHIaimsFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFingerprintFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFreeFormFractionalFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGAMESSOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGAMESSInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGaussianCubeFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGaussianOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGaussianInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGaussianZMatrixInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGenBankFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGhemicalFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGROFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGROMOS96Format)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGULPFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theHINFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theJaguarOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theJaguarInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theLMPDATFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMCDLFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMOLFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theSDFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&themmCIFFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMacroModFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMNAFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMOL2Format)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&moldenFormat__)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMolproOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMolproInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMolReportFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMOPACFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMOPACCARTFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMOPACINTFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMPDFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMPQCFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMPQCInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMSIFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&msmsFormat__)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theNulFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theNWChemOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theNWChemInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpenDXCubeFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePCModelFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePDBFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePDBQTFormat)->GetID());
#ifdef HAVE_LIBZ
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePNGFormat)->GetID());
#endif
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePointCloudFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePovrayFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePQRFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePQSFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePWscfFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theQChemOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theQChemInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theReportFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theSmiReactFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theRXNFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theShelXFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theSMIFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theSTLFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCANSMIFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFIXFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theSVGFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theTextFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theThermoFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theTinkerFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theTitleFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theTurbomoleFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theUniChemFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theVASPFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theViewMolFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theXEDFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theXSFFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theXYZFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theYOBFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theZINDOFormat)->GetID());
#ifdef HAVE_STATIC_LIBXML
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theChemDrawXMLFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCMLFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCMLReactFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thePubChemFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theXMLFormat)->GetID());
#endif
#ifdef HAVE_STATIC_INCHI
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theInChIFormat)->GetID());
#endif
#ifdef HAVE_RPC_XDR_H
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theXTCFormat)->GetID());
#endif
#ifdef HAVE_REGEX_H
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGAMESSUKInputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGAMESSUKOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOrcaOutputFormat)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOrcaInputFormat)->GetID());
#endif

  // descriptors
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theCanSmiles)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&dummyCmpFilter)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMWFilter)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&firstSmartsFilter)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&secondSmartsFilter)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theTitleFilter)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&TheFormulaDescriptor)->GetID());
  //plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFPCount)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theInChIFilter)->GetID());
  // smarts descriptors
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theHBD)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theHBA1)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theHBA2)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thenF)->GetID());
  // group contribution descriptors
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thelogP)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theTPSA)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMR)->GetID());

  // fingerprints
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&thefingerprint2)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&FP3PatternFP)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&FP4PatternFP)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theECFP0)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theECFP2)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theECFP4)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theECFP6)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theECFP8)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theECFP10)->GetID());

  // forcefields
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theForceFieldGaff)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theForceFieldGhemical)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theForceFieldMMFF94)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theForceFieldMMFF94s)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theForceFieldUFF)->GetID());

  // operations
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpAddInIndex)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpAddPolarH)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpAddNonPolarH)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpCanonical)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpDelPolarH)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpDelNonPolarH)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpFillUC)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpEnergy)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpMinimize)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpGen2D)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpGen3D)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpNewS)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpPartialCharge)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpReadConformers)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpSort)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpExtraOut)->GetID());
#ifdef HAVE_EIGEN
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpConformer)->GetID());
#endif
#ifdef HAVE_STATIC_INCHI
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theOpUnique)->GetID());
#endif

  // charges
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theGasteigerCharges)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theMMFF94Charges)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theNoCharges)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theFromFileCharges)->GetID());
#ifdef HAVE_EIGEN
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theQEqCharges)->GetID());
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theQTPIECharges)->GetID());
#endif
#ifdef HAVE_EIGEN3
  plugin_ids.push_back(reinterpret_cast<OBPlugin*>(&theEQEqCharges)->GetID());
#endif

  return plugin_ids;
}

#endif // USING_DYNAMIC_LIBS

/**
  /// @addtogroup plugins Plugins
  ///@{
 @page plugins Plugins
 Plugins are a way of extending OpenBabel without having to alter any of the
 existing code. They may be built as shared libraries (DLLs with an extension
 .obf or .so files in a specified location) and distributed separately, but
 plugin classes can also be in the main code. In both cases they are discovered
 at startup when a global instance of the plugin class is instantiated. It iss
 registered by its constructor and is added to a static record of all the
 plugins of its particular type that are currently loaded.

 There are two levels of plugin. The top layer (at the time of writing) are:
  formats descriptors fingerprints forcefields charges ops loaders
 but additional types can be added without disturbing the main API. At runtime
   obabel -L
 will list the top level of plugins. They typically are abstract classes with
 virtual functions that define an interface for that type. Classes derived
 from these are the second layer of plugins, and can be listed at runtime like,
 for instance:
   obabel -L formats cml
 where formats is the top level of plugin and cml is the id of a derived class
 of this type.

 The top level of plugins will usually have their interfaces declared in header
 files compiled with the main API. The second level of plugin will typically
 not be known to the API at compile time, usually will not have a header file
 and must be accessed indirectly, to allow for the possibility that they may
 not be loaded:
 @code
   OBOp* pOp = OBOp::FindType("gen3D");
   if(!pOp)
     ...report error
   pOp->Do(mol);
 @endcode
 This retrieves the global instance of the plugin. This is usually adequate but
 making a new instance may be appropriate in some cases.

 Instances of some plugin classes can be constructed at startup from information
 in a text file and used in the same way as those defined in code. See OBDefine.
 This is appropriate for some classes that differ only by the datafile or
 SMARTS strings they use.
*/
///@}
}//namespace

//! \file plugin.cpp
//! \brief Simplify 'plugin' classes to be discovered and/or loaded at runtime
