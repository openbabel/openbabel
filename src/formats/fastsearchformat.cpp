/**********************************************************************
fastsearchformat.cpp: Preparation and searching of fingerprint-based index files
Copyright (C) 2005-2006 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/fingerprint.h>

using namespace std;
namespace OpenBabel {

  /// \brief Prepares and searches of fingerprint-based index files. See FastSearch class for details
  class FastSearchFormat : public OBFormat
  {
  public:
    //Register this format type ID
    FastSearchFormat() : fsi(NULL) 
    {
      OBConversion::RegisterFormat("fs",this);
      //Specify the number of option taken by options
      OBConversion::RegisterOptionParam("S", this, 1, OBConversion::GENOPTIONS);
      OBConversion::RegisterOptionParam("S", this, 1, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("f", this, 1);
      OBConversion::RegisterOptionParam("N", this, 1);
      OBConversion::RegisterOptionParam("u", this, 0);
      OBConversion::RegisterOptionParam("t", this, 1, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("l", this, 1, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("a", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("e", this, 0, OBConversion::INOPTIONS);
    }
	
    virtual const char* Description() //required
    { return
        "FastSearching\n"
        "Uses molecular fingerprints in an index file.\n"
"Writing to the fs format makes an index (a very slow process)\n"
"  babel datafile.xxx index.fs\n"
"Reading from the fs format does a fast search for:\n"
"  Identical molecule\n"
"    babel index.fs -sSMILES outfile.yyy -ae  or\n"
"    babel datafile.xxx -ifs -sSMILES outfile.yyy -ae\n"
"  Substructure\n"
"    babel index.fs -sSMILES outfile.yyy   or\n"
"    babel datafile.xxx -ifs -sSMILES outfile.yyy\n"
"  Molecular similarity based on Tanimoto coefficient\n"
"    babel index.fs -sSMILES outfile.yyy -t0.7  (Tanimoto >0.7)\n"
"    babel index.fs -sSMILES outfile.yyy -t15   (best 15 molecules)\n"
"  The structure spec can be a molecule from a file: -Spatternfile.zzz\n\n"
"Note that the parameter of the -s option needs to be a valid SMILES\n"
"molecule when using fastsearch. You can use the more versatile SMARTS\n"
"in a normal substructure search.\n\n"

"Write Options (when making index) e.g. -xfFP3\n"
" f# Fingerprint type\n"
" N# Fold fingerprint to # bits\n"
" u  Update an existing index\n\n"
"Read Options (when searching) e.g. -at0.7\n"
" t# Do similarity search: #mols or # as min Tanimoto\n"
" a  Add Tanimoto coeff to title in similarity search\n"
" l# Maximum number of candidates. Default<4000>\n"
" e  Exact match\n"
" S\"filename\"  Structure spec in a file:\n"
" h  SMARTS uses explicit H in pattern file\n\n"
;
    };

    virtual unsigned int Flags(){return READBINARY | READONEONLY | WRITEBINARY;};

  public:
    virtual bool ReadChemObject(OBConversion* pConv);
    virtual bool WriteChemObject(OBConversion* pConv);

  private:
      bool ObtainTarget(OBConversion* pConv, OBMol& patternMol, const string& indexname);

  private:
    ///big data structure which will remain in memory after it is loaded
    //until the program ends.
    FastSearch fs;
    FastSearchIndexer* fsi;
    streampos LastSeekpos; //used during update
    OBStopwatch sw; //used when preparing index
  };

  ///////////////////////////////////////////////////////////////
  //Make an instance of the format class
  FastSearchFormat theFastSearchFormat;

  ///////////////////////////////////////////////////////////////
  bool FastSearchFormat::ReadChemObject(OBConversion* pConv)
  {
    //Searches index file for structural matches
    //This function is called only once per search

    std::string auditMsg = "OpenBabel::Read fastsearch index ";
    std::string description(Description());
    auditMsg += description.substr(0,description.find('\n'));
    obErrorLog.ThrowError(__FUNCTION__,
                          auditMsg,
                          obAuditMsg);

    //Derive index name
    string indexname = pConv->GetInFilename();
    string::size_type pos=indexname.find_last_of('.');
    if(pos!=string::npos)
      {
        indexname.erase(pos);
        indexname += ".fs";
      }

    //Have to open input stream again because needs to be in binary mode
    ifstream ifs;
    stringstream errorMsg;
    if(!indexname.empty())
      ifs.open(indexname.c_str(),ios::binary);
    if(!ifs)
      {
        errorMsg << "Couldn't open " << indexname << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }

    string datafilename = fs.ReadIndex(&ifs);
    if(datafilename.empty())
      {
        errorMsg << "Difficulty reading from index " << indexname << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
    
    OBMol patternMol;
    bool doSubset = pConv->IsOption("s",OBConversion::INOPTIONS)!=NULL;// -as option
    bool exactmatch = pConv->IsOption("e",OBConversion::INOPTIONS)!=NULL;// -ae option
    if(!doSubset)
    {
      //Similarity or substructure
      if(!ObtainTarget(pConv, patternMol, indexname))
        return false;
    }
    
    //Open the datafile and put it in pConv
    //datafile name derived from index file probably won't have a file path
    //but indexname may. Derive a full datafile name
    string path;
    pos = indexname.find_last_of("/\\");
    if(pos==string::npos)
      path = datafilename;
    else
      path = indexname.substr(0,pos+1) + datafilename;
	
    ifstream datastream(path.c_str());
    if(!datastream)
      {
        errorMsg << "Difficulty opening " << path << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
    pConv->SetInStream(&datastream);
	
    //Input format is currently fs; set it appropriately
    if(!pConv->SetInAndOutFormats(pConv->FormatFromExt(datafilename.c_str()),pConv->GetOutFormat()))
			return false;
    pConv->AddOption("b",OBConversion::GENOPTIONS);


    //Now do searching
    const char* p = pConv->IsOption("t",OBConversion::INOPTIONS);
    if(p)
      {
        //Do a similarity search
        multimap<double, unsigned int> SeekposMap;
        string txt=p;
        if(txt.find('.')==string::npos)
          {
            //Finds n molecules with largest Tanimoto
            int n = atoi(p);
            fs.FindSimilar(&patternMol, SeekposMap, n);
          }
        else
          {
            //Finds molecules with Tanimoto > MinTani
            double MinTani = atof(txt.c_str());
//            if(doSubset)
//              fs.FindSubset(SeekposMap, MinTani);
//            else
            fs.FindSimilar(&patternMol, SeekposMap, MinTani);
          }
		
        //Don't want to filter through SMARTS filter
        pConv->RemoveOption("s", OBConversion::GENOPTIONS);
		
        multimap<double, unsigned int>::reverse_iterator itr;
        for(itr=SeekposMap.rbegin();itr!=SeekposMap.rend();++itr)
          {
            datastream.seekg(itr->second);

            if(pConv->IsOption("a", OBConversion::INOPTIONS))
              {
                //Adds Tanimoto coeff to title
                //First remove any previous value
                pConv->RemoveOption("addtotitle", OBConversion::GENOPTIONS);
                stringstream ss;
                ss << " " << itr->first;
                pConv->AddOption("addtotitle",OBConversion::GENOPTIONS, ss.str().c_str());
			
              }
            pConv->SetOneObjectOnly();
            if(itr != --SeekposMap.rend())
              pConv->SetMoreFilesToCome();//so that not seen as last on output 
            pConv->Convert(NULL,NULL);
          }
      }

    else
      {
        //Structure search
        int MaxCandidates = 4000;
        p = pConv->IsOption("l",OBConversion::INOPTIONS);
        if(p && atoi(p))
          MaxCandidates = atoi(p);
		
        vector<unsigned int> SeekPositions;
        if(exactmatch)
        {
          //Find mols where all fingerprint bits are the same as the target
          fs.FindMatch(&patternMol, SeekPositions, MaxCandidates);
          // ensure that SMARTS filter in transform.cpp looks only for an exact match
          // by setting an option with the number of heavy atoms in the pattern mol included.
          stringstream ss;
          ss << patternMol.NumHvyAtoms();
          pConv->AddOption("exactmatch", OBConversion::GENOPTIONS, ss.str().c_str());
        }

        else
        {
          //Do a substructure search
          fs.Find(&patternMol, SeekPositions, MaxCandidates);
          clog << SeekPositions.size() << " candidates from fingerprint search phase" << endl;
        }

          //Output the candidate molecules 
          //filtering through s filter, unless the fingerprint type does not require it
          if(fs.GetFingerprint()->Flags() & OBFingerprint::FPT_UNIQUEBITS)
            pConv->RemoveOption("s",OBConversion::GENOPTIONS);

          vector<unsigned int>::iterator itr;
          for(itr=SeekPositions.begin();itr!=SeekPositions.end();itr++)
            {
              datastream.seekg(*itr);
              //	datastream.seekg(*itr - datastream.tellg(), ios_base::cur); //Avoid retrieving start

              //debugging kludge to output all candidates directly
              if(pConv->IsOption("c",OBConversion::GENOPTIONS))
                {
                  string ln;
                  getline(datastream,ln);
                  datastream.seekg(*itr);
                  *pConv->GetOutStream() << "** " << ln << endl;
                }
              pConv->SetOneObjectOnly();
              if(itr+1 != SeekPositions.end())
                pConv->SetMoreFilesToCome();//so that not seen as last on output 
              pConv->Convert(NULL,NULL);
            }
      }
    return false;	//To finish	
  }

  /////////////////////////////////////////////////////
  bool FastSearchFormat::WriteChemObject(OBConversion* pConv)
  {
    //Prepares or updates an index file. Called for each molecule indexed

    bool update = pConv->IsOption("u")!=NULL;

    ostream* pOs = pConv->GetOutStream();// with named index it is already open
    bool NewOstreamUsed=false;
    if(fsi==NULL)
      {
        //First pass sets up FastSearchIndexer object
        string mes("prepare an");
        if(update)
          mes = "update the";
        clog << "This will " << mes << " index of " << pConv->GetInFilename()
             <<  " and may take some time..." << flush;

        if(!pConv->IsLastFile())
        {
          obErrorLog.ThrowError(__FUNCTION__,
            "There should not be multiple input files. A .fs file is an index of a single datafile.",
            obError);
          return false;
        }

        std::string auditMsg = "OpenBabel::Write fastsearch index ";
        std::string description(Description());
        auditMsg += description.substr( 0, description.find('\n') );
        obErrorLog.ThrowError(__FUNCTION__,auditMsg,obAuditMsg);

        sw.Start();
		
        FptIndex* pidx; //used with update

        if(pOs==&cout)
          {
            //No index filename specified
            //Derive index name from datafile name
            string indexname=pConv->GetInFilename();
            string::size_type pos=indexname.find_last_of('.');
            if(pos!=string::npos)
              indexname.erase(pos);
            indexname += ".fs";

            bool idxok=true;
            if(update)
              {
                LastSeekpos = 0;

                //Read in existing index
                idxok=false;
                ifstream ifs(indexname.c_str(),ifstream::binary);
                if(ifs.good())
                  {
                    pidx = new FptIndex;
                    idxok = pidx->Read(&ifs);
                  }
              }//ifs closed here

            pOs = new ofstream(indexname.c_str(),ofstream::binary);

            if(!pOs->good() || !idxok)
              {
                stringstream errorMsg;
                errorMsg << "Trouble opening or reading " << indexname << endl;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
                return false;
              }
            NewOstreamUsed=true;
          }
        else // not cout
          {
            if(update)
              {
                obErrorLog.ThrowError(__FUNCTION__,
                  "Currently, updating	can only be done with index files that \
have the same name as the datafile. Use the form:\n \
	babel datafile.xxx -ofs -xu", obError);
                return false;
              }
          }

        int nbits = 0;
        const char* p = pConv->IsOption("N");
        if(p)
          nbits = atoi(p);

        string fpid; //fingerprint type
        p=pConv->IsOption("f");
        if(p)
          fpid=p;

        //Prepare name without path
        string datafilename = pConv->GetInFilename();
        if(datafilename.empty())
          {
            obErrorLog.ThrowError(__FUNCTION__, "No datafile!", obError);
            return false;
          }
        string::size_type pos = datafilename.find_last_of("/\\");
        if(pos!=string::npos)
          datafilename=datafilename.substr(pos+1);

        if(update)
          {
            fsi = new FastSearchIndexer(pidx,pOs);//using existing index

            //Seek to position in datafile of last of old objects
            LastSeekpos = *(pidx->seekdata.end()-1);
            pConv->GetInStream()->seekg(LastSeekpos);
          }
        else
          fsi = new FastSearchIndexer(datafilename, pOs, fpid, nbits);
		
        obErrorLog.StopLogging();
      }

    //All passes provide an object for indexing
    OBBase* pOb = pConv->GetChemObject();
    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    if(pmol)
      pmol->ConvertDativeBonds();//use standard form for dative bonds
	
    streampos seekpos = pConv->GetInPos();
    if(!update || seekpos>LastSeekpos) 
      fsi->Add(pOb, seekpos );
    else
      //Don't index old objects during update. Don't increment pConv->Index.
      pConv->SetOutputIndex(pConv->GetOutputIndex()-1);

    if(pConv->IsLast())
      {
        //Last pass 
        delete fsi; //saves index file
        if(NewOstreamUsed)
          delete pOs;

        //return to starting conditions
        fsi=NULL;

        obErrorLog.StartLogging();

        double secs = sw.Elapsed();
        if(secs>150)
          clog << "\n It took " << secs/60 << " minutes" << endl;
        else
          clog << "\n It took " << secs << " seconds" << endl;
      }
    delete pOb;
    return true;
  }

///////////////////////////////////////////////////////////////
  bool FastSearchFormat::ObtainTarget(OBConversion* pConv, OBMol& patternMol, const string& indexname)
  {
    //Obtains an OBMol
    //   either from the SMARTS string in the -s option
    //   or by converting the file in the -S option
    //or, if neither option is provided, displays information on the index file.

    stringstream smiles(stringstream::out);
    ifstream patternstream;
    OBConversion PatternConv(&patternstream,&smiles);

    const char* p = pConv->IsOption("s",OBConversion::GENOPTIONS);
    string txt;
    if(p) 
      {
        // Use the -s option
        txt=p;
        stringstream smarts(txt, stringstream::in);		
        OBConversion Convsm(&smarts);
        if(!Convsm.SetInFormat("smi")) return false;
        Convsm.Read(&patternMol);

        if(patternMol.Empty())
        {
          obErrorLog.ThrowError(__FUNCTION__, 
            "Could not make a molecule from " + smarts.str()
            + "\nThis needs to be valid SMILES when using fastsearch."
            "You can use the more versatile SMARTS in a normal substructure search." , obError);
            return false;
        }
      }
    else
      {
      // or Make OBMol from file in -S option or -aS option	
      p = pConv->IsOption("S",OBConversion::GENOPTIONS);
      if(!p)
        p = pConv->IsOption("S",OBConversion::INOPTIONS);//for GUI mainly
      }

    if(!p)
      {
        //neither -s or -S options provided. Output info rather than doing search
        const FptIndexHeader& header = fs.GetIndexHeader();
        string id(header.fpid);
        if(id.empty())
          id = "default";
        clog << indexname << " is an index of\n " << header.datafilename 
             << ".\n It contains " << header.nEntries 
             << " molecules. The fingerprint type is " << id << " with "
             << OBFingerprint::Getbitsperint() * header.words << " bits.\n"
             << "Typical usage for a substructure search:\n"
             << "babel indexfile.fs -osmi -sSMILES" << endl;
        return false;
      }

    if(p && patternMol.Empty())
      {
        txt=p;
        string::size_type pos = txt.find_last_of('.');
        if(pos==string::npos)
          {
            obErrorLog.ThrowError(__FUNCTION__, "Filename of pattern molecule in -S option must\n"
              "have an extension to define its format.", obError);
            return false;
          }
        patternstream.open(txt.c_str());
        if(!patternstream)
          {
            stringstream errorMsg;
		  
            errorMsg << "Cannot open " << txt << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
            return false;
          }

        PatternConv.SetOneObjectOnly();
        if(PatternConv.SetInFormat(txt.substr(pos+1).c_str()))
          PatternConv.Read(&patternMol);
      }

    if(patternMol.Empty())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Cannot derive a molecule from the -s or -S options", obWarning);
        return false;
      }
    patternMol.ConvertDativeBonds();//use standard form for dative bonds

    //If the -s option is not already present, generate one by converting to SMILES
    if(!pConv->IsOption("s",OBConversion::GENOPTIONS))
    {
      if(!PatternConv.SetOutFormat("smi"))
        return false;
      //if -ah option ensure explicit H remains as such in SMARTS phase
      if(pConv->IsOption("h",OBConversion::INOPTIONS))
        PatternConv.AddOption("h",OBConversion::OUTOPTIONS);
      PatternConv.Write(&patternMol);
      //remove name to leave smiles string
      string smilesstr(smiles.str());
      string::size_type pos = smilesstr.find_first_of(" \t\r\n");
      if(pos!=string::npos)
        smilesstr = smilesstr.substr(0,pos);
      pConv->AddOption("s", OBConversion::GENOPTIONS, smilesstr.c_str());
    }

    return true;
  }

}//Openbabel

//! \file fastsearchformat.cpp
//! \brief Preparation and searching of fingerprint-based index files
