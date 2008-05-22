/**********************************************************************
Copyright (C) 2007 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <openbabel/babelconfig.h>
#include <openbabel/format.h>
#include <openbabel/obconversion.h>
#include <zlib.h>

using namespace std;
namespace OpenBabel
{

class PNGFormat : public OBFormat
{
public:
  PNGFormat()
  {
    OBConversion::RegisterFormat("png",this);
    OBConversion::RegisterOptionParam("y", this, 1, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("y", this, 1, OBConversion::OUTOPTIONS);
  }

  virtual const char* Description()
  {
    return
    "PNG files with embedded data\n"
    "Extract chemical structure data embedded in PNG image files\n"
    "Read Options e.g. -ay\n"
    " y <additional chunk ID> Look also in chunks with specified ID\n\n"

    "Write Options e.g. -xO cml\n"
    " O <format ID> Format of embedded text\n"
    " y <additional chunk ID> Write to a chunk with specified ID\n\n";
  };

  virtual const char* TargetClassDescription()
  {
    static string txt;
    txt = " PNG_files\n"; //so reports "n PNG_files converted"
    txt += OBFormat::TargetClassDescription(); //to display OBMol options in GUI
    return txt.c_str();
  }

  virtual unsigned int Flags()
  {
      return READONEONLY | READBINARY | WRITEBINARY;
  };

  virtual bool ReadChemObject(OBConversion* pConv)
  {
    bool ret = ReadMolecule(NULL, pConv);
    pConv->GetChemObject(); //increments output index
    return ret;
  };
  virtual bool WriteChemObject(OBConversion* pConv)
  {
    OBBase* pOb = pConv->GetChemObject();
    bool ret = WriteMolecule(pOb, pConv);
    return ret;
  };

virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:
  int count; //number of chemical objects converted
  vector<char> CopyOfInput;
  unsigned bytesToIEND; //number of bytes upto but not including the IEND chunk.
  unsigned origBytesToIEND; //saved between WriteMolecule calls

  //Read and write number consisting of 4 bytes with most significant bytes first.
  //Should be independent of compiler and platform.
  unsigned long Read32(istream& ifs)
  {
    char ch;
    unsigned long val=0;
    for(int i=0; i<4; ++i)
    {
      if(!ifs.get(ch))
        return 0;
      val = val * 0x100 + (unsigned char)ch;
    }
    return val;
  }

  void Write32(unsigned long val, ostream& ofs)
  {
    char p[4];
    for(int i=0; i<4; ++i)
    {
      p[3-i] = (char)val % 0x100;
      val /= 0x100;
    }
    ofs.write(p, 4);
  }
};  
  ////////////////////////////////////////////////////

//Make an instance of the format class
PNGFormat thePNGFormat;

/////////////////////////////////////////////////////////////////

bool PNGFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  istream& ifs = *pConv->GetInStream();
  if(pConv->IsFirstInput())
    count=0;

  const char pngheader[] = {-119,80,78,71,13,10,26,10,0};
  char readbytes[9];
  ifs.read(readbytes, 8);
  
  if(!equal(pngheader, pngheader+8, readbytes))
  {
    obErrorLog.ThrowError("PNG Format","Not a PNG file", obError);
     return false;
  }

  //Loop through all the chunks
  while(ifs)
  {
    unsigned int len = Read32(ifs);
    ifs.read(readbytes,4);
    string chunkid(readbytes, readbytes+4);
    if(chunkid=="IEND")
    {
      bytesToIEND = ifs.tellg();
      bytesToIEND -= 8;
      break;
    }
    streampos pos = ifs.tellg();

    const char* altid = pConv->IsOption("y",OBConversion::INOPTIONS);
    if(chunkid=="tEXt" || chunkid=="zTXt" || (altid && chunkid==altid))
    {
      string keyword;
      getline(ifs, keyword, '\0');
      unsigned int datalength = len - keyword.size()-1;

      //remove "file" from end of keyword
      transform(keyword.begin(),keyword.end(),keyword.begin(),::tolower);
      string::size_type pos = keyword.find("file");
      if(pos!=string::npos)
        keyword.erase(pos);

      OBFormat* pFormat = OBConversion::FindFormat(keyword.c_str());
      if(pFormat)
      {
        //We have found embedded text that we need to extract
        stringstream ss;
        if(chunkid[0]!='z')
        {
          //Copy it to a stringstream
          istreambuf_iterator<char> initer(ifs);
          ostreambuf_iterator<char> outiter(ss);
          for(int i=0; i<datalength; ++i)
            *outiter++ = *initer++;
        }

        else
        {
          //Needs to be uncompressed first
          Bytef* pCompTxt = new Bytef[datalength];
          ifs.read((char*)pCompTxt, datalength);
          --datalength; //for compression method byte
          uLongf uncompLen;
          Bytef* pUncTxt = new Bytef[datalength*6];//guess uncompressed length. NASTY!
          if(*pCompTxt!=0 /*compression method*/
            || uncompress(pUncTxt, &uncompLen, pCompTxt+1, datalength)!=Z_OK)
          {
            obErrorLog.ThrowError("PNG Format","Errors in decompression", obError);
            return false;
          }
          pUncTxt[uncompLen] = '\0';
          ss.str((char*)pUncTxt);
          delete[] pUncTxt;
          delete[] pCompTxt;
        }

        //Use a new OBConversion object to convert embedded text
        OBConversion conv2(&ss, pConv->GetOutStream());
        conv2.CopyOptions(pConv);
        conv2.SetInAndOutFormats(pFormat, pConv->GetOutFormat());
        count += conv2.Convert();

        ifs.ignore(4);//CRC
        continue; //already at the end of the chunk
      }
    }
    //Move to end of chunk
    ifs.seekg(pos);
    ifs.ignore(len+4); //data + CRC
  }


  //if we will be writing a png file, read and save the whole input file.
  CopyOfInput.clear();
  if(pConv->GetOutFormat()==this)
  {
    ifs.seekg(0);
    copy(istreambuf_iterator<char>(ifs), istreambuf_iterator<char>(),back_inserter(CopyOfInput));
  }

  if(pConv->IsLastFile())
  {
    pConv->ReportNumberConverted(count); //report the number of chemical objects
    pConv->SetOutFormat(this); //so that number of files is reported as "PNG_files"
  }

  return true;
}

/////////////////////////////////////////////////////////////////
bool PNGFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  ostream& ofs = *pConv->GetOutStream();

  //Copy the saved input, except the IEND chunk, to the output
  if(!CopyOfInput.empty() && bytesToIEND>0)//check for input PNG not yet written
  {
    ostreambuf_iterator<char> outiter(pConv->GetOutStream()->rdbuf());
    //In Windows the output stream needs to be in binary mode to avoid extra CRs here
    copy(CopyOfInput.begin(), CopyOfInput.begin()+bytesToIEND, outiter);
    origBytesToIEND = bytesToIEND;
    bytesToIEND=0;//to ensure not copied again
  }

  //Convert pOb and write it to a tEXt chunk
  const char* formatID = pConv->IsOption("O", OBConversion::OUTOPTIONS);
  OBConversion conv2;
  conv2.CopyOptions(pConv); //So that can use commandline options in this conversion
  if(!formatID)
  {
    formatID="inchi";
    obErrorLog.ThrowError("PNG Format","Embedding in InChI format.\n"
      "Use the -xO (uppercase O) option for a different format", obWarning);
  }
  if(!conv2.SetOutFormat(OBConversion::FindFormat(formatID)))
  {
    obErrorLog.ThrowError("PNG Format","Format not found", obError);
    return false;
  
  }
  //Write new chunk
  stringstream ss;
  ss.str("");
  const char* pid = pConv->IsOption("y");
  if(pid && strlen(pid)==4)
    ss << pid;
  else
    ss  << "tEXt";
  ss  << formatID  << '\0'; 
  bool ret = conv2.Write(pOb, &ss);
  if(ret)
  {
    unsigned long len = ss.str().size() - 4; //don't count length of tEXt
    Write32(len, ofs);
    ofs << ss.str();

    //ss has type, keyword and data 
    uLong crc = crc32(0L, Z_NULL, 0);
    crc       = crc32(crc, (unsigned char*)ss.str().c_str(), ss.str().size());
    Write32(crc, ofs);

  }
  else
    obErrorLog.ThrowError("PNG Format","Failed when converting the molecule", obError);
    
  if(pConv->IsLast())
  {
    //Write the IEND chunk
    ostreambuf_iterator<char> outiter(pConv->GetOutStream()->rdbuf());
    copy(CopyOfInput.begin()+origBytesToIEND, CopyOfInput.end(), outiter);
    CopyOfInput.clear();

    // Decrement output index to not count the PNG file
    pConv->SetOutputIndex(pConv->GetOutputIndex()-1);
    // and report as output objects, not "PNG_files"
    pConv->SetOutFormat(formatID);
  }

  return ret;
}

/*
Reading
PNGFormat extracts chemical information that is embedded in PNG files.
The data can be in chunks of types tEXt, zTXt or, if in any type specified
with the -aa option. If the first letter of the type is 'z' the data is
decompressed.
The keyword in the chunk should be an OpenBabel Format ID, optionally with file added,
e.g. cml, InChI, molfile.
There can be multiple molecules in each chunk, multiple chunks with
chemical info and multiple png files can be read together.

Writing
This embeds chemical information into an existing PNG file.
A PNG file should be the first input file, followed by one or more chemical
files in any format. Each can contain multiple molecules. Each molecule is output
in a separate chunk in a format specified by the -xO option. the default with no
option is InChI. The chunk ID is normally tEXt but can be specified in the -xa option.
For example
  babel OrigImg.png Firstmol.smi Secondmol.mol2 OutImg.png -xO "cml" -xa "chEm"

It should be possible to embed into png filesusing the API.
The following is simplified and UNTESTED:

OBConversion conv;
conv.SetInAndOutFormats("png","png");
stringstream ss;
ifstream ifs("img.png");
ofstream ofs("img_with_chem.png");
OBMol mol;
conv.Read(&mol, &ifs); //Reads the embedded molecule
...manipulate mol
Note that the content of the PNG file is stored in PNGFormat, so do
not input from another PNG file until this one is written.

//Set the format of the embedded molecule on output
conv.AddOption("O",OBConversion::OUTOPTIONS,"smi"); 
conv.Write(&mol, ofs);

*/
} //namespace OpenBabel

