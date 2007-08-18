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
    OBConversion::RegisterFormat("PNG",this);
    OBConversion::RegisterOptionParam("a", this, 0, OBConversion::INOPTIONS);
  }

  virtual const char* Description()
  {
    return
    "PNG format\n"
    "Extract chemical structure data embedded in PNG image files\n"
    "Read Options e.g. -aa\n"
    " a <additional chunk ID>Look also in chunks with specified ID\n\n";
  };

  virtual const char* TargetClassDescription(){return " PNG_files";}

  virtual unsigned int Flags()
  {
      return NOTWRITABLE | READONEONLY | READBINARY;
  };

  virtual bool ReadChemObject(OBConversion* pConv)
  {
    bool ret = ReadMolecule(NULL, pConv);
    pConv->GetChemObject(); //increments output index
    return ret;
  };

virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

private:
  int count; //number of chemical objects converted
  struct chunkheader
  {
    long length;
    char type[4];
  };
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

  char pngheader[] = {-119,80,78,71,13,10,26,10};
  char readbytes[9];
  ifs.read(readbytes, 8);
  
  if(!equal(pngheader, pngheader+8, readbytes)) //This gives a warning about a "MS deprecated" function, in spite of the warning being supressed.
  {
    obErrorLog.ThrowError("PNG Format","Not a PNG file", obError);
     return false;
  }

  //Loop through all the chunks
  while(ifs)
  {
    chunkheader chh;
    ifs.read((char*)&chh, 8);
    char* p = (char*)&chh.length;
    #ifndef WORDS_BIGENDIAN
      swap(p[0], p[3]);
      swap(p[1], p[2]);
    #endif
    string chunkid(p+4, p+8);
    if(chunkid=="IEND")
      break;
    streampos pos = ifs.tellg();

    const char* altid = pConv->IsOption("a",OBConversion::INOPTIONS);
    if(chunkid=="tEXt" || chunkid=="zTXt" || (altid && chunkid==altid))
    {
      string keyword;
      getline(ifs, keyword, '\0');
      int datalength = chh.length - keyword.size()-1;

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
    ifs.ignore(chh.length+4); //data + CRC
  }

  if(pConv->IsLastFile())
  {
    pConv->ReportNumberConverted(count); //report the number of chemical objects
    pConv->SetOutFormat(this); //so that number of files is reported as "PNG_files"
  }
  return true;
}

////////////////////////////////////////////////////////////////


} //namespace OpenBabel

