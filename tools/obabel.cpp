/**********************************************************************
main.cpp - Main conversion program, command-line handling.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004-2010 by Chris Morley

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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <map>
#if HAVE_CONIO_H
	#include <conio.h>
#endif
#include <cstdlib> // for exit() on Linux

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <openbabel/obconversion.h>
#include <openbabel/plugin.h>

using namespace std;
using namespace OpenBabel;

void DoOption(const char* p, OBConversion& Conv, OBConversion::Option_type typ,
	      int& arg, int argc, char *argv[]); 
void usage();
void help();

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{
  OBConversion Conv(&cin, &cout); //default input and output are console 

  OBFormat* pInFormat = NULL;
  OBFormat* pOutFormat = NULL;
  bool outGzip = false;
  vector<string> FileList, OutputFileList;
  string OutputFileName;

  // Parse commandline
  bool gotInType = false, gotOutType = false;
  bool SplitOrBatch=false;

  char *oext = NULL;
  char *iext = NULL;

  //Save name of program without its path (and .exe)
  string pn(argv[0]);
  string::size_type pos;
#ifdef _WIN32
  pos = pn.find(".exe");
  if(pos!=string::npos)
    argv[0][pos]='\0';
#endif
  pos = pn.find_last_of("/\\");
  if(pos==string::npos)
    program_name=argv[0];
  else
    program_name=argv[0]+pos+1;

  const char* p;
  int arg;
  for (arg = 1; arg < argc; ++arg)
    {
      if (argv[arg])
        {
          if (argv[arg][0] == '-')
            {
              char opchar[2]="?";
              opchar[0]=argv[arg][1];
              switch (opchar[0])
                {

                case 'V':
                  {
                    cout << "Open Babel " << BABEL_VERSION << " -- " 
                         << __DATE__ << " -- " << __TIME__ << endl;
                    exit(0);
                  }

                case 'i':
                  //Parameter is the input format which overrides any file extensions
                  gotInType = true;
                  iext = argv[arg] + 2;
                  if(!*iext)
                    iext = argv[++arg]; //space left after -i: use next argument

                  if (iext && strncasecmp(iext, "MIME", 4) == 0)
                    {
                      // get the MIME type from the next argument
                      iext = argv[++arg];
                      pInFormat = Conv.FormatFromMIME(iext);
                    }
                  else
                      pInFormat = Conv.FindFormat(iext);
                  if(pInFormat==NULL)
                    {
                      cerr << program_name << ": cannot read input format!" << endl;
                      usage();
                      exit(1);
                    }
                  break;

                case 'o':
                  //Parameter is the output format which overrides any file extension
                  gotOutType = true;
                  oext = argv[arg] + 2;
                  if(!*oext)
                    oext = argv[++arg]; //space left after -i: use next argument

                  if (oext && strncasecmp(oext, "MIME", 4) == 0)
                    {
                      // get the MIME type from the next argument
                      oext = argv[++arg];
                      pOutFormat = Conv.FormatFromMIME(oext);
                    }
                  else
                    pOutFormat = Conv.FindFormat(oext);

                  if(pOutFormat==NULL)
                    {
                      cerr << program_name << ": cannot write output format!" << endl;
                      usage();
                      exit(1);
                    }
                  break;

                case 'O':
                  OutputFileName = argv[arg] + 2;
                  if(OutputFileName.size()<3)
                    OutputFileName = argv[++arg]; //space left after -O: use next argument
                  break;

                case 'L': //display a list of plugin type or classes
                  {
                    const char* param=NULL;
                    if(argc>arg+1)
                      param = argv[arg+2];

                    // First assume first arg is a plugin type and
                    // param is a subtype, like babel -L ops gen3D
                    // or first arg is a plugin ID, like babel -L cml
                    OBPlugin* plugin;
                    if ((OBPlugin::GetPlugin("plugins", argv[arg+1]) &&
                         (plugin = OBPlugin::GetPlugin(argv[arg+1], param))) ||
                        (plugin = OBPlugin::GetPlugin(NULL, argv[arg+1])))
                    {
                      //Output details of subtype
                      string txt;
                      plugin->Display(txt, "verbose", argv[arg+1]);
                      cout << "One of the " << plugin->TypeID() << '\n' << txt << endl;
                      return 0;
                    }
                    //...otherwise assume it is a plugin type, like babel -L forcefields
                    //Output list of subtypes
                    OBPlugin::List(argv[arg+1], param);
                    return 0;
                  }

                case '?':
                case 'H':
                  if(isalnum(argv[arg][2]) || arg==argc-2)
                    {
                      if(strncasecmp(argv[arg]+2,"all",3))
                        {
                          OBFormat* pFormat
                            = (arg==argc-2) ? Conv.FindFormat(argv[arg+1]) : Conv.FindFormat(argv[arg]+2);
                          if(pFormat)
                            {
                              cout << argv[arg]+2 << "  " << pFormat->Description() << endl;
                              if(pFormat->Flags() & NOTWRITABLE)
                                cout << " This format is Read-only" << endl;
                              if(pFormat->Flags() & NOTREADABLE)
                                cout << " This format is Write-only" << endl;

                              if(strlen(pFormat->SpecificationURL()))
                                cout << "Specification at: " << pFormat->SpecificationURL() << endl;
                            }
                          else
                            cout << "Format type: " << argv[arg]+2 << " was not recognized" <<endl;
                        }
                      else
                        {
                          OBPlugin::List("formats","verbose");
                        }
                    }
                  else
                    help();
                  return 0;

                case '-': //long option --name text
                  {
                    //Option's text is in the next and subsequent args, until one starts with -
                    char* nam = argv[arg]+2;
                    if(!strcasecmp(nam, "help")) //special case handled here
                    {
                      help();
                      return 0;
                    }
                    if(*nam != '\0') //Do nothing if name is empty
                      {
                        string txt;
                        while(arg<argc-1 && *argv[arg+1]!='-')
                          {
                            //use text from subsequent args
                            if(!txt.empty())txt += ' '; //..space separated if more than one
                            txt += argv[++arg]; 
                          }

                        // If a API directive, e.g.---errorlevel
                        // send to the pseudoformat "obapi" (without any leading -)
                        if(*nam=='-')
                          {
                            OBConversion apiConv;
                            OBFormat* pAPI= OBConversion::FindFormat("obapi");
                            if(pAPI)
                              {
                                apiConv.SetOutFormat(pAPI);
                                apiConv.AddOption(nam+1, OBConversion::GENOPTIONS, txt.c_str());
                                apiConv.Write(NULL, &std::cout);
                              }
                          }
                        else
                          // Is a normal long option name, e.g --addtotitle
                          Conv.AddOption(nam,OBConversion::GENOPTIONS,txt.c_str());
                      }
                  }
                  break;
					
                case 'm': //multiple output files
                  SplitOrBatch=true;
                  break;
					
                case 'a': //single character input option
                  p = argv[arg]+2;
                  DoOption(p,Conv,OBConversion::INOPTIONS,arg,argc,argv);
                  break;

                case 'x': //single character output option
                  p = argv[arg]+2;
                  DoOption(p,Conv,OBConversion::OUTOPTIONS,arg,argc,argv);
                  break;

                //Not essential, but allows these options to be before input filenames
                //since we know they take one parameter, and are the most likely options to be misplaced
                case 'f':
                case 'l':
                  p = argv[arg] + 2;
                  if(!*p)
                    p = argv[++arg]; //space left after -f: use next argument
                  Conv.AddOption(opchar, OBConversion::GENOPTIONS, p);
                  break;
                
                case ':':
                  //e.g. -:c1ccccc1. SMILES passed as a file name and handled in OBConversion
                  FileList.push_back(argv[arg]);
                  break;

                default: //single character general option
                  p = argv[arg]+1;
                  DoOption(p,Conv,OBConversion::GENOPTIONS,arg,argc,argv);
                  break;
                }
            }
          else //filenames
              FileList.push_back(argv[arg]);
        }
    }

#if defined(_WIN32) && defined(USING_DYNAMIC_LIBS)
  //Expand wildcards in input filenames and add to FileList
  vector<string> tempFileList(FileList);
  FileList.clear();
  vector<string>::iterator itr;
  for(itr=tempFileList.begin();itr!=tempFileList.end();++itr)
  {
    if((*itr)[0]=='-')
      FileList.push_back(*itr);
    else
      DLHandler::findFiles (FileList, *itr);
  }
#endif
  
  if (!gotInType)
    {
      if(FileList.empty())
        {
          cerr << "No input file or format spec or possibly a misplaced option.\n"
            "Most options must come after the input files. (-i -o -O -m can be anywhwere.)\n" <<endl;
          usage();
          exit(1);
        }
    }

  if (!gotOutType)
    {
      //check there is a valid output format, but the extension will be re-interpreted in OBConversion
      pOutFormat = Conv.FormatFromExt(OutputFileName.c_str(), outGzip);
      if(OutputFileName.empty() || pOutFormat==NULL)
        {
          cerr << "Missing or unknown output file or format spec or possibly a misplaced option.\n"
            "Options, other than -i -o -O -m, must come after the input files.\n" <<endl;
          usage();
          exit(1);
        }
    }
  
    if(!Conv.SetInFormat(pInFormat)) //rely on autodetection for gzipped input
    {
      cerr << "Invalid input format" << endl;
      usage();
      exit(1);
    }
    if(!Conv.SetOutFormat(pOutFormat, outGzip))
    {
      cerr << "Invalid output format" << endl;
      usage();
      exit(1);
    }

  if(SplitOrBatch)
    {
      //Put * into output file name before extension (or ext.gz)
      if(OutputFileName.empty())
        {
          OutputFileName = "*.";
          OutputFileName += oext;
        }
      else
        {
          string::size_type pos = OutputFileName.rfind(".gz");
          if(pos==string::npos)
            pos = OutputFileName.rfind('.');
          else
            pos = OutputFileName.rfind('.',pos-1);
          if(pos==string::npos)
            OutputFileName += '*';
          else
            OutputFileName.insert(pos,"*");
        }
    }

  int count = Conv.FullConvert(FileList, OutputFileName, OutputFileList);
 
  Conv.ReportNumberConverted(count);

  if(OutputFileList.size()>1)
    {
      clog << OutputFileList.size() << " files output. The first is " << OutputFileList[0] <<endl;
    }

  //std::string messageSummary = obErrorLog.GetMessageSummary();
  //if (messageSummary.size())
  //  {
  //    clog << messageSummary << endl;
  //  }

#ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
#endif
  
  return 0;
}

void DoOption(const char* p, OBConversion& Conv,
	      OBConversion::Option_type typ, int& arg, int argc, char *argv[]) 
{
  //Unlike babel, cannot have multiple concatenated single char options
  //accepts: -sCCC -s CCC -s"CCC" -s CCC red -sCCC red
  char ch[2]="?";
  *ch = *p++;
  std::string txt;
  //Get the option text
  if(*p)
    txt = p; //use text immediately following the option letter, and keep looking

  while(arg<argc-1 && *argv[arg+1]!='-')
  {
    //use text from subsequent args
    if(!txt.empty())txt += ' '; //..space separated if more than one
    txt += argv[++arg]; 
  }
  Conv.AddOption(ch, typ, txt.c_str());
}

void usage()
{
  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage:\n" << program_name
       << " [-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]" << endl;
  cout << "Try  -H option for more information." << endl;

#ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
#endif
}

void help()
{
  cout << "Open Babel converts chemical structures from one file format to another"<< endl << endl;
  cout << "Usage: " << endl;
  cout << program_name << "[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]" << endl;
  cout << "The extension of a file decides the format, unless it is overridden" <<endl;
  cout << " by -i or -o options, e.g. -icml, or -o smi" << endl;
  cout << "See below for available format-types, which are the same as the " << endl;
  cout << "file extensions and are case independent." << endl; 
  cout << "If no input or output file is given stdin or stdout are used instead." << endl << endl; 
  cout << "More than one input file can be specified and their names can contain" <<endl;
  cout << "wildcard chars (* and ?). The format of each file can be different unless" <<endl;
  cout << "the -i option has been used, when they are all the same." <<endl;
  cout << "By default, the molecules are aggregated in the output file," << endl;
  cout << " but see -m option, Splitting, below.\n" << endl;

  cout << "Options, other than -i -o -O -m, must come after the input files.\n" <<endl;
  cout << OBConversion::Description(); // Conversion options
  cout << "-H Outputs this help text" << endl;
  cout << "-Hxxx (xxx is file format ID e.g. -Hcml) gives format info" <<endl; 
  cout << "-Hall Outputs details of all formats" <<endl; 
  cout << "-V Outputs version number" <<endl; 
  cout << "-L <category> Lists plugin classes of this category, e.g. <formats>" << endl;
  cout << "   Use just -L for a list of plugin categories." << endl; 
  cout << "   Use -L <ID> e.g. -L sdf for details of a format or other plugin." << endl; 
  cout << "-m Produces multiple output files, to allow:" <<endl;
  cout << "    Splitting: e.g.        " << program_name << " infile.mol -O new.smi -m" <<endl;
  cout << "      puts each molecule into new1.smi new2.smi etc" <<endl;
  cout << "    Batch conversion: e.g. " << program_name << " *.mol -osmi -m" <<endl;
  cout << "      converts each input file to a .smi file" << endl;
#ifdef _WIN32
  cout << "   In Windows these can also be done using the forms" <<endl;
  cout << "     " << program_name << " infile.mol -O new*.smi and " << program_name << " *.mol -O *.smi respectively.\n" <<endl;
#endif
  
  OBFormat* pDefault = OBConversion::GetDefaultFormat();
  if(pDefault)
    cout << pDefault->TargetClassDescription();// some more options probably for OBMol
  
  OBFormat* pAPI= OBConversion::FindFormat("obapi");
  if(pAPI)
    cout << pAPI->Description();
  
  cout << "To see a list of recognized file formats use\n  babel -L formats [read] [write]\n"
       << "To see details and specific options for a particular format, e.g CML, use\n  babel -L cml\n"
       << endl;
  //cout << "The following file formats are recognized:" << endl;
  //OBPlugin::List("formats");
  //cout << "\nSee further specific info and options using -H<format-type>, e.g. -Hcml" << endl;
}

/* OpenBabel man page*/
/** \page babel a converter for chemistry and molecular modeling data files
*
* \n
* \par SYNOPSIS
*
* \b babel [-H<help-options>] [-V] [-m] [-d] [-h] [-p] [-s<SMARTS-pattern>] [-v<SMARTS-pattern>] [-f<#> -l<#>] [-c] [-x<format-options>] [-i<input-type>] \<infile\> [-o<output-type>] -O\<outfile\>
*
* \par DESCRIPTION
*
* Open Babel is a program designed to interconvert a number of 
* file formats currently used in molecular modeling software. \n\n
*
* Note that Open Babel can also be used as a library to interconvert
* many file formats and to provide standard chemistry software routines.
* See the Open Babel web pages (http://openbabel.org) for more
* information.
*
* \par OPTIONS
*
* If only input and output files are given, Open Babel will guess 
* the file type from the filename extension. \n\n
*
* \b -V :
*     Output version number and exit \n\n
* \b -H :
*     Output usage information \n\n
* \b -H\<format-ID\> :
*     Output formatting information and options for format specified\n\n
* \b -Hall :
*     Output formatting information and options for all formats\n\n
* \b -i :
*     Specifies input format, see below for the available formats \n\n
* \b -o :
*     Specifies output format, see below for the available formats \n\n
* \b -m :
*     Produce multiple output files, to allow:\n
*      * Splitting one input file - put each molecule into consecutively numbered output files \n
*      * Batch conversion - convert each of multiple input files into a specified output format \n
*     See examples below \n\n
* \b -d : 
*     Delete Hydrogens \n\n
* \b -h : 
*     Add Hydrogens \n\n
* \b -p : 
*     Add Hydrogens appropriate for pH (use transforms in phmodel.txt)  \n\n
* \b -t :
*     All input files describe a single molecule \n\n
* \b -f\<#\> : 
*     For multiple entries input, start import at molecule # \n\n
* \b -l\<#\> : 
*     For multiple entries input, stop import at molecule # \n\n
* \b -c : 
*     Center atomic coordinates at (0,0,0) \n\n
* \b -s\<SMARTS\> :
*     Convert only molecules matching the SMARTS pattern specified \n\n
* \b -v\<SMARTS\> :
*     Convert only molecules \b NOT matching SMARTS pattern specified \n\n
*
* \par FILE FORMATS
*
* The following formats are currently supported by Open Babel:
*  \n    alc -- Alchemy format
*  \n    bgf -- BGF format     
*  \n    box -- Dock 3.5 Box format
*  \n    bs -- Ball and Stick format
*  \n    c3d1 -- Chem3D Cartesian 1 format
*  \n    c3d2 -- Chem3D Cartesian2 format
*  \n    caccrt -- Cacao format
*  \n    cache -- CAChe format [Writeonly]
*  \n    cacint -- CacaoInternal format [Writeonly]
*  \n    car -- MSI Biosym/Insight II CAR format [Readonly]
*  \n    ccc -- CCC format [Readonly]
*  \n    cht -- ChemTool format [Writeonly]
*  \n    cml -- Chemical Markup Language
*  \n    com -- Gaussian Cartesian Input [Writeonly]
*  \n    crk2d -- Chemical Resource Kit diagram format (2D)
*  \n    crk3d -- Chemical Resource Kit 3D format
*  \n    csr -- CSR format [Writeonly]
*  \n    cssr -- CSD CSSR format [Writeonly]
*  \n    ct -- ChemDraw Connection Table format 
*  \n    dmol -- DMol3 coordinates format
*  \n    ent -- Protein Data Bank format
*  \n    feat -- Feature format
*  \n    fh -- Fenske-Hall Z-Matrix format [Writeonly]
*  \n    fix -- FIX format [Writeonly]
*  \n    g03 -- Gaussian 98/03 Output [Readonly]
*  \n    g98 -- Gaussian 98/03 Output [Readonly]
*  \n    gam -- GAMESS Output
*  \n    gamout -- GAMESS Output
*  \n    gau -- Gaussian Cartesian Input [Writeonly]
*  \n    gjc -- Gaussian Cartesian Input [Writeonly]
*  \n    gjf -- Gaussian Cartesian Input [Writeonly]
*  \n    gpr -- Ghemical format
*  \n    gr96 -- GROMOS96 format [Writeonly]
*  \n    gzmat -- Gaussian Z-Matrix Input
*  \n    hin -- HyperChem Input format
*  \n    ins -- ShelX format [Readonly]
*  \n    jout -- Jaguar output format
*  \n    mdl -- MDL MOL format
*  \n    mm1gp -- Ghemical format
*  \n    mm3 -- MM3 format [Writeonly]
*  \n    mmd -- MacroMod format
*  \n    mmod -- MacroMod format
*  \n    mol -- MDL MOL format
*  \n    mol2 -- Sybyl Mol2 format
*  \n    mopcrt -- MOPAC Cartesian format
*  \n    mopout -- MOPAC Output format [Readonly]
*  \n    mpqc -- MPQC format [Readonly]
*  \n    nwo -- NWChem format
*  \n    pdb -- Protein Data Bank format
*  \n    pov -- POV-Ray input format [Writeonly]
*  \n    pqs -- Parallel Quantum Solutions format
*  \n    prep -- Amber Prep format [Readonly]
*  \n    qcout -- QChem output format
*  \n    qm1gp -- Ghemical format
*  \n    report -- Report format [Writeonly]
*  \n    res -- ShelX format [Readonly]
*  \n    rxn -- MDL RXN format
*  \n    sd -- MDL MOL format
*  \n    sdf -- MDL MOL format
*  \n    smi -- SMILES format
*  \n    tmol -- TurboMole Coordinate format
*  \n    txyz -- Tinker format [Writeonly]
*  \n    unixyz -- UniChem XYZ format
*  \n    vmol -- ViewMol format
*  \n    xed -- XED format [Writeonly]
*  \n    xyz -- XYZ cartesian coordinates format
*  \n    zin -- ZINDO input format [Writeonly]
*
* \par FORMAT OPTIONS
*  Individual file formats may have additional formatting options. \n
*  Input format options are preceded by 'a', e.g. -as \n
*  Output format options are preceded by 'x', e.g. -xn \n
*    For further specific information and options, use -H<format-type> \n
*    e.g., -Hcml
* 
* \par EXAMPLES
*  - Standard conversion \n
*     babel -ixyz ethanol.xyz -opdb ethanol.pdb \n
*  - Conversion from a SMI file in STDIN to a Mol2 file written to STDOUT \n
*     babel -ismi -omol2 \n
*  - Split a multi-molecule file into new1.smi, new2.smi, etc. \n
*     babel infile.mol new.smi -m \n
*
* \par AUTHORS
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at http://openbabel.org/
**/
