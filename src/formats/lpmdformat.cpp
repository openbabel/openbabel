/**********************************************************************
Copyright (C) 2012 GNM http://www.gnm.cl

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
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <cstdlib>

using namespace std;
namespace OpenBabel
{

class LpmdFormat : public OBMoleculeFormat
// Derive directly from OBFormat for objects which are not molecules.
{
 public:
  //Register this format type ID in the constructor
  LpmdFormat()
  {
   file_line=0;
   OBConversion::RegisterFormat("lpmd",this);
  }
  virtual const char* Description() //required
  {
   return
    "LPMD format\n"
    "Read and write LPMD's atomic configuration file\n\n"
    "Read Options e.g. -ab\n"
    "  s  Output single bonds only\n"
    "  b  Disable bonding entirely\n\n"
    "Write Options e.g. -xf 1\n"
    "  f# Indicate the level of the output file: 0 (default), 1 or 2.\n"
    "  m# Indicate the mode for level 2 output files\n"
    "        0 (default) is for accelerations and 1 for forces\n"
    "  c <vectorcells> Set the cell vectors if not present\n"
    "        Example: ``-xc 10.0,0,0,0.0,10.0,0.0,0.0,0.0,20.0``\n"
    "  e Add the charge to the output file\n\n"
    ;
  };

  //Optional URL where the file format is specified
  virtual const char* SpecificationURL()
  {return "http://www.lpmd.cl/index.php/documentation/the-lpmd-format";};

  //Optional
  virtual const char* GetMIMEType()  { return "chemical/lpmd"; };

  virtual int SkipObjects(int n, OBConversion* pConv)
  {
   return 0;
  };

  ////////////////////////////////////////////////////
  /// Declarations for the "API" interface functions. Definitions are below

  //Converts a string to a numerical type
  //This purloined from: http://www.codeguru.com/forum/showthread.php?t=231054
  template <class T>
  bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
  {
   std::istringstream iss(s);
   return !(iss >> f >> t).fail();
  }

  bool ReadHeader(std::istream &ifs, OBMol &mol);
  bool ReadAtom(std::istream &ifs, OBMol &mol);
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

 private:
  char buffer[BUFF_SIZE];
  std::vector<std::string> tokens;
  std::vector<std::string> headers; //Header information.
  int N; //the number of atoms
  int file_line;
  std::vector< vector3 > forces;
  std::vector< vector3 > accele;
  std::vector< vector3 > veloci;
};

LpmdFormat theLpmdFormat;

bool LpmdFormat::ReadHeader( std::istream &ifs, OBMol &mol )
{
 //Header Line
 if( ! ifs.getline(buffer,BUFF_SIZE) )
 {
  obErrorLog.ThrowError(__FUNCTION__,"Problem reading header line",obWarning);
  return false;
 }
 tokenize(tokens, buffer," ");
 if(tokens.size() == 0)
 {
  obErrorLog.ThrowError(__FUNCTION__,"The initial line it is empty!!! non LPMD format",obError);
  return false;
 }
 if(tokens.at(0).compare("LPMD") != 0 || tokens.at(1).compare("2.0") != 0)
 {
  obErrorLog.ThrowError(__FUNCTION__,"The start line, doesn't identify this file like a lpmd 2.0 file",obError);
  return false;
 }
 if(tokens.size()==3 && tokens.at(2).compare("Z")==0)
 {
  obErrorLog.ThrowError(__FUNCTION__,"There is not support for zipped files yet.",obError);
  return false;
 }
 if ( ! ifs.getline(buffer,BUFF_SIZE) )
 {
  obErrorLog.ThrowError(__FUNCTION__,"Problem reading header line",obError);
  return false;
 }
 tokenize(headers, buffer, " ");
 if(headers.size()<=1 || headers.at(0).compare("HDR")!=0)
 {
  obErrorLog.ThrowError(__FUNCTION__,"Problem reading header, check the HDR line",obError);
 }
 file_line = 2;
 return true;
} // End ReadHeader
/////////////////////////////////////////////////////////////////
bool LpmdFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
 OBMol* pmol = pOb->CastAndClear<OBMol>();
 if(pmol==NULL) return false;

 N=0;
 std::istream &ifs = *pConv->GetInStream();
 OBMol &mol = *pmol;

 //Read the header only once.
 if(file_line==0)
 {
  if ( ! ReadHeader( ifs, mol ) ) return false;
 }
 std::stringstream ErrMsg;
 //Read the cell vectors.
 double ax=0.0e0,ay=0.0e0,az=0.0e0;
 double bx=0.0e0,by=0.0e0,bz=0.0e0;
 double cx=0.0e0,cy=0.0e0,cz=0.0e0;
 //Reading the number of atoms.
 ifs.getline(buffer,BUFF_SIZE);
 file_line++;
 tokenize(tokens, buffer, " ");
 if(tokens.size()!=1)
 {
  obErrorLog.ThrowError(__FUNCTION__,"The number of atoms line was not correctly readed",obError);
 }
 from_string<int>(N, tokens.at(0), std::dec);

 tokens.clear();
 ifs.getline(buffer,BUFF_SIZE);
 file_line++;
 tokenize(tokens, buffer, " ");
 from_string<double>(ax, tokens.at(0), std::dec);
 from_string<double>(ay, tokens.at(1), std::dec);
 from_string<double>(az, tokens.at(2), std::dec);
 from_string<double>(bx, tokens.at(3), std::dec);
 from_string<double>(by, tokens.at(4), std::dec);
 from_string<double>(bz, tokens.at(5), std::dec);
 from_string<double>(cx, tokens.at(6), std::dec);
 from_string<double>(cy, tokens.at(7), std::dec);
 from_string<double>(cz, tokens.at(8), std::dec);

 vector3 vx = vector3(ax,ay,az);
 vector3 vy = vector3(bx,by,bz);
 vector3 vz = vector3(cx,cy,cz);
 OBUnitCell * unitcell = new OBUnitCell();
 unitcell->SetData( vx, vy, vz );

 //////////////////////////////////////////////////////
 //Start molecule modification.////////////////////////
 //////////////////////////////////////////////////////
 mol.BeginModify();
 mol.SetData( unitcell );

 for(int i=0 ; i<N ; ++i)
 {
  OBAtom *atom  = mol.NewAtom();
  ifs.getline(buffer,BUFF_SIZE);
  file_line++;
  tokenize(tokens, buffer, " ");
  if(int(headers.size()-1)!=tokens.size())
  {
   ErrMsg << "There was a problem reading an atomic configuration,  "
          << "the line # " << file_line << " doesn't have the number "
	  << "of columns indicated in the HDR (second line). ";
   obErrorLog.ThrowError(__FUNCTION__, ErrMsg.str(), obError);
  }
  //Based in the OBAtom class, only X Y Z VX VY VZ FX FY FZ
  //FCH PCH CHG are the first candidates to use from a LPMD file.
  double X=0.0,Y=0.0,Z=0.0, VX=0.0, VY=0.0, VZ=0.0;
  double AX=0.0, AY=0.0, AZ=0.0, FX=0.0, FY=0.0, FZ=0.0;
  double CHG=0.0,FCH=0.0;
  std::string symbol = "Xx";
  int ma=0,mf=0; //mode for acceleration or forces.
  for(int i=1 ; i < headers.size() ; ++i) //the first element is HDR
  {
   //Basic information position and atomic symbol.
   if(headers.at(i).compare("X")==0) from_string<double>(X, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("Y")==0) from_string<double>(Y, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("Z")==0) from_string<double>(Z, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("SYM")==0) symbol = tokens.at(i-1);
   //Velocities
   if(headers.at(i).compare("VX")==0) from_string<double>(VX, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("VY")==0) from_string<double>(VY, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("VZ")==0) from_string<double>(VZ, tokens.at(i-1), std::dec);
   //Accelerations
   if(headers.at(i).compare("AX")==0) {from_string<double>(AX, tokens.at(i-1), std::dec);ma++;}
   if(headers.at(i).compare("AY")==0) {from_string<double>(AY, tokens.at(i-1), std::dec);ma++;}
   if(headers.at(i).compare("AZ")==0) {from_string<double>(AZ, tokens.at(i-1), std::dec);ma++;}
   //Forces
   if(headers.at(i).compare("FX")==0) {from_string<double>(FX, tokens.at(i-1), std::dec);mf++;}
   if(headers.at(i).compare("FY")==0) {from_string<double>(FY, tokens.at(i-1), std::dec);mf++;}
   if(headers.at(i).compare("FZ")==0) {from_string<double>(FZ, tokens.at(i-1), std::dec);mf++;}
   //Charges
   if(headers.at(i).compare("CHG")==0) from_string<double>(CHG, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("PHG")==0) from_string<double>(CHG, tokens.at(i-1), std::dec);
   if(headers.at(i).compare("FHG")==0) from_string<double>(FCH, tokens.at(i-1), std::dec);
  }
  atom->SetVector(unitcell->FractionalToCartesian(vector3(X,Y,Z)));
  int atomicNum = OBElements::GetAtomicNum(symbol.c_str());
  atom->SetAtomicNum(atomicNum);
  //Conditional or zero??
  if( CHG!=0.0e0 ) atom->SetPartialCharge(CHG);
  //Always fill it?
  if(ma!=0 && mf ==0)
  {
   double mass = atom->GetAtomicMass();
   double FX = AX/mass; double FY = AY/mass; double FZ = AZ/mass;
  }
  forces.push_back(vector3(FX,FY,FZ));
  veloci.push_back(vector3(VX,VY,VZ));
 }

 OBConformerData * conformer = new OBConformerData();
 std::vector< std::vector< vector3 > > forceslist;
 std::vector< std::vector< vector3 > > velocilist;
 forceslist.push_back( forces );
 velocilist.push_back( veloci );
 conformer->SetForces( forceslist );
 conformer->SetVelocities( velocilist );
 mol.SetData(conformer);

 while(ifs.peek() != EOF && ifs.good() &&
       (ifs.peek() == '\n' || ifs.peek() == '\r'))
   ifs.getline(buffer,BUFF_SIZE);

 if (!pConv->IsOption("b",OBConversion::INOPTIONS))
    mol.ConnectTheDots();
 if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
    mol.PerceiveBondOrders();

 mol.EndModify();
 return true;
}

////////////////////////////////////////////////////////////////

bool LpmdFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
 OBMol* pmol = dynamic_cast<OBMol*>(pOb);
 if(pmol==NULL) return false;

 ostream& ofs = *pConv->GetOutStream();
 OBMol &mol = *pmol;
 OBUnitCell myUC;
 OBUnitCell *uc = NULL;

 std::vector< std::vector< vector3 > > forceslist;
 std::vector< std::vector< vector3 > > velocilist;

 std::stringstream ErrMsg;
 std::vector < vector3 > v;

 //Vectors read from command line.
 const char *c = pConv->IsOption("c");
 if(c)
 {
  double ax=0.0,ay=0.0,az=0.0;
  double bx=0.0,by=0.0,bz=0.0;
  double cx=0.0,cy=0.0,cz=0.0;
  tokenize(tokens, c, ",");
  from_string<double>(ax, tokens.at(0), std::dec);
  from_string<double>(ay, tokens.at(1), std::dec);
  from_string<double>(az, tokens.at(2), std::dec);
  from_string<double>(bx, tokens.at(3), std::dec);
  from_string<double>(by, tokens.at(4), std::dec);
  from_string<double>(bz, tokens.at(5), std::dec);
  from_string<double>(cx, tokens.at(6), std::dec);
  from_string<double>(cy, tokens.at(7), std::dec);
  from_string<double>(cz, tokens.at(8), std::dec);
  v.push_back(vector3(ax,ay,az));
  v.push_back(vector3(bx,by,bz));
  v.push_back(vector3(cx,cy,cz));
  myUC.SetData(v[0],v[1],v[2]);
  uc = &myUC;
 }

 if(mol.HasData(OBGenericDataType::UnitCell))
 {
  v.clear();
  uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  v = uc -> GetCellVectors();
 }
 else if(v.size()!=3)
 {
  ErrMsg << "\n*******************************************************************\n"
   << "Molecule -> " << pConv->GetOutputIndex()<< " : \n"
   << "The original file doesn't have the information about the unitcell\n"
   << "*******************************************************************\n";
  obErrorLog.ThrowError(__FUNCTION__, ErrMsg.str(), obError);
  return false;
 }

 //Define level and mode [(f)orce or (a)cceleration] and extra [(c)harge].
 int level=0, mode=0;
 bool charge = false;
 const char* p = pConv->IsOption("f");
 const char* q = pConv->IsOption("m");
 const char* r = pConv->IsOption("e");
 if(p) level = atoi(p);
 if(q) mode  = atoi(q);
 if(r) charge = true;

 if(pConv->GetOutputIndex()==1)
 {
  ofs << "LPMD 2.0 L\n" ;
  ofs << "HDR X Y Z ";
  if(level>=1) ofs << "VX VY VZ ";
  if(level>=2 && mode==0) ofs << "AX AY AZ ";
  if(level>=2 && mode==1) ofs << "FX FY FZ ";
  if(charge) ofs << "CHG ";
  ofs << '\n';
 }

 snprintf(buffer, BUFF_SIZE, "%d\n", mol.NumAtoms());
 ofs << buffer;

 //Fill velocities and forces.
 if(mol.HasData(OBGenericDataType::ConformerData))
 {
  OBConformerData * conformer = (OBConformerData*)mol.GetData(OBGenericDataType::ConformerData);
  forceslist = conformer->GetForces();
  velocilist = conformer->GetVelocities();
 }
 else
 {
  std::vector < vector3 > empty;
  empty.push_back(vector3(0,0,0));
  empty.push_back(vector3(0,0,0));
  empty.push_back(vector3(0,0,0));
  for(int i=0; i < mol.NumAtoms() ; ++i)
  {
   forceslist.push_back(empty);
   velocilist.push_back(empty);
  }
 }

 //Cell vectors
 snprintf(buffer, BUFF_SIZE, "%-10.3f%-10.3f%-10.3f", v[0].x(),v[0].y(),v[0].z());
 ofs << buffer;
 snprintf(buffer, BUFF_SIZE, "%-10.3f%-10.3f%-10.3f", v[1].x(),v[1].y(),v[1].z());
 ofs << buffer;
 snprintf(buffer, BUFF_SIZE, "%-10.3f%-10.3f%-10.3f\n", v[2].x(),v[2].y(),v[2].z());
 ofs << buffer;

 //Iteration over atoms
 for(int i=0;i<mol.NumAtoms();++i)
 {
  OBAtom *atom = mol.GetAtom(i + 1);
  vector3 tmp=uc->CartesianToFractional(vector3(atom->GetX(),atom->GetY(),atom->GetZ()));
  snprintf(buffer, BUFF_SIZE, "%-3s%15.5f%15.5f%15.5f",OBElements::GetSymbol(atom->GetAtomicNum()),
   tmp.GetX(),
   tmp.GetY(),
   tmp.GetZ());
  ofs << buffer;
  if(level>=1)
  {
   vector3 v = vector3(velocilist[0][i].x(),velocilist[0][i].y(),velocilist[0][i].z());
   snprintf(buffer, BUFF_SIZE, "%15.5f%15.5f%15.5f",v.x(),v.y(),v.z());
   ofs << buffer;
  }
  if(level>=2 && mode==0)
  {
   double mass = atom -> GetAtomicMass();
   vector3 a = vector3(forceslist[0][i].x()/mass,forceslist[0][i].y()/mass,forceslist[0][i].z()/mass);
   snprintf(buffer, BUFF_SIZE, "%15.5f%15.5f%15.5f",a.x(),a.y(),a.z());
   ofs << buffer;
  }
  if(level>=2 && mode==1)
  {
   vector3 force = vector3(forceslist[0][i].x(),forceslist[0][i].y(),forceslist[0][i].z());
   snprintf(buffer, BUFF_SIZE, "%15.5f%15.5f%15.5f",force.x(),force.y(),force.z());
   ofs << buffer;
  }
  if(charge)
  {
   snprintf(buffer, BUFF_SIZE, "%15.5f", atom -> GetPartialCharge());
   ofs << buffer;
  }
  ofs << '\n';
 }

 return(true);
}

} //namespace OpenBabel
