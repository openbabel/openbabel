/**********************************************************************
fileformat.cpp - Read and write file formats.

Copyright (C) 2000-2003 by Geoffrey Hutchison

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

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "fileformat.h"
#include "mol.h"

using namespace std;

namespace OpenBabel {

bool OBFileFormat::ReadMolecule(istream &ifs, OBMol &mol, const char *title)
{
  bool result;

  if (!ifs)
    return false;
  
  switch(mol.GetInputType())
    {
    case ALCHEMY:   result = ReadAlchemy(ifs,mol,title);	break;
    case BALLSTICK: result = ReadBallAndStick(ifs,mol,title);	break;
    case BGF:	    result = ReadBGF(ifs,mol,title);		break;
    case BIOSYM:    result = ReadBiosymCAR(ifs,mol,title);	break;
    case BOX:       result = ReadBox(ifs,mol,title);		break;
    case CACAO:	    result = ReadCaccrt(ifs,mol,title);	       	break;
    case CCC:       result = ReadCCC(ifs,mol,title);	       	break;
    case CHEM3D1:   result = ReadChem3d1(ifs,mol,title);       	break;
    case CHEM3D2:   result = ReadChem3d2(ifs,mol,title);       	break;
    case CML:       result = ReadCML(ifs,mol,title);		break;
    case CRK2D:     result = ReadCRK2D(ifs,mol,title);          break;
    case CRK3D:     result = ReadCRK3D(ifs,mol,title);          break;
    case DMOL:      result = ReadDMol(ifs,mol,title);		break;
    case FEATURE:   result = ReadFeat(ifs,mol,title);		break;
    case GAMESSOUT: result = ReadGAMESS(ifs,mol,title);		break;
    case GHEMICAL:  result = ReadGhemical(ifs,mol,title);	break; 
    case HIN:	    result = ReadHIN(ifs,mol,title);		break;
    case NWCHEMOUT: result = ReadNWChem(ifs, mol, title);	break;
    case MMD:       result = ReadMacroModel(ifs,mol,title);	break;
    case MMADS:     result = ReadMmads(ifs,mol,title);		break;
    case MOL2:      result = ReadMol2(ifs,mol,title);		break;
    case MOPACOUT:  result = ReadMOPAC(ifs,mol,title);		break;
    case MOPACCART: result = ReadMOPACCartesian(ifs,mol,title);	break;
    case MPQC:      result = ReadMPQC(ifs,mol,title);		break;
    case OEBINARY:  result = ReadBinary(ifs,mol);		break;
    case PDB:       result = ReadPDB(ifs,mol,title);		break;
    case PREP:	    result = ReadAmberPrep(ifs,mol,title);	break;
    case PQS:       result = ReadPQS(ifs,mol,title);            break;
    case JAGUAROUT: result = ReadJaguar(ifs,mol,title);		break;
    case QCHEMOUT:  result = ReadQChem(ifs,mol,title);		break;
    case SDF:       result = ReadSDFile(ifs,mol,title);		break;
    case SMI:       result = ReadSmiles(ifs,mol,title);		break;
    case SHELX:	    result = ReadShelX(ifs,mol,title);		break;
    case UNICHEM:   result = ReadUnichem(ifs,mol,title);	break;
    case VIEWMOL:   result = ReadViewMol(ifs,mol,title);	break;
    case XYZ:	    result = ReadXYZ(ifs,mol,title);		break;
    default:
      ThrowError("Input type not defined");
      return false;
    }
  
  return result;
}

bool OBFileFormat::WriteMolecule(ostream &ofs,OBMol &mol, 
				 const char *dimension, const char *options)
{
  switch(mol.GetOutputType())
    {
    case ALCHEMY:   WriteAlchemy(ofs,mol);		break;
    case BALLSTICK: WriteBallAndStick(ofs,mol);		break;
    case BGF:	    WriteBGF(ofs,mol);			break;
    case BOX:	    WriteBox(ofs,mol,1.0);		break;
    case CACAO:     WriteCaccrt(ofs,mol);		break;
    case CACAOINT:  WriteCacaoInternal(ofs,mol);	break;
    case CACHE:     WriteCache(ofs,mol);		break;
    case CHEMDRAW:  WriteChemDraw(ofs,mol);		break;
    case CHEMTOOL:  WriteCHT(ofs,mol);			break;
    case CHEM3D1:   WriteChem3d1(ofs,mol);		break;
    case CHEM3D2:   WriteChem3d2(ofs,mol);		break;
    case CML:       WriteCML(ofs,mol,dimension, options);break;
    case CRK2D:     WriteCRK2D(ofs,mol);                break;
    case CRK3D:     WriteCRK3D(ofs,mol);                break;
    case CSR:       WriteCSR(ofs,mol);			break;
    case CSSR:      WriteCSSR(ofs,mol);			break;
    case DMOL:      WriteDMol(ofs,mol);			break;
    case DELPDB:    WriteDelphiPDB(ofs,mol);  		break;
    case FEATURE:   WriteFeat(ofs,mol);			break;
    case FH:	    WriteFenskeZmat(ofs,mol);		break;
    case FIX:       WriteFixFile(ofs,mol);    		break;
    case GAMESSIN:  WriteGAMESS(ofs,mol);		break;
    case GHEMICAL:  WriteGhemical(ofs,mol);   		break;
    case GROMOS96A: WriteGromos96A(ofs,mol);		break;
    case GROMOS96N: WriteGromos96N(ofs,mol);		break;
    case GAUSSIANCART:WriteGaussianCart(ofs,mol);	break;
    case HIN:	    WriteHIN(ofs,mol);			break;
    case JAGUARIN:  WriteJaguar(ofs,mol);		break;
    case OEBINARY:  WriteBinary(ofs,mol);     		break;
    case NWCHEMIN:  WriteNWChem(ofs, mol);		break;
    case MMD:       WriteMacroModel(ofs,mol); 		break;
    case MMADS:     WriteMmads(ofs,mol); 		break;
    case MOL2:      WriteMol2(ofs,mol,dimension);  	break;
    case MOPACCART: WriteMOPACCartesian(ofs,mol);	break;
    case PDB:       WritePDB(ofs,mol);			break;
#ifndef WIN32
    case POV:	    WritePovray(ofs,mol,options);       break;
#endif
    case PQS:       WritePQS(ofs,mol);                  break;
    case QCHEMIN:   WriteQChem(ofs,mol);		break;
    case REPORT:    WriteReport(ofs,mol);		break;
    case SDF:       WriteSDFile(ofs,mol,dimension);     break;
    case SMI:       WriteSmiles(ofs,mol);     		break;
    case TINKER:    WriteTinker(ofs,mol);		break;
    case TITLE:	    WriteTitles(ofs,mol); 	        break;
    case UNICHEM:   WriteUnichem(ofs,mol);		break;
    case VIEWMOL:   WriteViewMol(ofs,mol);		break;
    case XED:	    WriteXED(ofs,mol);			break;
    case XYZ:	    WriteXYZ(ofs,mol);			break;
    case ZINDO:	    WriteZINDO(ofs,mol);		break;

    default:
      ThrowError("Output type not defined");
    }

  return((ofs) ? true : false);
}

}
