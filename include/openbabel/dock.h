/**********************************************************************
dock.h - OBDock class.
 
Copyright (C) 2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#ifndef OB_DOCK_H
#define OB_DOCK_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  // score functions
#define OBSCORE_GHEMICAL	(1 << 0)   //!< score the poses using the Ghemical force field (poor results)
#define OBSCORE_MMFF94		(2 << 0)   //!< score the poses using the MMFF94 force field (this is probably what you want)
#define OBSCORE_UFF		(3 << 0)   //!< score the poses using the UFF force field (don't know about results)
#define OBSCORE_XSCORE 		(4 << 0)   //!< score the poses using the XScore empirical scoring function
  // conformer generation
#define OBDOCK_INITIAL		(1 << 0)   //!< use the initial,optimized (for optimize_steps) ligand coordinates, 
#define OBDOCK_GEN_SYST		(2 << 0)   //!< generate coordinates for the ligand using systematic rotor search
#define OBDOCK_GEN_RAND		(2 << 0)   //!< generate coordinates for the ligand using random rotor search
#define OBDOCK_GEN_WGHT		(2 << 0)   //!< generate coordinates for the ligand using weighted rotor search
  
  //! \class OBDock
  //! \brief Docking Class
  class OBDock {
    public:
      //! constructor
      OBDock();
      //! destructor
      ~OBDock();
      /*! Create a ligand from all atoms of residue named resname.
       *  \param mol OBMol object to use as input
       *  \param resname residue name of the ligand 
       *  \return vector with atom indices for the binding pocket atoms
       */
      static std::vector<int> CreateLigandFromResidue(OBMol &mol, std::string &resname); 
      /*! Create a binding pocket from all atoms within radius from atoms of the ligand.
       *  \param mol OBMol object to use as input
       *  \param resname residue name of the ligand 
       *  \param radius radius
       *  \return vector with atom indices for the binding pocket atoms
       */
      static std::vector<int> CreatePocketFromLigand(OBMol &mol, std::vector<int> &ligand , double radius);
      /*! Create a binding pocket from all atoms within radius from atoms of residue named resname.
       *  \param mol OBMol object to use as input
       *  \param resname residue name of the ligand 
       *  \param radius radius
       *  \return vector with atom indices for the binding pocket atoms
       */
      static std::vector<int> CreatePocketFromResidue(OBMol &mol, std::string &resname , double radius);
      
      //////////////////////////////////////////////////////////////////////////
      //
      //  Grid docking
      //
      //////////////////////////////////////////////////////////////////////////
      //! \name Grid Docking
      //@{
      /*! Initialize a grid docking for the ligand in the binding pocket. 
       *
       *  \param mol OBMol object that contains ligand and pocket (could also contain more ignored atoms)
       *  \param ligand indices for the ligand
       *  \param pocket indices for the binding pocket
       *  \param optimize_steps number of geometry optimization steps for each generated conformer or pose
       *  \param translate_step the translate step size
       *  \param rotate_step the rotate step size
       */
      void GridDockInitialize(OBMol &mol, std::vector<int> &ligand, std::vector<int> &pocket, 
          int optimize_steps, double translate_step, double rotate_step);
      /*! Generate the next pose, assign a score and store it. 
       *  \return true if there are more poses
       */
      bool GridDockNextPose();
      //@}
      //////////////////////////////////////////////////////////////////////////
      //
      //  Get informeration on poses
      //
      //////////////////////////////////////////////////////////////////////////
      //! \name Grid Docking
      //@{
      //! \return the number of generated poses.
      int NumPoses();
      //! \return the number of generated poses.
      //@}
 
    private:
      OBMol 			_mol; //!< OBMol object that holds the atoms for the ligand, binding pocket and possibly more ignored atoms
      std::vector<int>		_ligand; //!< indices for the ligand
      std::vector<int>		_pocket; //!< indices for the binding pocket
      vector3 			_pocket_dimentions; //!< pocket dimentions, constant
      vector3 			_ligand_dimentions; //!< ligand dimentions, conformer and rotation specific
      int			_optimize_steps; //!< pose & conformer optimization steps
      double			_translate_step; //!< translate step size
      double			_rotate_step; //!< rotate step size
      // iterate variables
      vector3			_translate; //!< hold current translation
      vector3			_rotate; //!< hold the current rotation
      int  			_conformer; //!< hold the current conformer
      //////////////////////////////////////////////////////////////////////////
      //
      //  Privare method for internal use
      //
      //////////////////////////////////////////////////////////////////////////
      //! \name Private methods for internal use
      //@{
      /*! Get the dimentions for the selected indices.
       *  \return the number of generated poses.
       */
      vector3 GetDimentions(std::vector<int> &selection);
      //@}
 
  }; // class OBDock

}// namespace OpenBabel

#endif   // OB_DOCK_H

//! \file dock.h
//! \brief Docking Class
