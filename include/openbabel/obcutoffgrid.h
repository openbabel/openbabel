/*   Zodiac - molecular modelling environment. www.zeden.org
    Copyright (C) 2008  Nicola Zonta (nicola.zonta(at)zeden.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef OB_CUTOFFGRID_H
#define OB_CUTOFFGRID_H

#include <vector>
#include <Eigen/Core>

namespace OpenBabel
{

  template <class T>
  struct staticCutoffNode {
    std::vector<T> cellObjects;
    std::vector<T>* objects;
  };

  template <class T>
  class OBCutoffGrid
  {
    public:
      OBCutoffGrid(std::vector<T>& objects, float cutRadius);
      ~OBCutoffGrid();
      std::vector<T>* GetNbrs(const Eigen::Vector3d &pos);
      std::vector<staticCutoffNode<T>*> grid;
      
      float gridResolution;
	
    private:
      int m_xDim, m_yDim, m_zDim;
      int m_xyDim;
      float cutoffRadius;
      Eigen::Vector3d m_gridMin;
      Eigen::Vector3d m_gridMax;
      unsigned int maxSID;
  };
 
};

#endif

