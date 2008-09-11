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

#include <openbabel/obcutoffgrid.h>

namespace OpenBabel
{
  template <class T>
  OBCutoffGrid<T>::OBCutoffGrid(std::vector<T>& objects, float cutRadius)
  {
    gridResolution = cutRadius;	
    cutoffRadius = cutRadius;
	
    if (!objects.size()) {
      m_gridMin = m_gridMax = Eigen::Vector3d::Zero();
      m_xDim = m_yDim = m_zDim = 0;
    } else {
      for (unsigned int i = 0; i < objects.size(); i++) {
        if (!i) {
          m_gridMin = m_gridMax = objects[i]->GetVector();
	} else {
          if (objects[i]->GetVector().x() < m_gridMin.x())
            m_gridMin.x() = objects[i]->GetVector().x();
          if (objects[i]->GetVector().y() < m_gridMin.y()) 
            m_gridMin.y() = objects[i]->GetVector().y();
          if (objects[i]->GetVector().z() < m_gridMin.z())
            m_gridMin.z() = objects[i]->GetVector().z();

          if (objects[i]->GetVector().x() > m_gridMax.x())
            m_gridMax.x() = objects[i]->GetVector().x();
          if (objects[i]->GetVector().y() > m_gridMax.y())
            m_gridMax.y() = objects[i]->GetVector().y();
          if (objects[i]->GetVector().z() > m_gridMax.z())
            m_gridMax.z() = objects[i]->GetVector().z();
	}
      }
    }

    m_gridMin -= 2.0 * cutRadius * Eigen::Vector3d::Ones();
    m_gridMax += 2.0 * cutRadius * Eigen::Vector3d::Ones();

    m_xDim = static_cast<int>((m_gridMax.x() - m_gridMin.x()) / gridResolution + 0.5);
    m_yDim = static_cast<int>((m_gridMax.y() - m_gridMin.y()) / gridResolution + 0.5);
    m_zDim = static_cast<int>((m_gridMax.z() - m_gridMin.z()) / gridResolution + 0.5);
		
    m_xyDim = m_xDim * m_yDim;

    grid.resize(m_xyDim * m_zDim, 0);

    std::vector<int*> nodes;

    for (unsigned int i = 0; i < objects.size(); i++) {
      Eigen::Vector3i gridPos;

      gridPos.x() = static_cast<int>((objects[i]->GetVector().x() - m_gridMin.x()) / gridResolution + 0.5);
      gridPos.y() = static_cast<int>((objects[i]->GetVector().x() - m_gridMin.x()) / gridResolution + 0.5);
      gridPos.z() = static_cast<int>((objects[i]->GetVector().x() - m_gridMin.x()) / gridResolution + 0.5);

      if (gridPos.x() >= 0 && gridPos.x() < m_xDim &&
          gridPos.y() >= 0 && gridPos.y() < m_yDim &&
          gridPos.z() >= 0 && gridPos.z() < m_zDim) 
      {
        int index = gridPos.x() + gridPos.y() * m_xDim + gridPos.z() * m_xyDim;
        if (!grid[index]) {
          grid[index] = new staticCutoffNode<T>;
          grid[index]->objects = new std::vector<T>;
          grid[index]->cellObjects.push_back(objects[i]);
          int *nodeIndex = new int[3];
          nodeIndex[0] = gridPos[0];
          nodeIndex[1] = gridPos[1];
          nodeIndex[2] = gridPos[2];

          nodes.push_back(nodeIndex);
        } else {	
          grid[index]->cellObjects.push_back(objects[i]);
        }
      }
    }

    int radius = 1;
    int x, y, z;

    std::vector<int*> sphere;

    for (x = -radius; x <= radius; x++)
      for (y = -radius; y <= radius; y++)
	for (z = -radius; z <= radius; z++) {
          int *p = new int[3];			
          p[0] = x;
          p[1] = y;
          p[2] = z;
          sphere.push_back(p);
        }
    for (unsigned int i = 0; i < nodes.size(); i++) {
      int *ni = nodes[i];
      int index = ni[0] + ni[1] * m_xDim + ni[2] * m_xyDim;
      for (unsigned int j = 0; j < sphere.size(); j++) {
        int *p = sphere[j];
        int newPoint[3];
        newPoint[0] = ni[0] + p[0];
        newPoint[1] = ni[1] + p[1];
        newPoint[2] = ni[2] + p[2];

        if (newPoint[0] >= 0 && newPoint[0] < m_xDim &&
            newPoint[1] >= 0 && newPoint[1] < m_yDim &&
            newPoint[2] >= 0 && newPoint[2] < m_zDim) 
        {
          int newIndex = newPoint[0] + newPoint[1] * m_xDim + newPoint[2] * m_xyDim;
          if (!grid[newIndex]) {
            grid[newIndex] = new staticCutoffNode<T>;
            grid[newIndex]->objects = new std::vector<T>;
          }
          for (unsigned int k = 0; k < grid[index]->cellObjects.size(); k++) {
            grid[newIndex]->objects->push_back(grid[index]->cellObjects[k]);
          }
        }
      }  
      
      delete[] nodes[i];
    }
    
    for (unsigned int i = 0; i < sphere.size(); i++) {
      if (sphere[i]) {
        delete[] sphere[i];
        sphere[i] = 0;
      }
    }
    sphere.clear();
  }


  template <class T>
  OBCutoffGrid<T>::~OBCutoffGrid()
  {
    for (unsigned int i = 0; i < grid.size(); i++) {
      if (grid[i]) {
        if (grid[i]->objects) {
          grid[i]->objects->clear();
          delete grid[i]->objects;
        }
       
        grid[i]->objects = 0;
        delete grid[i];
        grid[i] = 0;
      }
    }
    
    grid.clear();
  }

template <class T>
std::vector<T>* OBCutoffGrid<T>::GetNbrs(const Eigen::Vector3d &pos)
{
  int gridPos[3];
  
  gridPos[0] = (int)((pos.x() - m_gridMin[0]) / gridResolution + 0.5);
  gridPos[1] = (int)((pos.y() - m_gridMin[1]) / gridResolution + 0.5);
  gridPos[2] = (int)((pos.z() - m_gridMin[2]) / gridResolution + 0.5);
  
  if (gridPos[0] >= 0 && gridPos[0] < m_xDim &&
    gridPos[1] >= 0 && gridPos[1] < m_yDim &&
    gridPos[2] >= 0 && gridPos[2] < m_zDim) {
    
    int index = gridPos[0] + gridPos[1] * m_xDim + gridPos[2] * m_xyDim;
    if (grid[index]) {
      return grid[index]->objects;
    }
    else {
      return 0;
    }
  }
  return 0;
}

};

