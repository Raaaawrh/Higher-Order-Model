#include "Mesh.h"

// ---------- // 1D // ----------

Mesh<1>::Mesh(double _minX, double _maxX, size_t _sizeX)
    : m_sizeX{_sizeX},
      m_minX{_minX}, m_maxX{_maxX},
      m_lengthX{m_maxX - m_minX},
      m_valuesX{vector<double>(m_sizeX, 0.0)}
{
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        m_valuesX[iX] = m_minX + iX * m_lengthX / (_sizeX - 1);
    }
}

// ---------- // 2D // ----------

Mesh<2>::Mesh(double _minX, double _minY,
              double _maxX, double _maxY,
              size_t _sizeX, size_t _sizeY)
    : m_sizeX{_sizeX}, m_sizeY{_sizeY},
      m_minX{_minX}, m_minY{_minY},
      m_maxX{_maxX}, m_maxY{_maxY},
      m_lengthX{m_maxX - m_minX}, m_lengthY{m_maxY - m_minY},
      m_valuesX{vector<vector<double>>(m_sizeX, vector<double>(m_sizeY, 0.0))},
      m_valuesY { vector<vector<double>>(m_sizeX, vector<double>(m_sizeY, 0.0)) }
{
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            m_valuesX[iX][iY] = m_minX + iX * m_lengthX / (m_sizeX - 1);
            m_valuesY[iX][iY] = m_minY + iY * m_lengthY / (m_sizeY - 1);
        }
    }
}

// ---------- // 3D // ----------

Mesh<3>::Mesh(double _minX, double _minY, double _minZ,
              double _maxX, double _maxY, double _maxZ,
              size_t _sizeX, size_t _sizeY, size_t _sizeZ)
    : m_sizeX{_sizeX}, m_sizeY{_sizeY}, m_sizeZ{_sizeZ},
      m_minX{_minX}, m_minY{_minY}, m_minZ{_minZ},
      m_maxX{_maxX}, m_maxY{_maxY}, m_maxZ{_maxZ},
      m_lengthX{m_maxX - m_minX}, m_lengthY{m_maxY - m_minY}, m_lengthZ{m_maxZ - m_minZ},
      m_valuesX{vector<vector<vector<double>>>(m_sizeX, vector<vector<double>>(m_sizeY, vector<double>(m_sizeZ, 0.0)))},
      m_valuesY{vector<vector<vector<double>>>(m_sizeX, vector<vector<double>>(m_sizeY, vector<double>(m_sizeZ, 0.0)))},
      m_valuesZ{vector<vector<vector<double>>>(m_sizeX, vector<vector<double>>(m_sizeY, vector<double>(m_sizeZ, 0.0)))}
{
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
            {
                m_valuesX[iX][iY][iZ] = m_minX + iX * m_lengthX / (m_sizeX - 1);
                m_valuesY[iX][iY][iZ] = m_minY + iY * m_lengthY / (m_sizeY - 1);
                m_valuesZ[iX][iY][iZ] = m_minZ + iZ * m_lengthZ / (m_sizeZ - 1);
            }
        }
    }
}