#include "Field.h"

// ---------- // 1D // ----------

// Constructors and destructors
Field<1>::Field(size_t _sizeX)
    : m_sizeX{_sizeX},
      m_values{vector<double>(m_sizeX, 0.0)} {}

// Getters
size_t Field<1>::getSizeX() const { return m_sizeX; }

// Operators
double &
Field<1>::operator()(size_t _indexX) { return m_values[_indexX]; }

// ---------- // 2D // ----------

// Constructors and destructors
Field<2>::Field(size_t _sizeX, size_t _sizeY)
    : m_sizeX{_sizeX}, m_sizeY{_sizeY},
      m_values{vector<vector<double>>(m_sizeX, vector<double>(m_sizeY, 0.0))} {}

// Getters
size_t Field<2>::getSizeX() const { return m_sizeX; }
size_t Field<2>::getSizeY() const { return m_sizeY; }

// Operators
double &Field<2>::operator()(size_t _indexX, size_t _indexY) { return m_values[_indexX][_indexY]; }
double const Field<2>::operator()(size_t _indexX, size_t _indexY) const { return m_values[_indexX][_indexY]; }

// ---------- // 3D // ----------

// Constructors and destructors
Field<3>::Field(size_t _sizeX, size_t _sizeY, size_t _sizeZ)
    : m_sizeX{_sizeX}, m_sizeY{_sizeY}, m_sizeZ{_sizeZ},
      m_values{vector<vector<vector<double>>>(m_sizeX, vector<vector<double>>(m_sizeY, vector<double>(m_sizeZ, 0.0)))}
{
}

// Getters
size_t Field<3>::getSizeX() const { return m_sizeX; }
size_t Field<3>::getSizeY() const { return m_sizeY; }
size_t Field<3>::getSizeZ() const { return m_sizeZ; }

// Operators
double &Field<3>::operator()(size_t _indexX, size_t _indexY, size_t _indexZ) { return m_values[_indexX][_indexY][_indexZ]; }
double Field<3>::operator()(size_t _indexX, size_t _indexY, size_t _indexZ) const { return m_values[_indexX][_indexY][_indexZ]; }
