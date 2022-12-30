#include "Field.h"

//// 1D

template <typename T>
Field<T, 1>::Field()
    : Field<T, 1>(0){};

template <typename T>
Field<T, 1>::Field(size_t _sizeX)
    : m_Nx{_sizeX},
      m_vals{vector<T>(m_Nx, static_cast<T>(0))} {};

template <typename T>
size_t Field<T, 1>::sizeX() const { return m_Nx; }

template <typename T>
T &Field<T, 1>::operator()(size_t _indexX) { return m_vals[_indexX]; }

//// 2D

template <typename T>
Field<T, 2>::Field()
    : Field<T, 2>(0, 0){};

template <typename T>
Field<T, 2>::Field(size_t _sizeX, size_t _sizeY)
    : m_Nx{_sizeX}, m_Ny{_sizeY},
      m_vals{vector<vector<T>>(m_Nx, vector<T>(m_Ny, static_cast<T>(0)))} {};

template <typename T>
size_t Field<T, 2>::sizeX() const { return m_Nx; }

template <typename T>
size_t Field<T, 2>::sizeY() const { return m_Ny; }

template <typename T>
T &Field<T, 2>::operator()(size_t _indexX, size_t _indexY) { return m_vals[_indexX][_indexY]; }

//// 3D

template <typename T>
Field<T, 3>::Field()
    : Field<T, 3>(0, 0, 0){};

template <typename T>
Field<T, 3>::Field(size_t _sizeX, size_t _sizeY, size_t _sizeZ)
    : m_Nx{_sizeX}, m_Ny{_sizeY}, m_Nz{_sizeZ},
      m_vals{vector<vector<vector<T>>>(m_Nx, vector<vector<T>>(m_Ny, vector<T>(m_Nz, static_cast<T>(0))))} {};

template <typename T>
size_t Field<T, 3>::sizeX() const { return m_Nx; }

template <typename T>
size_t Field<T, 3>::sizeY() const { return m_Ny; }

template <typename T>
size_t Field<T, 3>::sizeZ() const { return m_Nz; }

template <typename T>
T &Field<T, 3>::operator()(size_t _indexX, size_t _indexY, size_t _indexZ) { return m_vals[_indexX][_indexY][_indexZ]; }
