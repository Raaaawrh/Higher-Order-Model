#if !defined(MESH_H)
#define MESH_H

#include <vector>

using std::size_t, std::vector;

// ---------- // Mesh class // ----------

template <size_t dimensions>
class Mesh
{
public:
    // Constructors and destructors
    Mesh(){};
    ~Mesh(){};
};

// ---------- // 1D // ----------

template <>
class Mesh<1>
{
public:
    // Constructors and destructors
    Mesh(){};
    Mesh(double _minX, double _maxX, size_t _sizeX /* reference to transform */);
    ~Mesh() {}

    // Getters
    std::size_t getSizeX() const { return m_sizeX; }
    double getMinX() const { return m_minX; }
    double getMaxX() const { return m_maxX; }
    double getLenghtX() const { return m_lengthX; }
    vector<double> getValuesX() const { return m_valuesX; }

    double getX(size_t _indexX) const { return m_valuesX[_indexX]; }

    double getDX(size_t _indexX) const { return m_valuesX[_indexX + 1] - m_valuesX[_indexX]; }

private:
    // Sizes: X
    size_t m_sizeX;

    // Limits: X
    double m_minX, m_maxX;

    // Length: X
    double m_lengthX;

    // Knots: X
    vector<double> m_valuesX;
};

// ---------- // 2D // ----------

template <>
class Mesh<2>
{
public:
    // Constructors and destructors
    Mesh(){};
    Mesh(double _minX, double _minY, double _maxX, double _maxY, size_t _sizeX, size_t _sizeY /* reference to transform */);
    ~Mesh(){};

    // Getters
    size_t getSizeX() const { return m_sizeX; }
    size_t getSizeY() const { return m_sizeY; }

    double getMinX() const { return m_minX; }
    double getMaxX() const { return m_maxX; }

    double getMinY() const { return m_minY; }
    double getMaxY() const { return m_maxY; }

    double getLengthX() const { return m_lengthX; }
    double getLengthY() const { return m_lengthY; }

    vector<vector<double>> getValuesX() const { return m_valuesX; }
    vector<vector<double>> getValuesY() const { return m_valuesY; }

    double getX(size_t _indexX, size_t _indexY) const { return m_valuesX[_indexX][_indexY]; }
    double getY(size_t _indexX, size_t _indexY) const { return m_valuesY[_indexX][_indexY]; }

    double getDX(size_t _indexX, size_t _indexY) const { return m_valuesX[_indexX + 1][_indexY] - m_valuesX[_indexX][_indexY]; }
    double getDY(size_t _indexX, size_t _indexY) const { return m_valuesY[_indexX][_indexY + 1] - m_valuesY[_indexX][_indexY]; }

private:
    // Sizes: X, Y
    size_t m_sizeX;
    size_t m_sizeY;

    // Limits: X, Y
    double m_minX, m_maxX;
    double m_minY, m_maxY;

    // Length: X, Y
    double m_lengthX;
    double m_lengthY;

    // Knots: X, Y
    vector<vector<double>> m_valuesX;
    vector<vector<double>> m_valuesY;
};

// ---------- // 3D // ----------

template <>
class Mesh<3>
{
public:
    // Constructors and destructors
    Mesh(){};
    Mesh(double _minX, double _minY, double _minZ, double _maxX, double _maxY, double _maxZ, size_t _sizeX, size_t _sizeY, size_t _sizeZ /* reference to transform */);
    ~Mesh(){};

    // Getters
    size_t getSizeX() const { return m_sizeX; }
    size_t getSizeY() const { return m_sizeY; }
    size_t getSizeZ() const { return m_sizeZ; }

    double getMinX() const { return m_minX; }
    double getMaxX() const { return m_maxX; }

    double getMinY() const { return m_minY; }
    double getMaxY() const { return m_maxY; }

    double getMinZ() const { return m_minZ; }
    double getMaxZ() const { return m_maxZ; }

    double getLengthX() const { return m_lengthX; }
    double getLengthY() const { return m_lengthY; }
    double getLengthZ() const { return m_lengthZ; }

    vector<vector<vector<double>>> getValuesX() const { return m_valuesX; }
    vector<vector<vector<double>>> getValuesY() const { return m_valuesY; }
    vector<vector<vector<double>>> getValuesZ() const { return m_valuesZ; }

    double getX(size_t _indexX, size_t _indexY, size_t _indexZ) const { return m_valuesX[_indexX][_indexY][_indexZ]; }
    double getY(size_t _indexX, size_t _indexY, size_t _indexZ) const { return m_valuesY[_indexX][_indexY][_indexZ]; }
    double getZ(size_t _indexX, size_t _indexY, size_t _indexZ) const { return m_valuesZ[_indexX][_indexY][_indexZ]; }

    double getDX(size_t _indexX, size_t _indexY, size_t _indexZ) const {return m_valuesX[_indexX + 1][_indexY][_indexZ] - m_valuesX[_indexX][_indexY][_indexZ];}
    double getDY(size_t _indexX, size_t _indexY, size_t _indexZ) const {return m_valuesY[_indexX][_indexY + 1][_indexZ] - m_valuesY[_indexX][_indexY][_indexZ];}
    double getDZ(size_t _indexX, size_t _indexY, size_t _indexZ) const {return m_valuesZ[_indexX][_indexY][_indexZ + 1] - m_valuesZ[_indexX][_indexY][_indexZ];}

private:
    // Sizes: X, Y, Z
    size_t m_sizeX;
    size_t m_sizeY;
    size_t m_sizeZ;

    // Limits: X, Y, Z
    double m_minX, m_maxX;
    double m_minY, m_maxY;
    double m_minZ, m_maxZ;

    // Length: X, Y, Z
    double m_lengthX;
    double m_lengthY;
    double m_lengthZ;

    // Knots: X, Y, Z
    vector<vector<vector<double>>> m_valuesX;
    vector<vector<vector<double>>> m_valuesY;
    vector<vector<vector<double>>> m_valuesZ;
};

#endif // MESH_H
