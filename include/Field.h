#if !defined(FIELD_H)
#define FIELD_H

#include <vector>

using std::size_t, std::vector;

template <size_t dimensions>
class Field
{
public:
    Field(){};
    ~Field(){};
};

// 1D //
template <>
class Field<1>
{
public:
    // Constructors and destructors
    Field(){};
    Field(size_t _sizeX);
    ~Field(){};

    // Getters
    size_t getSizeX() const;

    // Operators
    double &operator()(size_t _indexX);

private:
    // Sizes: X
    size_t m_sizeX;
    // Values
    vector<double> m_values;
};

// 2D //
template <>
class Field<2>
{
public:
    // Constructors and destructors
    Field(){};
    Field(size_t _sizeX, size_t _sizeY);
    ~Field(){};

    // Getters
    size_t getSizeX() const;
    size_t getSizeY() const;

    // Operators
    double &operator()(size_t _indexX, size_t _indexY);
    double const operator()(size_t _indexX, size_t _indexY) const;

private:
    // Sizes: X, Y
    size_t m_sizeX;
    size_t m_sizeY;
    // Values
    vector<vector<double>> m_values;
};

// 3D //
template <>
class Field<3>
{
public:
    // Constructors and destructors
    Field(){};
    Field(size_t _sizeX, size_t _sizeY, size_t _sizeZ);
    ~Field(){};

    // Getters
    size_t getSizeX() const;
    size_t getSizeY() const;
    size_t getSizeZ() const;

    // Operators
    double &operator()(size_t _indexX, size_t _indexY, size_t _indexZ);
    double operator()(size_t _indexX, size_t _indexY, size_t _indexZ) const;

private:
    // Sizes: X
    size_t m_sizeX;
    size_t m_sizeY;
    size_t m_sizeZ;
    // Values
    vector<vector<vector<double>>> m_values;
};
#endif // FIELD_H
