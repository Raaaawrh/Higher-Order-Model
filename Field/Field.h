#if !defined(FIELD_H)
#define FIELD_H

#include <cstddef>
#include <vector>

using std::size_t, std::vector;

template <typename T, size_t dim>
class Field
{
public:
private:
    Field(){};
    ~Field(){};
};

template <typename T>
class Field<T, 1>
{
public:
    Field();
    Field(size_t _sizeX);

    ~Field(){};

    T &operator()(size_t _indexX);
    size_t sizeX() const;

private:
    size_t m_Nx;

    vector<T> m_vals;
};

template <typename T>
class Field<T, 2>
{
public:
    Field();
    Field(size_t _sizeX, size_t _sizeY);
    ~Field(){};

    size_t sizeX() const;
    size_t sizeY() const;

    T& operator()(size_t _indexX, size_t _indexY);

private:
    size_t
        m_Nx,
        m_Ny;

    vector<vector<T>> m_vals;
};

template <typename T>
class Field<T, 3>
{
public:
    Field();
    Field(size_t _sizeX, size_t _sizeY, size_t _sizeZ);
    ~Field(){};

    size_t sizeX() const;
    size_t sizeY() const;
    size_t sizeZ() const;

    T& operator()(size_t _indexX, size_t _indexY, size_t _indexZ);

private:
    size_t
        m_Nx,
        m_Ny,
        m_Nz;

    vector<vector<vector<T>>> m_vals;
};

#endif // FIELD_H
