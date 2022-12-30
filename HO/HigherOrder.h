#if !defined(HIGHERORDER_H)
#define HIGHERORDER_H

#include "Field.h"

class HigherOrder
{
public:
    HigherOrder(/* args */);
    ~HigherOrder();

    void solve();

private:
    // Meshes Order
    double dx;
    double dt;

    size_t
        m_Nx,
        m_Ny,
        m_Nz;

    // Meshes
    Field<double, 1>
        m_meshX,
        m_meshY,
        m_meshZ;

    
    Field<double, 2>
        m_b,
        m_s,
        m_l;

    Field<double, 3>
        m_vX,
        m_vY,
        m_vZ;

    Field <
};

HigherOrder::HigherOrder(/* args */)
{
}

HigherOrder::~HigherOrder()
{
}

#endif // HIGHERORDER_H
