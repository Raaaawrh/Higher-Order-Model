#include "Model.h"

// Инициализация коэффициентов замен переменных

void Model::calc_aX()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                    m_aX(indexX, indexY, indexZ) = -1.0 / m_h(indexX, indexY) *
                                                   (m_dbX(indexX, indexY) +
                                                    m_meshXYZ.getZ(indexX, indexY, indexZ) *
                                                        m_dhX(indexX, indexY));
                else
                    m_aX(indexX, indexY, indexZ) = 0.0;
            }
        }
    }
}

void Model::calc_aY()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                    m_aY(indexX, indexY, indexZ) = -1.0 / m_h(indexX, indexY) *
                                                   (m_dbY(indexX, indexY) +
                                                    m_meshXYZ.getZ(indexX, indexY, indexZ) *
                                                        m_dhY(indexX, indexY));
                else
                    m_aY(indexX, indexY, indexZ) = 0.0;
            }
        }
    }
}

void Model::calc_aZ()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            if (m_h(indexX, indexY) != 0)
                m_aZ(indexX, indexY) = 1.0 / m_h(indexX, indexY);
            else
                m_aZ(indexX, indexY) = 0.0;
        }
    }
}

void Model::calc_bX()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                    m_bX(indexX, indexY, indexZ) = -1.0 / m_h(indexX, indexY) *
                                                   (m_ddbX(indexX, indexY) +
                                                    m_meshXYZ.getZ(indexX, indexY, indexZ) * m_ddhX(indexX, indexY) +
                                                    2 * m_aX(indexX, indexY, indexZ) * m_dhX(indexX, indexY));
                else
                    m_bX(indexX, indexY, indexZ) = 0.0;
            }
        }
    }
}

void Model::calc_bY()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                    m_bY(indexX, indexY, indexZ) = -1.0 / m_h(indexX, indexY) *
                                                   (m_ddbY(indexX, indexY) +
                                                    m_meshXYZ.getZ(indexX, indexY, indexZ) * m_ddhY(indexX, indexY) +
                                                    2 * m_aY(indexX, indexY, indexZ) * m_dhY(indexX, indexY));
                else
                    m_bY(indexX, indexY, indexZ) = 0.0;
            }
        }
    }
}

void Model::calc_cXY()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                    m_cXY(indexX, indexY, indexZ) = -1.0 / m_h(indexX, indexY) *
                                                    (m_ddbXY(indexX, indexY) +
                                                     m_aX(indexX, indexY, indexZ) * m_dhY(indexX, indexY) +
                                                     m_aY(indexX, indexY, indexZ) * m_dhX(indexX, indexY) +
                                                     m_meshXYZ.getZ(indexX, indexY, indexZ) * m_ddhXY(indexX, indexY));
                else
                    m_cXY(indexX, indexY, indexZ) = 0.0;
            }
        }
    }
}

void Model::calc_cXZ()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            if (m_h(indexX, indexY) != 0)
                m_cXZ(indexX, indexY) = -1.0 / m_h(indexX, indexY) / m_h(indexX, indexY) * m_dhX(indexX, indexY);
            else
                m_cXZ(indexX, indexY) = 0.0;
        }
    }
}

void Model::calc_cYZ()
{
    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            if (m_h(indexX, indexY) != 0)
                m_cYZ(indexX, indexY) = -1.0 / m_h(indexX, indexY) / m_h(indexX, indexY) * m_dhY(indexX, indexY);
            else
                m_cYZ(indexX, indexY) = 0.0;
        }
    }
}

// // Нужно вычислить все производные s, b, h
void Model::calc_coeff()
{
    calc_aX();
    calc_aY();
    calc_aZ();

    calc_bX();
    calc_bY();

    calc_cXY();
    calc_cXZ();
    calc_cYZ();
}