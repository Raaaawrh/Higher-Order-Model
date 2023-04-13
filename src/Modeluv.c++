#include "Model.h"

#include <iostream>
#include <cmath>
#include <fstream>

// Инициализация методов для обработки u, v
// А именно, производные для определения вязкости m_eta

/// @brief Вычисление du/dx, du/dy, du/dz
void Model::calc_uDerivatives()
{
    calc_dX(m_u, m_dux);
    calc_dY(m_u, m_duy);
    calc_dZ(m_u, m_duz);
}
/// @brief Вычисление dv/dx, dv/dy, dv/dz
void Model::calc_vDerivatives()
{
    calc_dX(m_v, m_dvx);
    calc_dY(m_v, m_dvy);
    calc_dZ(m_v, m_dvz);
}

// // Заполнение системы для решения уравнений на скорость

/// @brief Заполнение системы для u
void Model::fill_uSystem()
{
    m_uvMatrix.clear();

    // Границы по x
    for (size_t indexY{0}; indexY < m_sizeY; indexY++)
    {
        for (size_t indexZ{0}; indexY < m_sizeY; indexY++)
        {
            // x = minX
            // m_uvMatrix.coeffRef(getIndexXYZ(0, indexY, indexZ), getIndexXYZ(0, indexY, indexZ)) = 1;
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(0, indexY, indexZ),
                getIndexXYZ(0, indexY, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(0, indexY, indexZ)) = 0;

            // x = maxX
            // m_uvMatrix.coeffRef(getIndexXYZ(m_sizeX - 1, indexY, indexZ), getIndexXYZ(m_sizeX - 1, indexY, indexZ)) = 1;

            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(m_sizeX - 1, indexY, indexZ)) = 0;
        }
    }

    // Границы по y
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
        {
            // y = minY
            // m_uvMatrix.coeffRef(getIndexXYZ(indexX, 0, indexZ), getIndexXYZ(indexX, 0, indexZ)) = 1;
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, 0, indexZ),
                getIndexXYZ(indexX, 0, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(indexX, 0, indexZ)) = 0;

            // y = maxY
            // m_uvMatrix.coeffRef(getIndexXYZ(indexX, m_sizeY - 1, indexZ), getIndexXYZ(indexX, m_sizeY - 1, indexZ)) = 1;
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(indexX, m_sizeY - 1, indexZ)) = 0;
        }
    }

    // Границы по z

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            // z = 0
            // m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, 0), getIndexXYZ(indexX, indexY, 0)) = 1;
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, indexY, 0),
                getIndexXYZ(indexX, indexY, 0),
                1));

            m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, 0)) = 0;

            // z = 1 (free-stress surface)
            if (m_h(indexX, indexY) != 0)
            {
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX - 1, indexY, m_sizeZ - 1)) = -2.0 * m_dsX(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY - 1, m_sizeZ - 1)) = -0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY, m_sizeZ - 2)) = (-4 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) - m_aZ(indexX, indexY)) /
                //                                                                                                                          deltaZ(indexX, indexY, m_sizeZ - 2, 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (4 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) + m_aZ(indexX, indexY)) /
                //                                                                                                                          deltaZ(indexX, indexY, m_sizeZ - 2, 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY + 1, m_sizeZ - 1)) = 0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX + 1, indexY, m_sizeZ - 1)) = 2.0 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1);
                //
                //                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_v(indexX - 1, indexY, m_sizeZ - 1) +
                //                                                                                (m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_v(indexX, indexY - 1, m_sizeZ - 1) +
                //                                                                                (2 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_v(indexX, indexY, m_sizeZ - 2) +
                //                                                                                (-2 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_v(indexX, indexY, m_sizeZ - 1) +
                //                                                                                (-m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_v(indexX, indexY + 1, m_sizeZ - 1) +
                //                                                                                (-0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_v(indexX + 1, indexY, m_sizeZ - 1);

                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX - 1, indexY, m_sizeZ - 1),
                    -2.0 * m_dsX(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY - 1, m_sizeZ - 1),
                    -0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 2),
                    (-4 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) - m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    (4 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) + m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY + 1, m_sizeZ - 1),
                    0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX + 1, indexY, m_sizeZ - 1),
                    2.0 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));

                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_v(indexX - 1, indexY, m_sizeZ - 1) +
                                                                                (m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_v(indexX, indexY - 1, m_sizeZ - 1) +
                                                                                (2 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_v(indexX, indexY, m_sizeZ - 2) +
                                                                                (-2 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_v(indexX, indexY, m_sizeZ - 1) +
                                                                                (-m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_v(indexX, indexY + 1, m_sizeZ - 1) +
                                                                                (-0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_v(indexX + 1, indexY, m_sizeZ - 1);
            }
            else // Если нет вещества
            {
                // m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 1;
                // m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 0;

                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    1));
                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 0;
            }
        }
    }

    // Внутренние узлы

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{1}; indexZ < m_sizeZ - 1; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                {
                    // Левая часть
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX - 1, indexY, indexZ - 1)) = 8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX - 1, indexY, indexZ)) = 8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + 4 * betaXX(indexX, indexY, indexZ, -1, 0, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX - 1, indexY, indexZ + 1)) = 8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY - 1, indexZ - 1)) = 2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY - 1, indexZ)) = 2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + betaYY(indexX, indexY, indexZ, 0, -1, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY - 1, indexZ + 1)) = 2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ - 1)) = (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) +
                    //                                                                                                                        (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ)) = (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) +
                    //                                                                                                                    (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 0) +
                    //                                                                                                                    4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) +
                    //                                                                                                                    m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) +
                    //                                                                                                                    4 * betaXX(indexX, indexY, indexZ, 0, 0, 0) +
                    //                                                                                                                    betaYY(indexX, indexY, indexZ, 0, 0, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ + 1)) = (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) +
                    //                                                                                                                        (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY + 1, indexZ - 1)) = 2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY + 1, indexZ)) = 2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + betaYY(indexX, indexY, indexZ, 0, 1, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY + 1, indexZ + 1)) = 2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX + 1, indexY, indexZ - 1)) = 8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX + 1, indexY, indexZ)) = 8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + 4 * betaXX(indexX, indexY, indexZ, 1, 0, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX + 1, indexY, indexZ + 1)) = 8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1);

                    //

                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ - 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + 4 * betaXX(indexX, indexY, indexZ, -1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ + 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ - 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + betaYY(indexX, indexY, indexZ, 0, -1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ + 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ - 1),
                        (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) + (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ),
                        (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) + (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) + 4 * betaXX(indexX, indexY, indexZ, 0, 0, 0) + betaYY(indexX, indexY, indexZ, 0, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ + 1),
                        (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) + (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ - 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + betaYY(indexX, indexY, indexZ, 0, 1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ + 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ - 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + 4 * betaXX(indexX, indexY, indexZ, 1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ + 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)));

                    // Правая часть
                    m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, indexZ)) = (-3 * betaXY(indexX, indexY, indexZ, -1, -1, 0)) * m_v(indexX - 1, indexY - 1, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)) * m_v(indexX - 1, indexY, indexZ - 1) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) - m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) - betaXY(indexX, indexY, indexZ, -1, 0, 0)) * m_v(indexX - 1, indexY, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)) * m_v(indexX - 1, indexY, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, -1, 1, 0)) * m_v(indexX - 1, indexY + 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)) * m_v(indexX, indexY - 1, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) - 2 * m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) - betaXY(indexX, indexY, indexZ, 0, -1, 0)) * m_v(indexX, indexY - 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)) * m_v(indexX, indexY - 1, indexZ + 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, -1) - m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, -1) + m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, -1)) * m_v(indexX, indexY, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 0) - m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) - m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 0) - 2 * m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) + m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * betaXY(indexX, indexY, indexZ, 0, 0, 0)) * m_v(indexX, indexY, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 1) - m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 1) + m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 1)) * m_v(indexX, indexY, indexZ + 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)) * m_v(indexX, indexY + 1, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) - 2 * m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) - betaXY(indexX, indexY, indexZ, 0, 1, 0)) * m_v(indexX, indexY + 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)) * m_v(indexX, indexY + 1, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, 1, -1, 0)) * m_v(indexX + 1, indexY - 1, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)) * m_v(indexX + 1, indexY, indexZ - 1) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) - m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) - betaXY(indexX, indexY, indexZ, 1, 0, 0)) * m_v(indexX + 1, indexY, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)) * m_v(indexX + 1, indexY, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, 1, 1, 0)) * m_v(indexX + 1, indexY + 1, indexZ) +
                                                                               (rho * g * m_dsX(indexX, indexY));
                }
                else
                {
                    // m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ)) = 1;

                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ),
                        1));
                    m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, indexZ)) = 0;
                }
            }
        }
    }

    // std::cout << m_uvMatrix;
    // std::cout << m_uvVector;
}
/// @brief Заполнение системы для v
void Model::fill_vSystem()
{
    // m_uvMatrix.setZero();
    m_uvMatrix.clear();

    // Границы по x
    for (size_t indexY{0}; indexY < m_sizeY; indexY++)
    {
        for (size_t indexZ{0}; indexY < m_sizeY; indexY++)
        {
            // x = minX
            // m_uvMatrix.coeffRef(getIndexXYZ(0, indexY, indexZ), getIndexXYZ(0, indexY, indexZ)) = 1;

            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(0, indexY, indexZ),
                getIndexXYZ(0, indexY, indexZ),
                1));
            m_uvVector.coeffRef(getIndexXYZ(0, indexY, indexZ)) = 0;

            // x = maxX
            // m_uvMatrix.coeffRef(getIndexXYZ(m_sizeX - 1, indexY, indexZ), getIndexXYZ(m_sizeX - 1, indexY, indexZ)) = 1;

            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                1));
            m_uvVector.coeffRef(getIndexXYZ(m_sizeX - 1, indexY, indexZ)) = 0;
        }
    }

    // Границы по y
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
        {
            // y = minY
            // m_uvMatrix.coeffRef(getIndexXYZ(indexX, 0, indexZ), getIndexXYZ(indexX, 0, indexZ)) = 1;

            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, 0, indexZ),
                getIndexXYZ(indexX, 0, indexZ),
                1));
            m_uvVector.coeffRef(getIndexXYZ(indexX, 0, indexZ)) = 0;

            // y = maxY
            // m_uvMatrix.coeffRef(getIndexXYZ(indexX, m_sizeY - 1, indexZ), getIndexXYZ(indexX, m_sizeY - 1, indexZ)) = 1;

            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                1));
            m_uvVector.coeffRef(getIndexXYZ(indexX, m_sizeY - 1, indexZ)) = 0;
        }
    }

    // Границы по z

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            // z = 0
            // m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, 0), getIndexXYZ(indexX, indexY, 0)) = 1;

            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, indexY, 0),
                getIndexXYZ(indexX, indexY, 0),
                1));
            m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, 0)) = 0;

            // z = 1 (free-stress surface)
            if (m_h(indexX, indexY) != 0)
            {
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX - 1, indexY, m_sizeZ - 1)) = -0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY - 1, m_sizeZ - 1)) = -2.0 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY, m_sizeZ - 2)) = (-m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - 4 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) - m_aZ(indexX, indexY)) /
                //                                                                                                                          deltaZ(indexX, indexY, m_sizeZ - 2, 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + 4 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) + m_aZ(indexX, indexY)) /
                //                                                                                                                          deltaZ(indexX, indexY, m_sizeZ - 2, 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY + 1, m_sizeZ - 1)) = 2.0 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1);
                //                m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX + 1, indexY, m_sizeZ - 1)) = 0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1);

                //

                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX - 1, indexY, m_sizeZ - 1),
                    -0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY - 1, m_sizeZ - 1),
                    -2.0 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 2),
                    (-m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - 4 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) - m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    (m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + 4 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) + m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY + 1, m_sizeZ - 1),
                    2.0 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX + 1, indexY, m_sizeZ - 1),
                    0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));

                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_u(indexX - 1, indexY, m_sizeZ - 1) +
                                                                                (0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_u(indexX, indexY - 1, m_sizeZ - 1) +
                                                                                (m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + 2 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_u(indexX, indexY, m_sizeZ - 2) +
                                                                                (-m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - 2 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_u(indexX, indexY, m_sizeZ - 1) +
                                                                                (-0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_u(indexX, indexY + 1, m_sizeZ - 1) +
                                                                                (-m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_u(indexX + 1, indexY, m_sizeZ - 1);
            }
            else // Если нет вещества
            {
                // m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1), getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 1;

                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    1));
                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 0;
            }
        }
    }

    // Внутренние узлы

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{1}; indexZ < m_sizeZ - 1; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                {
                    // Левая часть
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX - 1, indexY, indexZ - 1)) = 2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX - 1, indexY, indexZ)) = 2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + betaXX(indexX, indexY, indexZ, -1, 0, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX - 1, indexY, indexZ + 1)) = 2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY - 1, indexZ - 1)) = 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY - 1, indexZ)) = 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, -1, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY - 1, indexZ + 1)) = 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ - 1)) = (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) +
                    //                                                                                                                        (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ)) = (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) +
                    //                                                                                                                    (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 0) +
                    //                                                                                                                    m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) +
                    //                                                                                                                    4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) +
                    //                                                                                                                    betaXX(indexX, indexY, indexZ, 0, 0, 0) +
                    //                                                                                                                    4 * betaYY(indexX, indexY, indexZ, 0, 0, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ + 1)) = (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) +
                    //                                                                                                                        (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY + 1, indexZ - 1)) = 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY + 1, indexZ)) = 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, 1, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY + 1, indexZ + 1)) = 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX + 1, indexY, indexZ - 1)) = 2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX + 1, indexY, indexZ)) = 2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + betaXX(indexX, indexY, indexZ, 1, 0, 0);
                    //                    m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX + 1, indexY, indexZ + 1)) = 2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1);

                    //

                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ - 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + betaXX(indexX, indexY, indexZ, -1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ + 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ - 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, -1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ + 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ - 1),
                        (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) + (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ),
                        (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) + (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) + betaXX(indexX, indexY, indexZ, 0, 0, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ + 1),
                        (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) + (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ - 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ), 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, 1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ + 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ - 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + betaXX(indexX, indexY, indexZ, 1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ + 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)));

                    // Правая часть
                    m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, indexZ)) = (-3 * betaXY(indexX, indexY, indexZ, -1, -1, 0)) * m_u(indexX - 1, indexY - 1, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)) * m_u(indexX - 1, indexY, indexZ - 1) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) - 2 * m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + betaXY(indexX, indexY, indexZ, -1, 0, 0)) * m_u(indexX - 1, indexY, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)) * m_u(indexX - 1, indexY, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, -1, 1, 0)) * m_u(indexX - 1, indexY + 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)) * m_u(indexX, indexY - 1, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) - m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + betaXY(indexX, indexY, indexZ, 0, -1, 0)) * m_u(indexX, indexY - 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)) * m_u(indexX, indexY - 1, indexZ + 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, -1) + m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, -1) - m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, -1)) * m_u(indexX, indexY, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 0) - 2 * m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) + m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 0) - m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) - m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * betaXY(indexX, indexY, indexZ, 0, 0, 0)) * m_u(indexX, indexY, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 1) + m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 1) - m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 1)) * m_u(indexX, indexY, indexZ + 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)) * m_u(indexX, indexY + 1, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) - m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + betaXY(indexX, indexY, indexZ, 0, 1, 0)) * m_u(indexX, indexY + 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)) * m_u(indexX, indexY + 1, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, 1, -1, 0)) * m_u(indexX + 1, indexY - 1, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)) * m_u(indexX + 1, indexY, indexZ - 1) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) - 2 * m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + betaXY(indexX, indexY, indexZ, 1, 0, 0)) * m_u(indexX + 1, indexY, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)) * m_u(indexX + 1, indexY, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, 1, 1, 0)) * m_u(indexX + 1, indexY + 1, indexZ) +
                                                                               (rho * g * m_dsY(indexX, indexY));
                }
                else
                {
                    // m_uvMatrix.coeffRef(getIndexXYZ(indexX, indexY, indexZ), getIndexXYZ(indexX, indexY, indexZ)) = 1;

                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ),
                        1));
                    m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, indexZ)) = 0;
                }
            }
        }
    }

    // std::cout << m_uvMatrix;
    // std::cout << m_uvVector;
}

// // Решение системы для скорости

/// @brief Решение системы для u. Решение -- m_uIter
void Model::solve_u()
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.setTolerance(1e-4);
    // solver.setMaxIterations(1000);
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(m_nXYZ, m_nXYZ);
    matrix.setFromTriplets(m_uvMatrix.begin(), m_uvMatrix.end());

    Eigen::Vector<double, Eigen::Dynamic> guess(m_nXYZ);

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                guess(getIndexXYZ(indexX, indexY, indexZ)) = m_u(indexX, indexY, indexZ);
            }
        }
    }

    solver.compute(matrix);
    Eigen::Vector<double, Eigen::Dynamic> solution = solver.solveWithGuess(m_uvVector, guess);

    m_alpha_u = 0.0;

    double alpha1{0.0};
    double alpha2{0.0};

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                alpha1 += m_u(indexX, indexY, indexZ) * m_u(indexX, indexY, indexZ);
                alpha2 += (solution(getIndexXYZ(indexX, indexY, indexZ)) - m_u(indexX, indexY, indexZ)) *
                          (solution(getIndexXYZ(indexX, indexY, indexZ)) - m_u(indexX, indexY, indexZ));
            }
        }
    }

    alpha1 = std::sqrt(alpha1);
    alpha2 = std::sqrt(alpha2);

    if (alpha2 != 0.0)
        m_alpha_u = alpha1 / alpha2;

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_uIter(indexX, indexY, indexZ) = solution(getIndexXYZ(indexX, indexY, indexZ));
            }
        }
    }

    std::cout << "U: " << solver.iterations() << std::endl;

    // std::cout << "Solution" << std::endl;
    // std::cout << solution;
}
/// @brief Решение системы для v. Решение -- m_vIter
void Model::solve_v()
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    // solver.setTolerance(1e-10);
    // solver.setMaxIterations(1000);
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(m_nXYZ, m_nXYZ);
    matrix.setFromTriplets(m_uvMatrix.begin(), m_uvMatrix.end());

    Eigen::Vector<double, Eigen::Dynamic> guess(m_nXYZ);

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                guess(getIndexXYZ(indexX, indexY, indexZ)) = m_v(indexX, indexY, indexZ);
            }
        }
    }

    solver.compute(matrix);
    Eigen::Vector<double, Eigen::Dynamic> solution = solver.solveWithGuess(m_uvVector, guess);

    m_alpha_v = 0.0;

    double alpha1{0.0};
    double alpha2{0.0};

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                alpha1 += m_v(indexX, indexY, indexZ) * m_v(indexX, indexY, indexZ);
                alpha2 += (solution(getIndexXYZ(indexX, indexY, indexZ)) - m_v(indexX, indexY, indexZ)) *
                          (solution(getIndexXYZ(indexX, indexY, indexZ)) - m_v(indexX, indexY, indexZ));
            }
        }
    }

    alpha1 = std::sqrt(alpha1);
    alpha2 = std::sqrt(alpha2);

    m_alpha_v = alpha1 / alpha2;

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_vIter(indexX, indexY, indexZ) = solution(getIndexXYZ(indexX, indexY, indexZ));
            }
        }
    }

    std::cout << "V: " << solver.iterations() << std::endl;

    // std::cout << "Solution" << std::endl;
    // std::cout << solution;
}

// // Вычисление отклонений скоростей на итерациях

/// @brief Вычисление отклонения u на итерациях
void Model::calc_uErr()
{
    m_uErr = std::numeric_limits<double>::lowest();

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_uErr = std::max(std::fabs(m_u(indexX, indexY, indexZ) - m_uIter(indexX, indexY, indexZ)), m_uErr);
            }
        }
    }
}

/// @brief Вычисление отклонения v на итерациях
void Model::calc_vErr()
{
    m_vErr = std::numeric_limits<double>::lowest();

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_vErr = std::max(std::fabs(m_v(indexX, indexY, indexZ) - m_vIter(indexX, indexY, indexZ)), m_vErr);
            }
        }
    }
}

// // Вычисление полей средней скорости

/// @brief Вычисление средней горизонтальной скорости u
void Model::calc_uAverage()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_uAverage(indexX, indexY) = 0.0;
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_uAverage(indexX, indexY) += m_u(indexX, indexY, indexZ);
            }
            m_uAverage(indexX, indexY) /= m_sizeZ;
        }
    }
}
/// @brief Вычисление средней горизонтальной скорости v
void Model::calc_vAverage()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_vAverage(indexX, indexY) = 0.0;
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_vAverage(indexX, indexY) += m_v(indexX, indexY, indexZ);
            }
            m_vAverage(indexX, indexY) /= m_sizeZ;
        }
    }
}

// // Запись полученной скорости в предыдущую итерацию

/// @brief Обновление скорости u на итерации
void Model::update_uIter()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_u(indexX, indexY, indexZ) = 0.5 * m_uIter(indexX, indexY, indexZ) +
                                              (1 - 0.5) * m_u(indexX, indexY, indexZ);
            }
        }
    }
}
/// @brief Обновление скорости v на итерации
void Model::update_vIter()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_v(indexX, indexY, indexZ) = 0.5 * m_vIter(indexX, indexY, indexZ) +
                                              (1 - 0.5) * m_v(indexX, indexY, indexZ);
            }
        }
    }
}

// // Вязкость

/// @brief Вычисление вязкости в центре ячейки с индексом (_indexX + _dirX/2, _indexY + _dirY/2, _indexZ + _dirZ/2)
double Model::getEtaoh(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};

    result = m_eta(_indexX, _indexY, _indexZ) +
             m_eta(_indexX, _indexY, _indexZ + _dirZ) +
             m_eta(_indexX, _indexY + _dirY, _indexZ) +
             m_eta(_indexX, _indexY + _dirY, _indexZ + _dirZ) +
             m_eta(_indexX + _dirX, _indexY, _indexZ) +
             m_eta(_indexX + _dirX, _indexY, _indexZ + _dirZ) +
             m_eta(_indexX + _dirX, _indexY + _dirY, _indexZ) +
             m_eta(_indexX + _dirX, _indexY + _dirY, _indexZ + _dirZ);
    return result / 8.0;
}

/// @brief Вычисление поля вязкости в каждой точке
void Model::calc_eta()
{

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_eta(indexX, indexY, indexZ) = 0.5 * std::pow(A(), -1.0 / n) *
                                                std::pow((m_dux(indexX, indexY, indexZ) + m_aX(indexX, indexY, indexZ) * m_duz(indexX, indexY, indexZ)) *
                                                                 (m_dux(indexX, indexY, indexZ) + m_aX(indexX, indexY, indexZ) * m_duz(indexX, indexY, indexZ)) +
                                                             (m_dvy(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_dvz(indexX, indexY, indexZ)) *
                                                                 (m_dvy(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_dvz(indexX, indexY, indexZ)) +
                                                             (m_dux(indexX, indexY, indexZ) + m_aX(indexX, indexY, indexZ) * m_duz(indexX, indexY, indexZ)) *
                                                                 (m_dvy(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_dvz(indexX, indexY, indexZ)) +
                                                             0.25 *
                                                                 (m_duy(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_duz(indexX, indexY, indexZ) + m_dvx(indexX, indexY, indexZ) + m_aX(indexX, indexY, indexZ) * m_dvz(indexX, indexY, indexZ)) *
                                                                 (m_duy(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_duz(indexX, indexY, indexZ) + m_dvx(indexX, indexY, indexZ) + m_aX(indexX, indexY, indexZ) * m_dvz(indexX, indexY, indexZ)) +
                                                             0.25 *
                                                                 (m_aZ(indexX, indexY) * m_duz(indexX, indexY, indexZ)) *
                                                                 (m_aZ(indexX, indexY) * m_duz(indexX, indexY, indexZ)) +
                                                             0.25 *
                                                                 (m_aZ(indexX, indexY) * m_dvz(indexX, indexY, indexZ)) *
                                                                 (m_aZ(indexX, indexY) * m_dvz(indexX, indexY, indexZ)) +
                                                             m_eta0,
                                                         (1.0 - n) / (2.0 * n));
            }
        }
    }
}

void Model::fill_uvSystem()
{
    ////// Часть для u
    m_uvMatrix.clear();
    m_uvVector.resize(2 * m_nXYZ);

    // Границы по x
    for (size_t indexY{0}; indexY < m_sizeY; indexY++)
    {
        for (size_t indexZ{0}; indexY < m_sizeY; indexY++)
        {
            // x = minX
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(0, indexY, indexZ),
                getIndexXYZ(0, indexY, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(0, indexY, indexZ)) = 0;

            // x = maxX
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(m_sizeX - 1, indexY, indexZ)) = 0;
        }
    }

    // Границы по y
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
        {
            // y = minY
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, 0, indexZ),
                getIndexXYZ(indexX, 0, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(indexX, 0, indexZ)) = 0;

            // y = maxY
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                1));

            m_uvVector.coeffRef(getIndexXYZ(indexX, m_sizeY - 1, indexZ)) = 0;
        }
    }

    // Границы по z

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            // z = 0
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                getIndexXYZ(indexX, indexY, 0),
                getIndexXYZ(indexX, indexY, 0),
                1));

            m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, 0)) = 0;

            // z = 1 (free-stress surface)
            if (m_h(indexX, indexY) != 0)
            {
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX - 1, indexY, m_sizeZ - 1),
                    -2.0 * m_dsX(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY - 1, m_sizeZ - 1),
                    -0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 2),
                    (-4 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) - m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    (4 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) + m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY + 1, m_sizeZ - 1),
                    0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX + 1, indexY, m_sizeZ - 1),
                    2.0 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));

                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_v(indexX - 1, indexY, m_sizeZ - 1) +
                                                                                (m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_v(indexX, indexY - 1, m_sizeZ - 1) +
                                                                                (2 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_v(indexX, indexY, m_sizeZ - 2) +
                                                                                (-2 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_v(indexX, indexY, m_sizeZ - 1) +
                                                                                (-m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_v(indexX, indexY + 1, m_sizeZ - 1) +
                                                                                (-0.5 * m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_v(indexX + 1, indexY, m_sizeZ - 1);
            }
            else // Если нет вещества
            {
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    1));
                m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 0;
            }
        }
    }

    // Внутренние узлы

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{1}; indexZ < m_sizeZ - 1; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                {
                    // Левая часть
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ - 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + 4 * betaXX(indexX, indexY, indexZ, -1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX - 1, indexY, indexZ + 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ - 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + betaYY(indexX, indexY, indexZ, 0, -1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY - 1, indexZ + 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ - 1),
                        (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) + (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ),
                        (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) + (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) + 4 * betaXX(indexX, indexY, indexZ, 0, 0, 0) + betaYY(indexX, indexY, indexZ, 0, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ + 1),
                        (4 * m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) + (4 * m_bX(indexX, indexY, indexZ) + m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ - 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) + m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + betaYY(indexX, indexY, indexZ, 0, 1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY + 1, indexZ + 1),
                        2 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ - 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) + 4 * m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + 4 * betaXX(indexX, indexY, indexZ, 1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX + 1, indexY, indexZ + 1),
                        8 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)));

                    // Правая часть
                    m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, indexZ)) = (-3 * betaXY(indexX, indexY, indexZ, -1, -1, 0)) * m_v(indexX - 1, indexY - 1, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)) * m_v(indexX - 1, indexY, indexZ - 1) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) - m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) - betaXY(indexX, indexY, indexZ, -1, 0, 0)) * m_v(indexX - 1, indexY, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)) * m_v(indexX - 1, indexY, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, -1, 1, 0)) * m_v(indexX - 1, indexY + 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)) * m_v(indexX, indexY - 1, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) - 2 * m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) - betaXY(indexX, indexY, indexZ, 0, -1, 0)) * m_v(indexX, indexY - 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)) * m_v(indexX, indexY - 1, indexZ + 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, -1) - m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, -1) + m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, -1)) * m_v(indexX, indexY, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 0) - m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) - m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 0) - 2 * m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) + m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * betaXY(indexX, indexY, indexZ, 0, 0, 0)) * m_v(indexX, indexY, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 1) - m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 1) + m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 1)) * m_v(indexX, indexY, indexZ + 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)) * m_v(indexX, indexY + 1, indexZ - 1) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) - 2 * m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) - betaXY(indexX, indexY, indexZ, 0, 1, 0)) * m_v(indexX, indexY + 1, indexZ) +
                                                                               (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)) * m_v(indexX, indexY + 1, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, 1, -1, 0)) * m_v(indexX + 1, indexY - 1, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)) * m_v(indexX + 1, indexY, indexZ - 1) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) - m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) - betaXY(indexX, indexY, indexZ, 1, 0, 0)) * m_v(indexX + 1, indexY, indexZ) +
                                                                               (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)) * m_v(indexX + 1, indexY, indexZ + 1) +
                                                                               (-3 * betaXY(indexX, indexY, indexZ, 1, 1, 0)) * m_v(indexX + 1, indexY + 1, indexZ) +
                                                                               (rho * g * m_dsX(indexX, indexY));
                }
                else
                {
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        getIndexXYZ(indexX, indexY, indexZ),
                        getIndexXYZ(indexX, indexY, indexZ),
                        1));
                    m_uvVector.coeffRef(getIndexXYZ(indexX, indexY, indexZ)) = 0;
                }
            }
        }
    }

    /////// Часть для v

    // Границы по x
    for (size_t indexY{0}; indexY < m_sizeY; indexY++)
    {
        for (size_t indexZ{0}; indexY < m_sizeY; indexY++)
        {
            // x = minX
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                m_nXYZ + getIndexXYZ(0, indexY, indexZ),
                m_nXYZ + getIndexXYZ(0, indexY, indexZ),
                1));
            m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(0, indexY, indexZ)) = 0;

            // x = maxX
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                m_nXYZ + getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                m_nXYZ + getIndexXYZ(m_sizeX - 1, indexY, indexZ),
                1));
            m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(m_sizeX - 1, indexY, indexZ)) = 0;
        }
    }

    // Границы по y
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
        {
            // y = minY
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                m_nXYZ + getIndexXYZ(indexX, 0, indexZ),
                m_nXYZ + getIndexXYZ(indexX, 0, indexZ),
                1));
            m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, 0, indexZ)) = 0;

            // y = maxY
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                m_nXYZ + getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                m_nXYZ + getIndexXYZ(indexX, m_sizeY - 1, indexZ),
                1));
            m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, m_sizeY - 1, indexZ)) = 0;
        }
    }

    // Границы по z

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            // z = 0
            m_uvMatrix.push_back(Eigen::Triplet<double>(
                m_nXYZ + getIndexXYZ(indexX, indexY, 0),
                m_nXYZ + getIndexXYZ(indexX, indexY, 0),
                1));
            m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, indexY, 0)) = 0;

            // z = 1 (free-stress surface)
            if (m_h(indexX, indexY) != 0)
            {
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX - 1, indexY, m_sizeZ - 1),
                    -0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX, indexY - 1, m_sizeZ - 1),
                    -2.0 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 2),
                    (-m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - 4 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) - m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    (m_aX(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + 4 * m_aY(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY) + m_aZ(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX, indexY + 1, m_sizeZ - 1),
                    2.0 * m_dsY(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX + 1, indexY, m_sizeZ - 1),
                    0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)));

                m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = (m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_u(indexX - 1, indexY, m_sizeZ - 1) +
                                                                                         (0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_u(indexX, indexY - 1, m_sizeZ - 1) +
                                                                                         (m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) + 2 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_u(indexX, indexY, m_sizeZ - 2) +
                                                                                         (-m_aY(indexX, indexY, m_sizeZ - 1) * m_dsX(indexX, indexY) - 2 * m_aX(indexX, indexY, m_sizeZ - 1) * m_dsY(indexX, indexY)) / deltaZ(indexX, indexY, m_sizeZ - 2, 1) * m_u(indexX, indexY, m_sizeZ - 1) +
                                                                                         (-0.5 * m_dsX(indexX, indexY) / m_meshXYZ.getDY(indexX, indexY, m_sizeZ - 1)) * m_u(indexX, indexY + 1, m_sizeZ - 1) +
                                                                                         (-m_dsY(indexX, indexY) / m_meshXYZ.getDX(indexX, indexY, m_sizeZ - 1)) * m_u(indexX + 1, indexY, m_sizeZ - 1);
            }
            else // Если нет вещества
            {
                m_uvMatrix.push_back(Eigen::Triplet<double>(
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1),
                    1));
                m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, indexY, m_sizeZ - 1)) = 0;
            }
        }
    }

    // Внутренние узлы

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            for (size_t indexZ{1}; indexZ < m_sizeZ - 1; indexZ++)
            {
                if (m_h(indexX, indexY) != 0)
                {
                    // Левая часть

                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX - 1, indexY, indexZ - 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX - 1, indexY, indexZ),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + betaXX(indexX, indexY, indexZ, -1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX - 1, indexY, indexZ + 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY - 1, indexZ - 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY - 1, indexZ),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, -1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY - 1, indexZ + 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ - 1),
                        (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) + (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) + (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) + betaXX(indexX, indexY, indexZ, 0, 0, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ + 1),
                        (m_aX(indexX, indexY, indexZ) * m_aX(indexX, indexY, indexZ) + 4 * m_aY(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) + m_aZ(indexX, indexY) * m_aZ(indexX, indexY)) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) + (m_bX(indexX, indexY, indexZ) + 4 * m_bY(indexX, indexY, indexZ)) * betaZ(indexX, indexY, indexZ, 0, 0, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY + 1, indexZ - 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY + 1, indexZ), 8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) + 4 * m_aY(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + 4 * betaYY(indexX, indexY, indexZ, 0, 1, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY + 1, indexZ + 1),
                        8 * m_aY(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX + 1, indexY, indexZ - 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX + 1, indexY, indexZ),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) + m_aX(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + betaXX(indexX, indexY, indexZ, 1, 0, 0)));
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX + 1, indexY, indexZ + 1),
                        2 * m_aX(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)));

                    // Правая часть
                    m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, indexY, indexZ)) = (-3 * betaXY(indexX, indexY, indexZ, -1, -1, 0)) * m_u(indexX - 1, indexY - 1, indexZ) +
                                                                                        (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, -1)) * m_u(indexX - 1, indexY, indexZ - 1) +
                                                                                        (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 0) - 2 * m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, -1, 0) + betaXY(indexX, indexY, indexZ, -1, 0, 0)) * m_u(indexX - 1, indexY, indexZ) +
                                                                                        (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, -1, 0, 1)) * m_u(indexX - 1, indexY, indexZ + 1) +
                                                                                        (-3 * betaXY(indexX, indexY, indexZ, -1, 1, 0)) * m_u(indexX - 1, indexY + 1, indexZ) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, -1)) * m_u(indexX, indexY - 1, indexZ - 1) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 0) - m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, -1, 0) + betaXY(indexX, indexY, indexZ, 0, -1, 0)) * m_u(indexX, indexY - 1, indexZ) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, -1, 1)) * m_u(indexX, indexY - 1, indexZ + 1) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, -1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, -1) + m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, -1) - m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, -1)) * m_u(indexX, indexY, indexZ - 1) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 0) - 2 * m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 0, 0) + m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 0) - m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 0, 0) - m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 0) - 3 * betaXY(indexX, indexY, indexZ, 0, 0, 0)) * m_u(indexX, indexY, indexZ) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * m_aY(indexX, indexY, indexZ) * betaZZ(indexX, indexY, indexZ, 0, 0, 1) - 3 * m_cXY(indexX, indexY, indexZ) * betaZ(indexX, indexY, indexZ, 0, 0, 1) + m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 0, 0, 1) - m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 0, 1)) * m_u(indexX, indexY, indexZ + 1) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, -1)) * m_u(indexX, indexY + 1, indexZ - 1) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 0) - m_aX(indexX, indexY, indexZ) * gammaZY(indexX, indexY, indexZ, 1, 0) + betaXY(indexX, indexY, indexZ, 0, 1, 0)) * m_u(indexX, indexY + 1, indexZ) +
                                                                                        (-3 * m_aX(indexX, indexY, indexZ) * betaYZ(indexX, indexY, indexZ, 0, 1, 1)) * m_u(indexX, indexY + 1, indexZ + 1) +
                                                                                        (-3 * betaXY(indexX, indexY, indexZ, 1, -1, 0)) * m_u(indexX + 1, indexY - 1, indexZ) +
                                                                                        (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, -1)) * m_u(indexX + 1, indexY, indexZ - 1) +
                                                                                        (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 0) - 2 * m_aY(indexX, indexY, indexZ) * gammaZX(indexX, indexY, indexZ, 1, 0) + betaXY(indexX, indexY, indexZ, 1, 0, 0)) * m_u(indexX + 1, indexY, indexZ) +
                                                                                        (-3 * m_aY(indexX, indexY, indexZ) * betaXZ(indexX, indexY, indexZ, 1, 0, 1)) * m_u(indexX + 1, indexY, indexZ + 1) +
                                                                                        (-3 * betaXY(indexX, indexY, indexZ, 1, 1, 0)) * m_u(indexX + 1, indexY + 1, indexZ) +
                                                                                        (rho * g * m_dsY(indexX, indexY));
                }
                else
                {
                    m_uvMatrix.push_back(Eigen::Triplet<double>(
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        m_nXYZ + getIndexXYZ(indexX, indexY, indexZ),
                        1));
                    m_uvVector.coeffRef(m_nXYZ + getIndexXYZ(indexX, indexY, indexZ)) = 0;
                }
            }
        }
    }
}

void Model::solve_uvSystem()
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.setTolerance(1e-4);
    // solver.setMaxIterations(1000);
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(2 * m_nXYZ, 2 * m_nXYZ);
    matrix.setFromTriplets(m_uvMatrix.begin(), m_uvMatrix.end());

    solver.compute(matrix);
    Eigen::Vector<double, Eigen::Dynamic> solution = solver.solve(m_uvVector);

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_uIter(indexX, indexY, indexZ) = solution(getIndexXYZ(indexX, indexY, indexZ));
            }
        }
    }

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                m_vIter(indexX, indexY, indexZ) = solution(m_nXYZ + getIndexXYZ(indexX, indexY, indexZ));
            }
        }
    }
}
