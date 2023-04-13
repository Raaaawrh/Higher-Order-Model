#include "Model.h"
#include <iostream>
#include <fstream>

// Инициализация методов для обработки s, b, h

// Граничные узлы пропускаются, кроме db/dx, db/dy, dh/dx, dh/dy
// так как они нужны для определения ax, ay, az для определени вязкости
// m_eta

// // // Вычисление производных s, b, h

/// @brief Вычисление ds/dx, ds/dy
void Model::calc_ds()
{
    calc_dX(m_s, m_dsX);
    calc_dY(m_s, m_dsY);
}

/// @brief Вычисление db/dx, db/dy
void Model::calc_db()
{
    calc_dX(m_b, m_dbX);
    calc_dY(m_b, m_dbY);

    calc_ddX(m_b, m_ddbX);
    calc_ddY(m_b, m_ddbY);

    calc_dY(m_dbX, m_ddbXY);
}

/// @brief Вычисление h
void Model::calc_h()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_h(indexX, indexY) = m_s(indexX, indexY) - m_b(indexX, indexY);
        }
    }
}

/// @brief Вычисление dh/dx, dh/dy
void Model::calc_dh()
{

    calc_dX(m_h, m_dhX);
    calc_dY(m_h, m_dhY);

    calc_ddX(m_h, m_ddhX);
    calc_ddY(m_h, m_ddhY);

    calc_dY(m_dhX, m_ddhXY);
}

/// @brief Вычисление разом всех производных s, b, h
void Model::calc_derivatives()
{
    calc_ds();
    calc_db();
    calc_dh();
}

// // // Система для решения s

/// @brief Заполнение системы для s
void Model::fill_sSystem()
{

    // m_sMatrix.setZero();
    m_sMatrix.clear();
    m_sVector.setZero();

    // x = minX && x = maxX
    for (size_t indexY{0}; indexY < m_sizeY; indexY++)
    {
        m_sMatrix.push_back(Eigen::Triplet<double>(
            getIndexXY(0, indexY),
            getIndexXY(0, indexY),
            1));
        m_sVector.coeffRef(getIndexXY(0, indexY)) = 0;

        m_sMatrix.push_back(Eigen::Triplet<double>(
            getIndexXY(m_sizeX - 1, indexY),
            getIndexXY(m_sizeX - 1, indexY),
            1));

        m_sVector.coeffRef(getIndexXY(m_sizeX - 1, indexY)) = 0;
    }

    // y = minY && y = maxY

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        m_sMatrix.push_back(Eigen::Triplet<double>(
            getIndexXY(indexX, 0),
            getIndexXY(indexX, 0),
            1));

        m_sVector.coeffRef(getIndexXY(indexX, 0)) = 0;

        m_sMatrix.push_back(Eigen::Triplet<double>(
            getIndexXY(indexX, m_sizeY - 1),
            getIndexXY(indexX, m_sizeY - 1),
            1));

        m_sVector.coeffRef(getIndexXY(indexX, m_sizeY - 1)) = 0;
    }

    for (size_t indexX{1}; indexX < m_sizeX - 1; indexX++)
    {
        for (size_t indexY{1}; indexY < m_sizeY - 1; indexY++)
        {
            m_sMatrix.push_back(Eigen::Triplet<double>(
                getIndexXY(indexX, indexY),
                getIndexXY(indexX - 1, indexY),
                alphaXX(indexX, indexY, -1, 0)));
            m_sMatrix.push_back(Eigen::Triplet<double>(
                getIndexXY(indexX, indexY),
                getIndexXY(indexX, indexY - 1),
                alphaYY(indexX, indexY, 0, -1)));
            m_sMatrix.push_back(Eigen::Triplet<double>(
                getIndexXY(indexX, indexY),
                getIndexXY(indexX, indexY),
                alphaXX(indexX, indexY, 0, 0) + alphaYY(indexX, indexY, 0, 0) + 1 / m_dt));
            m_sMatrix.push_back(Eigen::Triplet<double>(
                getIndexXY(indexX, indexY),
                getIndexXY(indexX, indexY + 1),
                alphaYY(indexX, indexY, 0, 1)));
            m_sMatrix.push_back(Eigen::Triplet<double>(
                getIndexXY(indexX, indexY),
                getIndexXY(indexX + 1, indexY),
                alphaXX(indexX, indexY, 1, 0)));

            m_sVector(getIndexXY(indexX, indexY)) = m_sPrev(indexX, indexY) / m_dt +
                                                    (m_b(indexX, indexY) - m_bPrev(indexX, indexY)) / m_dt +
                                                    (Physics::Ms(m_meshXY.getX(indexX, indexY), m_meshXY.getY(indexX, indexY))) / 2.0;
        }
    }
}

/// @brief Решение системы для s. Решение -- m_sIter.
void Model::solve_S()
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> matrix(m_nXY, m_nXY);
    matrix.setFromTriplets(m_sMatrix.begin(), m_sMatrix.end());

    Eigen::Vector<double, Eigen::Dynamic> guess(m_nXY);

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            guess(getIndexXY(indexX, indexY)) = m_s(indexX, indexY);
        }
    }

    solver.compute(matrix);
    Eigen::Vector<double, Eigen::Dynamic> solution = solver.solveWithGuess(m_sVector, guess);

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            size_t index = indexX * m_sizeY + indexY;
            m_sIter(indexX, indexY) = solution(index);
        }
    }

    // Отсеиваем s < b
    check_sIter();
}

/// @brief Проверка s >= b
void Model::check_sIter()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_sIter(indexX, indexY) = std::max(m_b(indexX, indexY), m_sIter(indexX, indexY));
        }
    }
}

/// @brief Вычисление отклонения на итерациях s
void Model::calc_sErr()
{
    m_sErr = std::numeric_limits<double>::lowest();

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_sErr = std::max(std::fabs(m_s(indexX, indexY) - m_sIter(indexX, indexY)), m_sErr);
        }
    }
}

/// @brief Обновление s на итерациях
void Model::update_sIter()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_s(indexX, indexY) = m_sIter(indexX, indexY);
        }
    }
}

/// @brief Обновление s на предыдущем временном слое
void Model::update_sPrev()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_sPrev(indexX, indexY) = m_s(indexX, indexY);
        }
    }
}

// // // Диффузивность

/// @brief Получение диффузивности Dx в центре ячейки с индексом (__indexX + _dirX/2, _indexY + _dirY/2)
/// @param _indexX
/// @param _indexY
/// @param _dirX
/// @param _dirY
/// @return
double Model::getDXoh(size_t _indexX, size_t _indexY, int _dirX, int _dirY)
{
    double result{0.0};
    result = m_DX(_indexX, _indexY) +
             m_DX(_indexX, _indexY + _dirY) +
             m_DX(_indexX + _dirX, _indexY) +
             m_DX(_indexX + _dirX, _indexY + _dirY);
    return result / 4.0;
}

/// @brief Получение диффузивности Dy в центре ячейки с индексом (_indexX + _dirX / 2, _indexY + _dirY / 2)
/// @param _indexX
/// @param _indexY
/// @param _dirX
/// @param _dirY
/// @return
double Model::getDYoh(size_t _indexX, size_t _indexY, int _dirX, int _dirY)
{
    double result{0.0};
    result = m_DY(_indexX, _indexY) +
             m_DY(_indexX, _indexY + _dirY) +
             m_DY(_indexX + _dirX, _indexY) +
             m_DY(_indexX + _dirX, _indexY + _dirY);
    return result / 4.0;
}

/// @brief Вычисление диффузивностей Dx, Dy в каждой точке.
void Model::calc_D()
{
    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            m_DX(indexX, indexY) = (m_dsX(indexX, indexY) == 0) ? 0 : m_uAverage(indexX, indexY) * m_h(indexX, indexY) / m_dsX(indexX, indexY);
            m_DY(indexX, indexY) = (m_dsY(indexX, indexY) == 0) ? 0 : m_vAverage(indexX, indexY) * m_h(indexX, indexY) / m_dsY(indexX, indexY);
        }
    }
}
