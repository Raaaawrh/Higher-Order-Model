#include "Model.h"
#include <eigen3/Eigen/Sparse>
#include <limits>
#include <cmath>
#include <fstream>

/////////////// PUBLIC ///////////////

// Конструкторы
Model::Model(size_t _sizeX, size_t _sizeY, size_t _sizeZ,
             double _minX, double _maxX, double _minY, double _maxY, double _minZ, double _maxZ,
             double _mint, double _dt, double _maxt)
    : m_sizeX{_sizeX}, m_sizeY{_sizeY}, m_sizeZ{_sizeZ},
      m_minX{_minX}, m_minY{_minY}, m_minZ{_minZ},
      m_maxX{_maxX}, m_maxY{_maxY}, m_maxZ{_maxZ},
      m_nXY{m_sizeX * m_sizeY}, m_nXYZ{m_sizeX * m_sizeY * m_sizeZ},

      m_meshXY{Mesh<2>(m_minX, m_minY, m_maxX, m_maxY, m_sizeX, m_sizeY)},
      m_meshXYZ{Mesh<3>(m_minX, m_minY, m_minZ, m_maxX, m_maxY, m_maxZ, m_sizeX, m_sizeY, m_sizeZ)},

      m_aX{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},
      m_aY{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},
      m_aZ{Field<2>(m_sizeX, m_sizeY)},

      m_bX{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},
      m_bY{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},

      m_cXY{Field<3>(Field<3>(m_sizeX, m_sizeY, m_sizeZ))},
      m_cXZ{Field<2>(m_sizeX, m_sizeY)},
      m_cYZ{Field<2>(m_sizeX, m_sizeY)},

      m_mint{_mint}, m_dt{_dt}, m_maxt{_maxt},
      m_t{m_mint}, m_tPrev{m_t},

      m_s{Field<2>(m_sizeX, m_sizeY)}, m_sIter{Field<2>(m_sizeX, m_sizeY)}, m_sPrev{Field<2>(m_sizeX, m_sizeY)},
      m_dsX{Field<2>(m_sizeX, m_sizeY)}, m_dsY{Field<2>(m_sizeX, m_sizeY)},

      m_b{Field<2>(m_sizeX, m_sizeY)}, m_bPrev{Field<2>(m_sizeX, m_sizeY)},
      m_dbX{Field<2>(m_sizeX, m_sizeY)}, m_dbY{Field<2>(m_sizeX, m_sizeY)},
      m_ddbX{Field<2>(m_sizeX, m_sizeY)}, m_ddbY{Field<2>(m_sizeX, m_sizeY)}, m_ddbXY{Field<2>(m_sizeX, m_sizeY)},

      m_h{Field<2>(m_sizeX, m_sizeY)}, m_hPrev{Field<2>(m_sizeX, m_sizeY)},
      m_dhX{Field<2>(m_sizeX, m_sizeY)}, m_dhY{Field<2>(m_sizeX, m_sizeY)},
      m_ddhX{Field<2>(m_sizeX, m_sizeY)}, m_ddhY{Field<2>(m_sizeX, m_sizeY)}, m_ddhXY{Field<2>(m_sizeX, m_sizeY)},

      m_DX{Field<2>(m_sizeX, m_sizeY)}, m_DY{Field<2>(m_sizeX, m_sizeY)},

      m_u{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_uIter{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_uPrev{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},
      m_uAverage{Field<2>(m_sizeX, m_sizeY)},
      m_dux{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_duy{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_duz{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},

      m_v{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_vIter{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_vPrev{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},
      m_vAverage{Field<2>(m_sizeX, m_sizeY)},
      m_dvx{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_dvy{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_dvz{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},

      m_w{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_wIter{Field<3>(m_sizeX, m_sizeY, m_sizeZ)}, m_wPrev{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},

      m_eta{Field<3>(m_sizeX, m_sizeY, m_sizeZ)},

      // m_sMatrix{Eigen::SparseMatrix<double>(m_nXY, m_nXY)},
      m_sMatrix{std::vector<Eigen::Triplet<double>>()},
      m_sVector{Eigen::Vector<double, Eigen::Dynamic>(m_nXY)},
      m_sTol{1e-3}, m_sErr{std::numeric_limits<double>::max()},

      // m_uvMatrix{Eigen::SparseMatrix<double>(m_nXYZ, m_nXYZ)},
      m_uvMatrix{std::vector<Eigen::Triplet<double>>()},
      m_uvVector{Eigen::Vector<double, Eigen::Dynamic>(m_nXYZ)},
      m_uvTol{1e-6}, m_uErr{std::numeric_limits<double>::max()}, m_vErr{std::numeric_limits<double>::max()},
      m_uvErr{1e-4}
{
    Eigen::setNbThreads(6);
}

// Деструкторы

Model::~Model() {}

// Запуск моделирования
void Model::startModeling()
{
    while (m_t <= m_maxt)
    {
        step();
    }
}

/////////////// PRIVATE ///////////////

double Model::A()
{
    return 10e-16;
}

// Шаг моделирования
void Model::step()
{
    // Итерации по получаемому полю поверхности s

    m_sErr = std::numeric_limits<double>::max();
    while (m_sTol < m_sErr)
    {
        // D
        calc_h();
        calc_derivatives();
        calc_uAverage();
        calc_vAverage();
        saveField(m_uAverage, "uAverage.txt");
        saveField(m_vAverage, "vAverage.txt");

        calc_D();

        saveField(m_DX, "DX.txt");
        saveField(m_DY, "DY.txt");

        fill_sSystem();
        solve_S();

        calc_sErr();
        saveField(m_sIter, "surface.txt");
        update_sIter();

        m_uErr = std::numeric_limits<double>::max();
        m_vErr = std::numeric_limits<double>::max();

        while (m_uvTol < m_uErr || m_uvTol < m_vErr)
        {
            calc_h();
            calc_derivatives();
            calc_coeff();

            calc_uDerivatives();
            calc_vDerivatives();

            calc_eta();

            fill_uSystem();
            solve_u();

            fill_vSystem();
            solve_v();

            calc_uErr();
            calc_vErr();

            update_uIter();
            update_vIter();
        }
    }

    // Обновление поверхности на предыдущем временном слое
    update_sPrev();

    // Смена временного слоя
    m_t += m_dt;

    saveField(m_s, "surface.txt");
}

void Model::saveField(Field<2> const &_field, std::string _filename)
{
    std::string filepath = m_resultsFilepath + _filename;
    std::ofstream file(filepath);
    file << m_sizeX << '\t'
         << m_minX << '\t'
         << m_maxX << '\t'
         << m_sizeY << '\t'
         << m_minY << '\t'
         << m_maxY << '\t'
         << m_t << std::endl;

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            file << _field(indexX, indexY) << '\t';
        }
        file << std::endl;
    }

    file.close();
}

void Model::saveField(Field<3> const &_field, std::string _filename)
{
    std::string filepath = m_resultsFilepath + _filename;
    std::ofstream file(filepath);
    file << m_sizeX << '\t'
         << m_minX << '\t'
         << m_maxX << '\t'
         << m_sizeY << '\t'
         << m_minY << '\t'
         << m_maxY << '\t'
         << m_sizeZ << '\t'
         << m_t << std::endl;

    for (size_t indexX{0}; indexX < m_sizeX; indexX++)
    {
        for (size_t indexY{0}; indexY < m_sizeY; indexY++)
        {
            for (size_t indexZ{0}; indexZ < m_sizeZ; indexZ++)
            {
                file << _field(indexX, indexY, indexZ) << '\t';
            }
            file << std::endl;
        }
    }

    file.close();
}

void Model::calc_dX(Field<2> const &_f, Field<2> &_dfX)
{
    // x = minX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        _dfX(0, iY) = (-3 * _f(0, iY) +
                       4 * _f(1, iY) -
                       _f(2, iY)) /
                      (2 *
                       m_meshXY.getDX(0, iY));
    }

    for (size_t iX{1}; iX < m_sizeX - 1; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            _dfX(iX, iY) = (_f(iX + 1, iY) - _f(iX - 1, iY)) /
                           (2 *
                            m_meshXY.getDX(iX, iY));
        }
    }

    // x = maxX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        _dfX(m_sizeX - 1, iY) = (3 * _f(m_sizeX - 1, iY) -
                                 4 * _f(m_sizeX - 2, iY) +
                                 _f(m_sizeX - 3, iY)) /
                                (2 *
                                 m_meshXY.getDX(m_sizeX - 2, iY));
    }
}

void Model::calc_dY(Field<2> const &_f, Field<2> &_dfY)
{
    // y = minY
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        _dfY(iX, 0) = (-3 * _f(iX, 0) +
                       4 * _f(iX, 1) -
                       _f(iX, 2)) /
                      (2 *
                       m_meshXY.getDY(iX, 0));
    }

    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{1}; iY < m_sizeY - 1; iY++)
        {
            _dfY(iX, iY) = (_f(iX, iY + 1) - _f(iX, iY - 1)) /
                           (2 *
                            m_meshXY.getDY(iX, iY));
        }
    }

    // x = maxX
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        _dfY(iX, m_sizeY - 1) = (3 * _f(iX, m_sizeY - 1) -
                                 4 * _f(iX, m_sizeY - 2) +
                                 _f(iX, m_sizeY - 3)) /
                                (2 *
                                 m_meshXY.getDY(iX, m_sizeY - 2));
    }
}

void Model::calc_dX(Field<3> const &_f, Field<3> &_dfX)
{
    // x = minX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfX(0, iY, iZ) = (-3 * _f(0, iY, iZ) +
                               4 * _f(1, iY, iZ) -
                               _f(2, iY, iZ)) /
                              (2 *
                               m_meshXYZ.getDX(0, iY, iZ));
        }
    }

    for (size_t iX{1}; iX < m_sizeX - 1; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
            {
                _dfX(iX, iY, iZ) = (_f(iX + 1, iY, iZ) - _f(iX - 1, iY, iZ)) /
                                   (2 *
                                    m_meshXYZ.getDX(iX, iY, iZ));
            }
        }
    }

    // x = maxX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfX(m_sizeX - 1, iY, iZ) = (3 * _f(m_sizeX - 1, iY, iZ) -
                                         4 * _f(m_sizeX - 2, iY, iZ) +
                                         _f(m_sizeX - 3, iY, iZ)) /
                                        (2 *
                                         m_meshXYZ.getDX(m_sizeX - 2, iY, iZ));
        }
    }
}

void Model::calc_dY(Field<3> const &_f, Field<3> &_dfY)
{
    // y = minY
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfY(iX, 0, iZ) = (-3 * _f(iX, 0, iZ) +
                               4 * _f(iX, 1, iZ) -
                               _f(iX, 2, iZ)) /
                              (2 *
                               m_meshXYZ.getDY(iX, 0, iZ));
        }
    }

    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{1}; iY < m_sizeY - 1; iY++)
        {
            for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
            {
                _dfY(iX, iY, iZ) = (_f(iX, iY + 1, iZ) - _f(iX, iY - 1, iZ)) /
                                   (2 *
                                    m_meshXYZ.getDY(iX, iY, iZ));
            }
        }
    }

    // x = maxX
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfY(iX, m_sizeY - 1, iZ) = (3 * _f(iX, m_sizeY - 1, iZ) -
                                         4 * _f(iX, m_sizeY - 2, iZ) +
                                         _f(iX, m_sizeY - 3, iZ)) /
                                        (2 *
                                         m_meshXYZ.getDY(iX, m_sizeY - 2, iZ));
        }
    }
}

void Model::calc_dZ(Field<3> const &_f, Field<3> &_dfZ)
{
    // z = minZ
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            _dfZ(iX, iY, 0) = _f(iX, iY, 0) *
                                  ((-deltaZ(iX, iY, 1, -1) - 2 * deltaZ(iX, iY, 1, 0)) /
                                   (deltaZ(iX, iY, 1, -1) * 2 * deltaZ(iX, iY, 1, 0))) +
                              _f(iX, iY, 1) *
                                  ((2 * deltaZ(iX, iY, 1, 0)) /
                                   (deltaZ(iX, iY, 1, -1) * deltaZ(iX, iY, 1, 1))) +
                              _f(iX, iY, 2) *
                                  ((-deltaZ(iX, iY, 1, -1)) /
                                   (deltaZ(iX, iY, 1, 1) * 2 * deltaZ(iX, iY, 1, 0)));
        }
    }

    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            for (size_t iZ{1}; iZ < m_sizeZ - 1; iZ++)
            {
                _dfZ(iX, iY, iZ) = _f(iX, iY, iZ - 1) *
                                       (-(deltaZ(iX, iY, iZ, 1)) /
                                        (2 * deltaZ(iX, iY, iZ, 0) * deltaZ(iX, iY, iZ, -1))) +
                                   _f(iX, iY, iZ) *
                                       ((deltaZ(iX, iY, iZ, 1) - deltaZ(iX, iY, iZ, -1)) /
                                        (deltaZ(iX, iY, iZ, -1) * deltaZ(iX, iY, iZ, 1))) +
                                   _f(iX, iY, iZ + 1) *
                                       ((deltaZ(iX, iY, iZ, -1)) /
                                        (2 * deltaZ(iX, iY, iZ, 0) * deltaZ(iX, iY, iZ, 1)));
            }
        }
    }

    // x = maxX
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{0}; iY < m_sizeZ; iY++)
        {
            _dfZ(iX, iY, m_sizeZ - 1) = _f(iX, iY, m_sizeZ - 1) *
                                            ((deltaZ(iX, iY, m_sizeZ - 2, 1) + 2 * deltaZ(iX, iY, m_sizeZ - 2, 0)) /
                                             (deltaZ(iX, iY, m_sizeZ - 2, 1) * 2 * deltaZ(iX, iY, m_sizeZ - 2, 0))) +
                                        _f(iX, iY, m_sizeZ - 2) *
                                            ((-2 * deltaZ(iX, iY, m_sizeZ - 2, 0)) /
                                             (deltaZ(iX, iY, m_sizeZ - 2, 1) * deltaZ(iX, iY, m_sizeZ - 2, -1))) +
                                        _f(iX, iY, m_sizeZ - 3) *
                                            ((deltaZ(iX, iY, m_sizeZ - 2, -1)) /
                                             (deltaZ(iX, iY, m_sizeZ - 2, -1) * 2 * deltaZ(iX, iY, m_sizeZ - 2, 0)));
        }
    }
}

void Model::calc_ddX(Field<2> const &_f, Field<2> &_dfX)
{
    // x = minX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        _dfX(0, iY) = (2 * _f(0, iY) -
                       5 * _f(1, iY) +
                       4 * _f(2, iY) -
                       _f(3, iY)) /
                      (m_meshXY.getDX(0, iY) *
                       m_meshXY.getDX(0, iY) *
                       m_meshXY.getDX(0, iY));
    }

    for (size_t iX{1}; iX < m_sizeX - 1; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            _dfX(iX, iY) = (_f(iX + 1, iY) - 2 * _f(iX, iY) + _f(iX - 1, iY)) /
                           (m_meshXY.getDX(iX, iY) *
                            m_meshXY.getDX(iX, iY));
        }
    }

    // x = maxX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        _dfX(m_sizeX - 1, iY) = (2 * _f(m_sizeX - 1, iY) -
                                 5 * _f(m_sizeX - 2, iY) +
                                 4 * _f(m_sizeX - 3, iY) -
                                 _f(m_sizeZ - 4, iY)) /
                                (m_meshXY.getDX(m_sizeX - 2, iY) *
                                 m_meshXY.getDX(m_sizeX - 2, iY) *
                                 m_meshXY.getDX(m_sizeX - 2, iY));
    }
}

void Model::calc_ddY(Field<2> const &_f, Field<2> &_dfY)
{
    // y = minY
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        _dfY(iX, 0) = (2 * _f(iX, 0) -
                       5 * _f(iX, 1) +
                       4 * _f(iX, 2) -
                       _f(iX, 3)) /
                      (m_meshXY.getDY(iX, 0) *
                       m_meshXY.getDY(iX, 0) *
                       m_meshXY.getDY(iX, 0));
    }

    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{1}; iY < m_sizeY - 1; iY++)
        {
            _dfY(iX, iY) = (_f(iX, iY + 1) - 2 * _f(iX, iY) + _f(iX, iY - 1)) /
                           (m_meshXY.getDY(iX, iY) *
                            m_meshXY.getDY(iX, iY));
        }
    }

    // x = maxX
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        _dfY(iX, m_sizeY - 1) = (2 * _f(iX, m_sizeY - 1) -
                                 5 * _f(iX, m_sizeY - 2) +
                                 4 * _f(iX, m_sizeY - 3) -
                                 _f(iX, m_sizeY - 4)) /
                                (m_meshXY.getDY(iX, m_sizeY - 2) *
                                 m_meshXY.getDY(iX, m_sizeY - 2) *
                                 m_meshXY.getDY(iX, m_sizeY - 2));
    }
}

void Model::calc_ddX(Field<3> const &_f, Field<3> &_dfX)
{
    // x = minX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfX(0, iY, iZ) = (2 * _f(0, iY, iZ) -
                               5 * _f(1, iY, iZ) +
                               4 * _f(2, iY, iZ) -
                               _f(3, iY, iZ)) /
                              (m_meshXYZ.getDX(0, iY, iZ) *
                               m_meshXYZ.getDX(0, iY, iZ) *
                               m_meshXYZ.getDX(0, iY, iZ));
        }
    }

    for (size_t iX{1}; iX < m_sizeX - 1; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
            {
                _dfX(iX, iY, iZ) = (_f(iX + 1, iY, iZ) - 2 * _f(iX, iY, iZ) + _f(iX - 1, iY, iZ)) /
                                   (m_meshXYZ.getDX(iX, iY, iZ) *
                                    m_meshXYZ.getDX(iX, iY, iZ));
            }
        }
    }

    // x = maxX
    for (size_t iY{0}; iY < m_sizeY; iY++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfX(m_sizeX - 1, iY, iZ) = (2 * _f(m_sizeX - 1, iY, iZ) -
                                         5 * _f(m_sizeX - 2, iY, iZ) +
                                         4 * _f(m_sizeX - 3, iY, iZ) -
                                         _f(m_sizeZ - 4, iY, iZ)) /
                                        (m_meshXYZ.getDX(m_sizeX - 2, iY, iZ) *
                                         m_meshXYZ.getDX(m_sizeX - 2, iY, iZ) *
                                         m_meshXYZ.getDX(m_sizeX - 2, iY, iZ));
        }
    }
}

void Model::calc_ddY(Field<3> const &_f, Field<3> &_dfY)
{
    // y = minY
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfY(iX, 0, iZ) = (2 * _f(iX, 0, iZ) -
                               5 * _f(iX, 1, iZ) +
                               4 * _f(iX, 2, iZ) -
                               _f(iX, 3, iZ)) /
                              (m_meshXYZ.getDY(iX, 0, iZ) *
                               m_meshXYZ.getDY(iX, 0, iZ) *
                               m_meshXYZ.getDY(iX, 0, iZ));
        }
    }

    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{1}; iY < m_sizeY - 1; iY++)
        {
            for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
            {
                _dfY(iX, iY, iZ) = (_f(iX, iY + 1, iZ) - 2 * _f(iX, iY, iZ) + _f(iX, iY - 1, iZ)) /
                                   (m_meshXYZ.getDY(iX, iY, iZ) *
                                    m_meshXYZ.getDY(iX, iY, iZ));
            }
        }
    }

    // x = maxX
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iZ{0}; iZ < m_sizeZ; iZ++)
        {
            _dfY(iX, m_sizeY - 1, iZ) = (2 * _f(iX, m_sizeY - 1, iZ) -
                                         5 * _f(iX, m_sizeY - 2, iZ) +
                                         4 * _f(iX, m_sizeY - 3, iZ) -
                                         _f(iX, m_sizeY - 4, iZ)) /
                                        (m_meshXYZ.getDY(iX, m_sizeY - 2, iZ) *
                                         m_meshXYZ.getDY(iX, m_sizeY - 2, iZ) *
                                         m_meshXYZ.getDY(iX, m_sizeY - 2, iZ));
        }
    }
}

void Model::calc_ddZ(Field<3> const &_f, Field<3> &_ddfZ)
{
    for (size_t iX{0}; iX < m_sizeX; iX++)
    {
        for (size_t iY{0}; iY < m_sizeY; iY++)
        {
            for (size_t iZ{1}; iZ < m_sizeZ - 1; iZ++)
            {
                _ddfZ(iX, iY, iZ) = _f(iX, iY, iZ) *
                                        ((1) /
                                         (deltaZ(iX, iY, iZ, -1) * deltaZ(iX, iY, iZ, 0))) +
                                    _f(iX, iY, iZ) *
                                        (-(2) /
                                         (deltaZ(iX, iY, iZ, -1) * deltaZ(iX, iY, iZ, 1))) +
                                    _f(iX, iY, iZ) *
                                        ((1) /
                                         (deltaZ(iX, iY, iZ, 1) * deltaZ(iX, iY, iZ, -1)));
            }
        }
    }
}