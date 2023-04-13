#if !defined(MODEL_H)
#define MODEL_H

#include "Field.h"
#include "Mesh.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include "Physics.h"
#include <string>
#include <vector>

class Model
{
public:
    Model(size_t _sizeX, size_t _szieY, size_t _sizeZ,
          double _minX, double _maxX, double _minY, double _maxY, double _minZ, double _maxZ,
          double _mint, double _dt, double _maxt);
    ~Model();

    void step();
    void startModeling();

private:
    // Параметры сетки и сама сетка
    size_t m_sizeX;
    size_t m_sizeY;
    size_t m_sizeZ;

    double m_minX, m_minY, m_minZ;
    double m_maxX, m_maxY, m_maxZ;

    size_t m_nXY;  // Число узлов в горизонтали XY
    size_t m_nXYZ; // Число узлов в сетке

    // Геттеры для получения индексов в перенумерованном формате
    size_t getIndexXY(size_t _indexX, size_t _indexY) const { return _indexX * m_sizeY + _indexY; }
    size_t getIndexXYZ(size_t _indexX, size_t _indexY, size_t _indexZ) const { return _indexX * m_sizeY * m_sizeZ + _indexY * m_sizeZ + _indexZ; }

    Mesh<2> m_meshXY;
    Mesh<3> m_meshXYZ;

    // // Коэффициенты замен переменных
    Field<3> m_aX;
    Field<3> m_aY;
    Field<2> m_aZ;

    Field<3> m_bX;
    Field<3> m_bY;

    Field<3> m_cXY;
    Field<2> m_cXZ;
    Field<2> m_cYZ;

    // // // Определения в файле ModelMesh.c++

    void calc_aX();
    void calc_aY();
    void calc_aZ();
    void calc_bX();
    void calc_bY();
    void calc_cXY();
    void calc_cXZ();
    void calc_cYZ();

    // Временная компонента
    double m_mint;  // Начальное значение времени
    double m_dt;    // Шаг по времени
    double m_maxt;  // Конечное значение времени
    double m_t;     // Текущее значение времени
    double m_tPrev; // Предыдущее значение времени

    // Поверхность
    Field<2> m_s;     // Значения Z координаты поверхности
    Field<2> m_sIter; // Значения на итераиции решателя
    Field<2> m_sPrev; // Значения на предыдущем временном слое

    Field<2> m_dsX; // Значения частных производных
    Field<2> m_dsY; // Значения частных производных
    void calc_sErr();
    void update_sIter();
    void update_sPrev();
    void calc_ds();

    void fill_sSystem();
    void solve_S();
    void check_sIter();

    // Координаты Z подложки
    Field<2> m_b;
    Field<2> m_bPrev;

    Field<2> m_dbX; // Первые производные
    Field<2> m_dbY;

    Field<2> m_ddbX; // Вторые производные
    Field<2> m_ddbY;
    Field<2> m_ddbXY;

    void calc_db();
    
    // Значения мощности
    Field<2> m_h;
    Field<2> m_hPrev;

    void calc_h(); // Вычисление H=s-b
    void calc_dh();
    
    Field<2> m_dhX; // Первые производные
    Field<2> m_dhY;

    Field<2> m_ddhX; // Вторые производные
    Field<2> m_ddhY;
    Field<2> m_ddhXY;

    // // Диффузивности
    Field<2> m_DX;
    Field<2> m_DY;
    void calc_D();

    // // // Получение их значений в смещенной сетке
    double getDXoh(size_t _indexX, size_t _indexY, int _dirX, int _dirY);
    double getDYoh(size_t _indexX, size_t _indexY, int _dirX, int _dirY);

    // Скорости
    // // В направлении X
    Field<3> m_u;
    Field<3> m_uIter;
    Field<3> m_uPrev;

    Field<3> m_dux;
    Field<3> m_duy;
    Field<3> m_duz;

    void calc_uDerivatives();
    void calc_uErr();

    Field<2> m_uAverage;
    void calc_uAverage();
    void update_uIter();

    void fill_uSystem();
    void solve_u();

    // // В направлении Y
    Field<3> m_v;
    Field<3> m_vIter;
    Field<3> m_vPrev;

    Field<3> m_dvx;
    Field<3> m_dvy;
    Field<3> m_dvz;

    void calc_vDerivatives();
    void calc_vErr();

    Field<2> m_vAverage;
    void calc_vAverage();
    void update_vIter();

    void fill_vSystem();
    void solve_v();

    // // В направлении Z
    Field<3> m_w;
    Field<3> m_wIter;
    Field<3> m_wPrev;

    // // Вязкость
    Field<3> m_eta;
    double const m_eta0 = 1e-30;
    void calc_eta();
    double getEtaoh(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);

    // Выражение Аррениуса
    double A();
    double const n = 3.0;

    //
    double const rho = 910;
    double const g = 9.81;

    // Матрицы и векторы правой части для решателей

    // // Уравнение поверхности
    // Eigen::SparseMatrix<double> m_sMatrix;
    std::vector<Eigen::Triplet<double>> m_sMatrix;
    Eigen::Vector<double, Eigen::Dynamic> m_sVector;
    double m_sTol; // Максимальное отклонение на итерациях
    double m_sErr;

    // // Уравнения компонент скорости
    // Eigen::SparseMatrix<double> m_uvMatrix;
    std::vector<Eigen::Triplet<double>> m_uvMatrix;
    Eigen::Vector<double, Eigen::Dynamic> m_uvVector;
    double m_uvTol; // Максимальное отклонение на итерациях
    double m_uErr;
    double m_vErr;
    double m_uvErr;
    double m_alpha_u;
    double m_alpha_v;
    void fill_uvSystem();
    void solve_uvSystem();

    // Опеределения в файле ModelCoeffs.c++
    // Коэффициенты схемы
    // // 2D схемы
    double alphaXX(size_t _indexX, size_t _indexY, int _directionX, int _directionY);
    double alphaXY(size_t _indexX, size_t _indexY, int _directionX, int _directionY);
    double alphaYX(size_t _indexX, size_t _indexY, int _directionX, int _directionY);
    double alphaYY(size_t _indexX, size_t _indexY, int _directionX, int _directionY);

    // // 3D схемы
    // // // Вычисление шага на неравномерной сетке
    double deltaZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _direction);

    double betaXX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaXY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaXZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaYX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaYY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaYZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaZX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaZY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaZZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double betaZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ);
    double gammaZX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirZ);
    double gammaZY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirY, int _dirZ);

    // Общие методы

    void calc_dX(Field<2> const &_f, Field<2> &_dfX);
    void calc_dY(Field<2> const &_f, Field<2> &_dfY);

    void calc_dX(Field<3> const &_f, Field<3> &_dfX);
    void calc_dY(Field<3> const &_f, Field<3> &_dfY);
    void calc_dZ(Field<3> const &_f, Field<3> &_dfZ);

    void calc_ddX(Field<2> const &_f, Field<2> &_ddfX);
    void calc_ddY(Field<2> const &_f, Field<2> &_ddfY);

    void calc_ddX(Field<3> const &_f, Field<3> &_ddfX);
    void calc_ddY(Field<3> const &_f, Field<3> &_ddfY);
    void calc_ddZ(Field<3> const &_f, Field<3> &_ddfZ);

    void calc_derivatives(); // Вычисление всех производных s, b, h
    void calc_coeff();       // Вычисление всех коэффициентов

    std::string const m_resultsFilepath = "../Results/";
    void saveField(Field<2> const &_field, std::string _filename);
    void saveField(Field<3> const &_field, std::string _filename);
};

#endif // MODEL_H
