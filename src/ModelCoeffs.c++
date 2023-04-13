#include "Model.h"

// Коэффициенты схемы для целевых уравнений

// // 2D

double Model::alphaXX(size_t _indexX, size_t _indexY, int _directionX, int _directionY)
{
    double result = 0;
    switch (_directionX)
    {
    case -1:
        result = (getDXoh(_indexX, _indexY, -1, -1) +
                  getDXoh(_indexX, _indexY, -1, 1));
        break;
    case 0:
        result = (-getDXoh(_indexX, _indexY, -1, -1) +
                  -getDXoh(_indexX, _indexY, -1, 1) +
                  -getDXoh(_indexX, _indexY, 1, -1) +
                  -getDXoh(_indexX, _indexY, 1, 1));
        break;
    case 1:
        result = (getDXoh(_indexX, _indexY, 1, -1) +
                  getDXoh(_indexX, _indexY, 1, 1));
        break;
    }
    return result /
           (2.0 *
            m_meshXY.getDX(_indexX, _indexY) *
            m_meshXY.getDX(_indexX, _indexY));
}

double Model::alphaXY(size_t _indexX, size_t _indexY, int _directionX, int _directionY)
{
    double result = 0;
    switch (_directionX)
    {
    case -1:
        switch (_directionY)
        {
        case -1:
            result = 0;
            break;
        case 0:
            result = 0;
            break;
        case 1:
            result = 0;
            break;
        }

        break;
        ////
    case 0:
        switch (_directionY)
        {
        case -1:
            result = 0;
            break;
        case 0:
            result = 0;
            break;
        case 1:
            result = 0;
            break;
        }

        break;
        ////
    case 1:
        switch (_directionY)
        {
        case -1:
            result = 0;
            break;
        case 0:
            result = 0;
            break;
        case 1:
            result = 0;
            break;
        }
        break;
        ////
    }
    return result;
}

double Model::alphaYX(size_t _indexX, size_t _indexY, int _directionX, int _directionY)
{
    double result = 0;
    switch (_directionX)
    {

    case -1:
        switch (_directionY)
        {
        case -1:
            result = 0;
            break;
        case 0:
            result = 0;
            break;
        case 1:
            result = 0;
            break;
        }

        break;
        ////
    case 0:
        switch (_directionY)
        {
        case -1:
            result = 0;
            break;
        case 0:
            result = 0;
            break;
        case 1:
            result = 0;
            break;
        }

        break;
        ////
    case 1:
        switch (_directionY)
        {
        case -1:
            result = 0;
            break;
        case 0:
            result = 0;
            break;
        case 1:
            result = 0;
            break;
        }
        break;
        ////
    }
    return result;
}

double Model::alphaYY(size_t _indexX, size_t _indexY, int _directionX, int _directionY)
{
    double result = 0;
    switch (_directionY)
    {
    case -1:
        result = (getDYoh(_indexX, _indexY, -1, -1) +
                  getDYoh(_indexX, _indexY, 1, -1));
        break;
    case 0:
        result = (-getDYoh(_indexX, _indexY, -1, -1) +
                  -getDYoh(_indexX, _indexY, 1, -1) +
                  -getDYoh(_indexX, _indexY, -1, 1) +
                  -getDYoh(_indexX, _indexY, 1, 1));
        break;
    case 1:
        result = (getDYoh(_indexX, _indexY, -1, 1) +
                  getDYoh(_indexX, _indexY, 1, 1));
        break;
    }
    return result /
           (2.0 *
            m_meshXY.getDY(_indexX, _indexY) *
            m_meshXY.getDY(_indexX, _indexY));
}

// 3D коэффициенты схемы

double Model::deltaZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _direction)
{
    double result{0.0};

    switch (_direction)
    {
    case -1:
        result = m_meshXYZ.getZ(_indexX, _indexY, _indexZ) - m_meshXYZ.getZ(_indexX, _indexY, _indexZ - 1);
        break;
    case 0:
        result = (m_meshXYZ.getZ(_indexX, _indexY, _indexZ + 1) - m_meshXYZ.getZ(_indexX, _indexY, _indexZ - 1)) /
                 2.0;
        break;
    case 1:
        result = m_meshXYZ.getZ(_indexX, _indexY, _indexZ + 1) - m_meshXYZ.getZ(_indexX, _indexY, _indexZ);
        break;
    }
    return result;
}

double Model::betaXX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};

    switch (_dirX)
    {
    case -1:
        result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1);

        break;
    case 0:
        result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1)) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1)) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1)) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
        break;
    case 1:
        result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1);
        break;
    }
    return result /
           (2.0 *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            2.0 *
            deltaZ(_indexX, _indexY, _indexZ, 0));
}

double Model::betaXY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};
    switch (_dirX)
    {
    case -1:
        switch (_dirY)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1);
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1);
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            break;
        }
        break;
    case 0:
        switch (_dirY)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1);
            break;
        }
        break;
    case 1:
        switch (_dirY)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1);
            break;
        }
        break;
    }
    return result /
           (4.0 *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            2.0 *
            deltaZ(_indexX, _indexY, _indexZ, 0));
}
double Model::betaXZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};
    switch (_dirX)
    {
    case -1:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1));
            break;
        case 0:
            result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            ;
            break;
        case 1:
            result = -deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            break;
        }
        break;
    case 0:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    case 1:
        switch (_dirZ)
        {
        case -1:
            result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     -deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));

            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    }
    return result /
           (4.0 *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            deltaZ(_indexX, _indexY, _indexZ, 0) *
            2.0);
}
double Model::betaYX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};
    switch (_dirX)
    {
    case -1:
        switch (_dirY)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1);
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            break;
        }
        break;
    case 0:
        switch (_dirY)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1);
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    case 1:
        switch (_dirY)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1);
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1);
            break;
        }
        break;
    }
    return result /
           (4.0 *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            2.0 *
            deltaZ(_indexX, _indexY, _indexZ, 0));
}

double Model::betaYY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};

    switch (_dirY)
    {
    case -1:
        result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1);

        break;
    case 0:
        result = deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1)) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1)) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1)) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * (-getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
        break;
    case 1:
        result = deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                 deltaZ(_indexX, _indexY, _indexZ, 1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                 deltaZ(_indexX, _indexY, _indexZ, -1) * getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1);
        break;
    }
    return result /
           (2.0 *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            2.0 *
            deltaZ(_indexX, _indexY, _indexZ, 0));
}
double Model::betaYZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};
    switch (_dirY)
    {
    case -1:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1));
            break;
        case 0:
            result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        case 1:
            result = -deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        }
        break;
    case 0:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    case 1:
        switch (_dirZ)
        {
        case -1:
            result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     -deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));

            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    }
    return result /
           (4.0 *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            deltaZ(_indexX, _indexY, _indexZ, 0) *
            2.0);
}
double Model::betaZX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};
    switch (_dirX)
    {
    case -1:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1));
            break;
        case 0:
            result = 2.0 *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1)) +
                     -deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            break;
        case 1:
            result = -deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
            break;
        }
        break;
    case 0:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = -2.0 *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    case 1:
        switch (_dirZ)
        {
        case -1:
            result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = -2.0 *
                         (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     -deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));

            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    }
    return result /
           (4.0 *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            deltaZ(_indexX, _indexY, _indexZ, 0) *
            2.0);
}
double Model::betaZY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};
    switch (_dirY)
    {
    case -1:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1));
            break;
        case 0:
            result = 2.0 *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1)) +
                     -deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        case 1:
            result = -deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
            break;
        }
        break;
    case 0:
        switch (_dirZ)
        {
        case -1:
            result = deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = -2.0 *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (-getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    case 1:
        switch (_dirZ)
        {
        case -1:
            result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1));
            break;
        case 0:
            result = -2.0 *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1)) +
                     deltaZ(_indexX, _indexY, _indexZ, 1) /
                         deltaZ(_indexX, _indexY, _indexZ, -1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                     -deltaZ(_indexX, _indexY, _indexZ, -1) /
                         deltaZ(_indexX, _indexY, _indexZ, 1) *
                         (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                          getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));

            break;
        case 1:
            result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
            break;
        }
        break;
    }
    return result /
           (4.0 *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            deltaZ(_indexX, _indexY, _indexZ, 0) *
            2.0);
}
double Model::betaZZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};

    switch (_dirZ)
    {
    case -1:
        result = 2.0 *
                     deltaZ(_indexX, _indexY, _indexZ, 1) /
                     deltaZ(_indexX, _indexY, _indexZ, -1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                 (1 -
                  deltaZ(_indexX, _indexY, _indexZ, -1) /
                      deltaZ(_indexX, _indexY, _indexZ, 1)) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
        break;
    case 0:
        result = -(1 +
                   2 * deltaZ(_indexX, _indexY, _indexZ, 1) /
                       deltaZ(_indexX, _indexY, _indexZ, -1) +
                   -1 * deltaZ(_indexX, _indexY, _indexZ, -1) /
                       deltaZ(_indexX, _indexY, _indexZ, 1)) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                 -(1 +
                   2 * deltaZ(_indexX, _indexY, _indexZ, -1) /
                       deltaZ(_indexX, _indexY, _indexZ, 1) +
                   -1 * deltaZ(_indexX, _indexY, _indexZ, 1) /
                       deltaZ(_indexX, _indexY, _indexZ, -1)) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
        break;
    case 1:
        result = (1 -
                  deltaZ(_indexX, _indexY, _indexZ, 1) /
                      deltaZ(_indexX, _indexY, _indexZ, -1)) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) +
                 2.0 *
                     deltaZ(_indexX, _indexY, _indexZ, -1) /
                     deltaZ(_indexX, _indexY, _indexZ, 1) *
                     (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                      getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
        break;
    }
    result /= (4.0 * 2.0 * deltaZ(_indexX, _indexY, _indexZ, 0));
    return result;
}

double Model::betaZ(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirY, int _dirZ)
{
    double result{0.0};

    double g_ijk = deltaZ(_indexX, _indexY, _indexZ, -1) /
                       (2 * deltaZ(_indexX, _indexY, _indexZ, 0)) *
                       (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                        getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                        getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                        getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1)) /
                       4 +
                   2 * deltaZ(_indexX, _indexY, _indexZ, 1) /
                       (2 * deltaZ(_indexX, _indexY, _indexZ, 0)) *
                       (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                        getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                        getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                        getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1)) /
                       4;

    switch (_dirZ)
    {
    case -1:
        result = -deltaZ(_indexX, _indexY, _indexZ, 1) /
                 deltaZ(_indexX, _indexY, _indexZ, -1) *
                 g_ijk;
        break;
    case 0:
        result = 2 * deltaZ(_indexX, _indexY, _indexZ, 0) *
                 (deltaZ(_indexX, _indexY, _indexZ, 1) - deltaZ(_indexX, _indexY, _indexZ, -1)) /
                 (deltaZ(_indexX, _indexY, _indexZ, 1) * deltaZ(_indexX, _indexY, _indexZ, -1)) *
                 g_ijk;
        break;
    case 1:
        result = deltaZ(_indexX, _indexY, _indexZ, -1) /
                 deltaZ(_indexX, _indexY, _indexZ, 1) *
                 g_ijk;
        break;
    }

    return result / (2 * deltaZ(_indexX, _indexY, _indexZ, 0));
}

double Model::gammaZX(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirX, int _dirZ)
{
    double result{0.0};

    switch (_dirX)
    {
    case -1:
        result = 2 *
                 (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1));
        break;
    case 0:
        result = -2 *
                 (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
    case 1:
        result = -2 *
                 (getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
        break;
    }
    return result /
           (2 *
            4 *
            m_meshXYZ.getDX(_indexX, _indexY, _indexZ) *
            deltaZ(_indexX, _indexY, _indexZ, 0));
}
double Model::gammaZY(size_t _indexX, size_t _indexY, size_t _indexZ, int _dirY, int _dirZ)
{
    double result{0.0};

    switch (_dirY)
    {
    case -1:
        result = 2 *
                 (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1));
        break;
    case 0:
        result = -2 *
                 (getEtaoh(_indexX, _indexY, _indexZ, -1, -1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, -1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, -1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, -1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
    case 1:
        result = -2 *
                 (getEtaoh(_indexX, _indexY, _indexZ, -1, 1, -1) +
                  getEtaoh(_indexX, _indexY, _indexZ, 1, 1, -1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, -1, 1, 1) +
                  -getEtaoh(_indexX, _indexY, _indexZ, 1, 1, 1));
        break;
    }
    return result /
           (2 *
            4 *
            m_meshXYZ.getDY(_indexX, _indexY, _indexZ) *
            deltaZ(_indexX, _indexY, _indexZ, 0));
}
