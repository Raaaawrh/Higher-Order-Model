#include "Physics.h"

#include <cmath>

Physics::Physics() {}

Physics::~Physics() {}

double Physics::Ms(double _x, double _y)
{
    double result;
    result = std::min(Mmax, Sb * (Rel - std::sqrt((_x - x_sum) * (_x - x_sum) + (_y - y_sum) * (_y - y_sum))));
    return result;
}

double Physics::A(double _temp)
{
    return 1e-16;
}

double Physics::T(double _x, double _y)
{
    double result;
    result = Tmin + St * std::sqrt((_x - x_sum) * (_x - x_sum) + (_y - y_sum) * (_y - y_sum));
    return result;
}