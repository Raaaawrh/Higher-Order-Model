#if !defined(PHYSICS_H)
#define PHYSICS_H

class Physics
{
public:
    Physics();
    ~Physics();

    static constexpr double u_b{0.0}, v_b{0.0};
    static constexpr double Mmax{0.5};
    static constexpr double Sb{1e-2};
    static constexpr double Rel{450};
    static constexpr double Tmin{238.15};
    static constexpr double St{1.67e-2};
    static constexpr double x_sum{750}, y_sum{750};

    static double Ms(double _x, double _y);
    static double A(double _temp = Tmin);
    static double T(double _x, double _y);

private:
};

#endif // PHYSICS_H
