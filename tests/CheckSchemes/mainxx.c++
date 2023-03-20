#include <sciplot/sciplot.hpp>
#include <vector>
#include <cmath>
using namespace sciplot;

double func(double x)
{
    double res{0};
    res = sin(x);
    return res;
}

double g(double x)
{
    double res{0};
    res = sin(x);
    return res;
}

double g12(std::vector<double> &x, std::size_t i)
{
    double res{0};
    res = g(x[i]) + g(x[i + 1]);
    return res / 2;
}

int main(int argc, char **argv)
{
    std::size_t size{200};

    double xmin{-10}, xmax{10};
    double xLength{xmax - xmin};
    double dx{xLength / (size - 1)};

    std::vector<double> x(size), f(size), fExact(size);

    for (std::size_t i{0}; i < size; i++)
    {
        x[i] = xmin + dx * i;
        fExact[i] = cos(x[i]) * cos(x[i]) - sin(x[i]) * sin(x[i]);
    }

    for (std::size_t i{1}; i < size - 1; i++)
    {
        f[i] = func(x[i - 1]) * g12(x, i - 1) +
               func(x[i]) * (-g12(x, i - 1) - g12(x, i)) +
               func(x[i + 1]) * g12(x, i);
        f[i] /= (dx * dx);
    }

    Plot2D plot;

    plot.xrange(xmin, xmax);
    plot.yrange(-3, 3);

    // plot.drawCurve(xs, fs);
    plot.drawCurve(x, fExact);
    plot.drawCurve(x, f);
    // plot.drawCurve(xs, fs);
    plot.grid()
        .show();

    // Create figure to hold plot
    Figure fig = {{plot}};
    // Create canvas to hold figure
    Canvas canvas = {{fig}};

    // Show the plot in a pop-up window
    canvas.show();

    // Save the plot to a PDF file
    canvas.save("example-sine-functions.pdf");
}