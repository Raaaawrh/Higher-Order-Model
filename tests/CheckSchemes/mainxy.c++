#include <sciplot/sciplot.hpp>
#include <vector>
#include <cmath>
using namespace sciplot;

int main(int argc, char **argv)
{
    std::size_t xSize{200}, ySize{200};

    double xmin{-10}, xmax{10},
        ymin{-10}, ymax{10};

    double xLength{xmax - xmin}, yLength{ymax - ymin};

    double dx{xLength / (xSize - 1)}, dy{yLength / (ySize - 1)};

    std::vector<std::vector<double>> x(xSize, std::vector<double>(ySize, 0.0)),
        y(xSize, std::vector<double>(ySize, 0.0)),
        z(xSize, std::vector<double>(ySize, 0.0)),
        zScheme(xSize, std::vector<double>(ySize, 0.0));

    for (std::size_t i{0}; i < xSize; i++)
    {
        for (std::size_t j{0}; j < ySize; j++)
        {
            x[i][j] = xmin + i * dx;
            y[i][j] = ymin + j * dy;
        }
    }

    Plot3D plot;

    plot.xrange(xmin, xmax);
    plot.yrange(ymin, ymax);
    plot.zrange(-3, 3);

    // plot.drawCurve(xs, fs);
    plot.drawPoints(x, y, z).lineStyle(0);
    // plot.drawCurve(x, f);

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