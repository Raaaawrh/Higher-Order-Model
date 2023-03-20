#include <sciplot/sciplot.hpp>
#include <vector>

int main() {
    // Create data points
    std::vector<double> x = {1, 2, 3, 4};
    std::vector<double> y = {1, 2, 3};
    std::vector<std::vector<double>> z = {{1, 2, 3, 4}, {2, 4, 6, 8}, {3, 6, 9, 12}};

    // Create plot object
    sciplot::Plot3D plot;

    // Set the x, y, and z data
    plot.drawPoints(x, y, z);

    sciplot::Figure fig = {{plot}};

    // Show the plot
    sciplot::Canvas canvas = {{fig}};

    canvas.show();

    return 0;
}