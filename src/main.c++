#include "Field.h"
#include "Mesh.h"
#include "Model.h"

int main(int argc, char const *argv[])
{
    size_t sizeX{61};
    size_t sizeY{61};
    size_t sizeZ{41};

    double minX{0}, maxX{1500};
    double minY{0}, maxY{1500};
    double minZ{0}, maxZ{1};

    double mint{0}, dt{20}, maxt{1000};

    Model model(sizeX, sizeY, sizeZ, minX, maxX, minY, maxY, minZ, maxZ, mint, dt, maxt);
    
    //model.test();
    model.startModeling();

    return 0;
}
