#include "Field.h"
#include "Mesh.h"

int main(int argc, char const *argv[])
{
    Field<1> field1(11);
    Field<2> field2(7, 6);
    Field<3> field3(7, 8, 9);

    Mesh<1> mesh1(0, 1, 11);
    Mesh<2> mesh2(0, 0, 1, 2, 11, 21);
    Mesh<3> mesh3(0, 0, 0, 1, 2, 3, 11, 21, 31);
    return 0;
}
