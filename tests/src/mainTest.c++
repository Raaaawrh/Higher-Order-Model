#include </usr/include/gtest/gtest.h>
//#include <gtest.h>
//#include <gmock.h>

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::InitGoogleMock(&argc, argv);

    return RUN_ALL_TESTS();
}
