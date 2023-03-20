#include </usr/include/gtest/gtest.h>
//#include <gtest.h>
#include <type_traits>

#include "../../include/Field.h"

#include <vector>

template <typename T>
class FieldTest : public ::testing::Test
{
public:
protected:
    Field<T::value> *testField;

    FieldTest();
    virtual ~FieldTest() {}

    void SetUp();
    void TearDown();
private:
    size_t m_sizeX;
    size_t m_sizeY;
    size_t m_sizeZ;
};

using test_types = ::testing::Types<
    std::integral_constant<size_t, 1>,
    std::integral_constant<size_t, 2>,
    std::integral_constant<size_t, 3>>;