#include "../include/FieldTest.h"

template <typename T>
FieldTest<T>::FieldTest()
{
    m_sizeX = 10;
    m_sizeY = 4;
    m_sizeZ = 17;
}

template <typename T>
void FieldTest<T>::SetUp()
{
    switch (T::value)
    {
    case size_t(1):
        testField = new Field<T::value>(m_sizeX);
        break;
    case size_t(2):
        testField = new Field<T::value>(m_sizeX, m_sizeY);
        break;
    case size_t(3):
        testField = new Field<T::value>(m_sizeX, m_sizeY, m_sizeZ);
        break;
    default:
        break;
    }
}

template <typename T>
void FieldTest<T>::TearDown()
{
    delete testField;
}

TYPED_TEST_CASE(FieldTest, test_types);

TYPED_TEST(FieldTest, fieldConstructorTest)
{
    size_t sizeX, sizeY, sizeZ;

    switch (T::value)
    {
    case size_t(1)
        //sizeX = testField.getX();
        //ASSERT_EQ(sizeX, m_sizeX);
        break;
    
    default:
        break;
    }
}