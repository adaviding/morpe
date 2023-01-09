#include "test_morpe.h"

using namespace morpe::numerics;

TEST(I2_tests, pascal)
{
    EXPECT_EQ(2,   I2::pascal(1,1));

    EXPECT_EQ(3,   I2::pascal(1,2));
    EXPECT_EQ(3,   I2::pascal(2,1));

    EXPECT_EQ(6,   I2::pascal(2,2));

    EXPECT_EQ(15,  I2::pascal(2,4));
    EXPECT_EQ(15,  I2::pascal(4,2));

    EXPECT_EQ(210, I2::pascal(4,6));
    EXPECT_EQ(210, I2::pascal(6,4));

    EXPECT_EQ(495, I2::pascal(4,8));
    EXPECT_EQ(495, I2::pascal(8,4));
}
