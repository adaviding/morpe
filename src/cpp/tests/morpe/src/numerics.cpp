#include "test_morpe.h"

using namespace morpe::numerics;

TEST(IX_tests, pascal)
{
    EXPECT_EQ(2,   IX::pascal(1,1));

    EXPECT_EQ(3,   IX::pascal(1,2));
    EXPECT_EQ(3,   IX::pascal(2,1));

    EXPECT_EQ(6,   IX::pascal(2,2));

    EXPECT_EQ(15,  IX::pascal(2,4));
    EXPECT_EQ(15,  IX::pascal(4,2));

    EXPECT_EQ(210, IX::pascal(4,6));
    EXPECT_EQ(210, IX::pascal(6,4));

    EXPECT_EQ(495, IX::pascal(4,8));
    EXPECT_EQ(495, IX::pascal(8,4));
}
