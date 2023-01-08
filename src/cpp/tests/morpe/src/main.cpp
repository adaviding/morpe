#include "test_morpe.h"

int main(int nargs, char** args)
{
    ::testing::InitGoogleTest(&nargs, args);
    return RUN_ALL_TESTS();
}
