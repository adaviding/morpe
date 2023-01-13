#include "morpe.h"

namespace morpe { namespace numerics { namespace I2
{
    int pascal(int a, int b)
    {
        int32_t aa = std::max(a, b);
        int32_t bb = std::min(a, b);
        int32_t c = a + b;
        int32_t output = 1;
        int32_t i;

        for (i = aa + 1; i <= c; i++)
            output *= i;

        for (i = 2; i <= bb; i++)
            output /= i;

        return output;
    }
}}}
