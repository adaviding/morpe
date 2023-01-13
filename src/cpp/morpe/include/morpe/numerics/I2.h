// General numerical types and methods related to signed 32-bit integers in 2 spatial dimensions.

#pragma once

#include "../_internal.h"

namespace morpe { namespace numerics { namespace I2
{
    /// Returns the entry (a,b) of the Pascal matrix.  The parameters are interchangeable because the matrix is symmetric.
    /// The value returned is equal to:  (a+b)! / a! / b!
    /// @param a The row of the Pascal matrix.
    /// @param b The column of the Pascal matrix.
    int pascal(int a, int b);
}}}
