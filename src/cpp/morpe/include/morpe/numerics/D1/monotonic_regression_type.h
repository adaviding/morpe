#pragma once

#include "../../_internal.h"

namespace morpe { namespace numerics { namespace D1
{
    /// The type of monotonic regression.  This describes the quality of the output.
    enum class monotonic_regression_type
    {
        /// This is the recommended setting.  It is a simple blending of two other methods (#non_decreasing, #increasing).
        /// The output is monotonic.
        blended,

        /// This is like #non_decreasing, but flatness is removed by linear interpolation (through any flat spot).
        /// The output is monotonic.
        increasing,

        /// The output is non-decreasing, but it may have "flat" spots (where the slope is 0), in which case it is not
        /// truly monotonic.
        non_decreasing,
    };
}}}
