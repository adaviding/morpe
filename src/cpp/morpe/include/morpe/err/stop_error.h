#pragma once

#include "../_internal.h"

#include <stdexcept>

namespace morpe { namespace err
{
    /// This error is typically thrown when a stop is requested.
    class stop_error : public std::runtime_error
    {
    public:
        stop_error()
            : std::runtime_error("Operation canceled:  A stop was requested.")
        {
        };
    };
}}
