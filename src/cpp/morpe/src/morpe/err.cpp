#include "morpe.h"

namespace morpe { namespace err
{
    void throw_if_stopped(
            std::stop_token stop_token)
    {
        if (stop_token.stop_requested())
        {
            chuck(stop_error());
        }
    }
}}
