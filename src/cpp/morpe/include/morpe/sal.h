#pragma once

#include "_internal.h"

#if defined_WIN32 || defined _WIN64
#   include <sal.h>          // from windows
#else
#   include "mingw64/sal.h"  // from mingw64
#endif
