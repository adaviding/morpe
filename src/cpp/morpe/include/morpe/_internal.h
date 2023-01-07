#pragma once

// This header file is special:  It gets pulled early into almost every other header file in this project.

//  Helper functions for macros
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

//  This is used for logging and debugging.  Shows filename and line.  Can be passed to a function as string.
#define THIS_LINE __FILE__ ":" TOSTRING(__LINE__)

//  Shared library export and import annotations.  For example, when consuming a DLL:
//
//    #include "morpe/_first_header_file.h"
//    #define MORPE_API MORPE_IMPORT
//    #include "morpe.h"
#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__ || defined __clang__
    #ifdef __GNUC__
        #define MORPE_IMPORT __attribute__((dllimport))
        #define MORPE_EXPORT __attribute__((dllexport))
    #else
        #define MORPE_IMPORT __declspec(dllimport)
        #define MORPE_EXPORT __declspec(dllexport)
    #endif
    #define MORPE_HIDE
#else
    #define MORPE_EXPORT __attribute__((visibility ("default")))
    #define MORPE_IMPORT MORPE_EXPORT
    #define MORPE_HIDE   __attribute__((visibility ("hidden")))
#endif

// By default, the header files are used to build the shared library (not consume it).
#define MORPE_API MORPE_EXPORT
