#pragma once

#include "_internal.h"
#include "sal.h"
#include "err/stop_error.h"

#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>

#include <ostream>
#include <stop_token>

namespace morpe { namespace err
{
    /// The type for stack trace info (as per boost examples).
    typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> boost_trace_info;

    /// Throw an exception with the stack trace attached.  This preserves an existing stack trace if one is already attached.
    /// @param ex The exception.
    template<class E>
    void chuck(
            _In_ const E &ex)
    {
        if (nullptr == boost::get_error_info<boost_trace_info>(ex))
        {
            throw boost::enable_error_info(ex) << boost_trace_info(boost::stacktrace::stacktrace());
        }
        else
        {
            throw ex;
        }
    };

    /// Throw an exception with the stack trace attached.  This overwrites any existing stack trace.
    /// @param ex The exception.
    template<class E>
    void rechuck(
            _In_ const E &ex)
    {
        throw boost::enable_error_info(ex) << boost_trace_info(boost::stacktrace::stacktrace());
    };

    /// Throws a #stop_error if a stop was requested.
    void throw_if_stopped(
            std::stop_token stop_token)
    {
        if (stop_token.stop_requested())
        {
            chuck(stop_error());
        }
    };

    /// Write the stack trace for the given exception.
    /// @param ostream The output stream on which to write.
    /// @param ex The exception.
    template<class E>
    void write_trace(
            _Inout_ std::ostream &ostream,
            _In_    const E &ex)
    {
        const boost::stacktrace::stacktrace* st = boost::get_error_info<boost_trace_info>(ex);
        if (st == nullptr)
        {
            ostream << "[no stack trace]";
        }
        else
        {
            ostream << *st;
        }
    };

    /// Write the message and stack trace for the given exception.
    /// @param ostream The output stream on which to write.
    /// @param ex The exception.
    template<class E>
    void write_exception(
            _Inout_ std::ostream &ostream,
            _In_    const E &ex) {
        ostream << ex.what();

        const boost::stacktrace::stacktrace* st = boost::get_error_info<boost_trace_info>(ex);
        if (st != nullptr)
        {
            ostream << std::endl << *st;

        }

        write_trace(ostream, ex);
    };
}}

/// Throw an exception (with stack trace) if the expression is false.
/// @param expr Any valid C++ expression.
#define ChkTrue(expr)                                                                           \
    if (BOOST_UNLIKELY(!(expr)))                                                                \
    {                                                                                           \
        ::morpe::err::rechuck(std::runtime_error("False: " TOSTRING(expr)));                    \
    }

/// Throw an exception (with stack trace) if the expression is false, and include a message.
/// @param expr Any valid C++ expression.
/// @param message  An additional message to provide in the exception.
#define ChkTrueMsg(expr, message)                                                               \
    if (BOOST_UNLIKELY(!(expr)))                                                                \
    {                                                                                           \
        ::morpe::err::rechuck(std::runtime_error(message "  False: " TOSTRING(expr)));          \
    }

/// Throw an exception (with stack trace) if the expression is true.
/// @param expr Any valid C++ expression.
#define ThrowIf(expr)                                                                           \
    if (BOOST_UNLIKELY(expr))                                                                   \
    {                                                                                           \
        ::morpe::err::rechuck(std::runtime_error("True: " TOSTRING(expr)));                     \
    }

/// Throw an exception (with stack trace) if the expression is true, and include a message.
/// @param expr Any valid C++ expression.
/// @param message  An additional message to provide in the exception.
#define ThrowIfMsg(expr, message)                                                               \
    if (BOOST_UNLIKELY(expr))                                                                   \
    {                                                                                           \
        ::morpe::err::rechuck(std::runtime_error(message "  True: " TOSTRING(expr)));           \
    }
