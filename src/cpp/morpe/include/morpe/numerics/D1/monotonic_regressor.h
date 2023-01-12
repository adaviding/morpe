#pragma once

#include "../../_internal.h"
#include "monotonic_regression_type.h"

#include <mutex>
#include <stop_token>
#include <vector>

namespace morpe { namespace numerics { namespace D1
{
    /// A class which encapsulates monotonic regression.
    ///
    /// This just exposes a single #run method.  The class is needed to recycle heap resources (for efficiency)
    /// in a scenario where the caller wants to execute #run repeatedly.
    class monotonic_regressor
    {
    public:
        int32_t run(
                _In_    std::stop_token stop_token,
                _In_    monotonic_regression_type type,
                _In_    std::vector<double>& input,
                _Inout_ std::vector<double>& output);

    private:

        /// Intermediate storage for the derivative, and it gets altered as the original unfolds.
        std::vector<double> derivative;

        /// A mutex that ensures that #run cannot be executed concurrently since each execution relies
        /// on some arrays which are allocated only once for the sake of efficiency.
        std::recursive_mutex mutex;

        /// Intermediate storage for an intermediate result:  The non-decreasing function.
        std::vector<double> non_decreasing_values;

        static int32_t annihilate_decreasing_energy(
                _In_    std::stop_token stop_token,
                _Inout_ double& dy_min,
                _In_    double dy_min_threshold,
                _In_    int32_t limit_num_trips_through_loop,
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy);

        static double calculate_blend_factor(
                _In_    const std::vector<double>& non_decreasing_values,
                _In_    int32_t length);

        static void convert_non_decreasing_to_increasing(
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy,
                _In_    double ymean,
                _In_    double ymin,
                _In_    double ymax,
                _In_    double dymin);

        static void eliminate_decreasing_energy(
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy,
                _In_    double ymean,
                _In_    double ymin,
                _In_    double ymax);

        static void repel_decreasing_energy_from_index(
                _In_    int32_t idx,
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy);

        static double update_derivative_and_find_min(
                _In_    std::vector<double>& y,
                _Inout_ std::vector<double>& dy);
    };
}}}
