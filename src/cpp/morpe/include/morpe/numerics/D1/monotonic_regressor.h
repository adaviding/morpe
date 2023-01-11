#pragma once

#include "../../_internal.h"
#include "monotonic_regression_type.h"

#include <mutex>
#include <stop_token>
#include <vector>

namespace morpe { namespace numerics { namespace D1
{
    /// A class which encapsulates monotonic regression.
    class monotonic_regressor
    {
    public:
        int32_t run(
                _In_    std::stop_token stop_token,
                _In_    monotonic_regression_type type,
                _In_    std::vector<double>& input,
                _Inout_ std::vector<double>& output);

    private:

        /// Intermediate storage for the derivative and the altered derivative.
        std::vector<double> dy;

        /// A mutex that ensures that 'run' cannot be executed concurrently..
        std::recursive_mutex mutex;

        /// Intermediate storage for an intermediate result (the non-decreasing function).
        std::vector<double> ynd;

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
