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
        /// Performs a monotonic regression of a tabulated function.  This method calculates a monotonic function that is approximately equal
        /// to the tabulated function provided.
        /// @param stop_token If triggered, a #stop_error will be thrown quickly.
        /// @param type The type of monotonic regression.
        /// @param input The tabulated function.
        /// @param output The monotonic function.  This should be supplied as a vector of length identical to input.
        /// The following will be guaranteed:
        ///     mean(output) == mean(input)
        ///     min(output)  >= min(input)
        ///     max(output)  <= max(input)
        /// @return The number of trips through an inner loop.
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

        /// Repeatedly spread decreasing energy to immediate neighbors until a threshold is met, or until we have gone through the
        /// loop too many times.  Decreasing energy is annihilated by increasing energy.  Thus, it tends to be annihilated by this
        /// spreading action.
        ///
        /// This is the best approach to annihilating the decreasing energy, but this method may not eliminate all of it.  Even so,
        /// it typically gets us most of the way home.
        /// @param stop_token If triggered, a #stop_error will be thrown quickly.
        /// @param dy_min The minimum value of 'dy'.  This value will be updated with each trip through the loop.
        /// @param dy_min_threshold The loop may exit when 'dy_min' exceeds this value.  This should be a small negative number.
        /// @param limit_num_trips_through_loop The maximum number of trips through the loop.
        /// @param y The tabulated values of the function 'y'.
        /// @param dy The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        /// @return The number of trips taken through the loop.
        static int32_t annihilate_decreasing_energy(
                _In_    std::stop_token stop_token,
                _Inout_ double& dy_min,
                _In_    double dy_min_threshold,
                _In_    int32_t limit_num_trips_through_loop,
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy);

        /// This is needed in the context of calculating a #monotonic_regression_type::blended monotonic regression,
        /// an operation basically blends the outputs of the other operation types.
        /// * #monotonic_regression_type::non_decreasing
        /// * #monotonic_regression_type::increasing
        ///
        /// This method calculates a "blend factor", or basically a weighting used to mix the two other output types.
        /// @param non_decreasing_values The monotonic regression output for #monotonic_regression_type::non_decreasing.
        /// @param length The actual number of non-decreasing values, since the vector passed in may have extra padding.
        /// @return The proportion of weight given to the #monotonic_regression_type::increasing values.
        static double calculate_blend_factor(
                _In_    const std::vector<double>& non_decreasing_values,
                _In_    int32_t length);

        /// Modifies the non-decreasing values of 'y' to be strictly increasing, so that each value of 'y' is greater than (not equal to) the prior value.
        /// @param y The tabulated values of the function 'y'.  On input, we know that these values are non-decreasing, but not necessarily
        ///     increasing.  On output, these values will be modified and will contain the increasing values.
        /// @param dy The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].  These values will be modified.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        ///
        ///     On exit, these values will be meaningless:  They are not kept in sync with 'y'.
        /// @param ymean The mean of 'y'.  Although we will change the 'y' values, we will ensure that it still has this mean.
        /// @param ymin The minimum allowed value of 'y'.  No value of 'y' will be below this minimum.
        /// @param ymax The maximum allowed value of 'y'.  No value of 'y' will be above this maximum.
        /// @param dymin The current minimum value of the 'dy' vector.
        static void convert_non_decreasing_to_increasing(
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy,
                _In_    double ymean,
                _In_    double ymin,
                _In_    double ymax,
                _In_    double dymin);

        /// This is executed after <see cref="AnnihilateDecreasingEnergy"/> to ensure that all decreasing energy has been finally eliminated.  This method is
        /// not the best way to annihilate the decreasing energy, but it does guarantee the following:
        /// * 'y' will be non-decreasing.
        /// * 'dy' will be non-negative.
        /// * The mean of 'y' will be 'yMean'.
        /// * The min of 'y' shall not be below 'yMin'.
        /// * The max of 'y' shall not be above 'yMax'.
        /// @param y The tabulated values of the function 'y'.
        /// @param dy The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        /// @param ymean The mean of 'y'.  Although we will change the 'y' values, we will ensure that it still has this mean.
        /// @param ymin The minimum allowed value of 'y'.  No value of 'y' will be below this minimum.
        /// @param ymax The maximum allowed value of 'y'.  No value of 'y' will be above this maximum.
        static void eliminate_decreasing_energy(
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy,
                _In_    double ymean,
                _In_    double ymin,
                _In_    double ymax);

        /// If decreasing energy is present at the given index, then it is moved to the adjacent indices (idx-1, idx+1).
        /// If either of the adjacent indices has increasing energy, it will annihilate the decreasing energy in equal
        /// measure.
        ///
        /// @idx The zero-based index into 'y' where the adjustment is to be made.
        /// @param y The tabulated values of the function 'y'.
        /// @param dy The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        static void repel_decreasing_energy_from_index(
                _In_    int32_t idx,
                _Inout_ std::vector<double>& y,
                _Inout_ std::vector<double>& dy);

        /// Updates the values of 'dy' given the current values of 'y' and also finds the minimum of 'dy'.
        /// @param y The current values of 'y'.
        /// @param dy A container which holds the output values of 'dy'.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        /// @return The minimum value of 'dy'.
        static double update_derivative_and_find_min(
                _In_    std::vector<double>& y,
                _Inout_ std::vector<double>& dy);
    };
}}}
