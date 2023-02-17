#include "morpe.h"

namespace morpe { namespace numerics { namespace D1
{
    // --------------------------------
    // public functions
    // --------------------------------

    int32_t monotonic_regressor::run(
            _In_    std::stop_token stop_token,
            _In_    monotonic_regression_type type,
            _In_    const std::vector<double>& input,
            _Inout_ std::vector<double>& output)
    {
        ThrowIf(input.size() < 2);

        // Count the number of trips through an inner loop.
        int32_t num_trips_through_loop = 0;
        int limit_num_trips_through_loop = 100 + (int)(30.0 * std::log((double)input.size()));

        // Some shortcut varaibles.
        int n = input.size();

        // Ensure the output is the same size as the input.
        output.resize(n);

        // Prevent concurrent execution.  We are using member variables to store intermediate answers,
        // and concurrency would create catastrophic interference.
        std::unique_lock<std::recursive_mutex> lock(this->mutex);

        // Ensure that intermediate memory is allocated.
        if (this->derivative.size() < n)
        {
            // Allocate a bit extra to reduce unnecessary allocations in the future.
            this->derivative.resize((size_t)(n * 1.2 + 0.5));
            this->non_decreasing_values.resize((size_t)(n * 1.2 + 0.5));
        }

        int i, j;

        //    Compute min, max, mean, etc.
        double ymin = input[0];
        double ymax = input[0];
        double ymean = input[0];
        double dy_min = input[1] - input[0];
        output[0] = input[0];
        for (i = 1; i < n; i++)
        {
            if (input[i] < ymin) ymin = input[i];
            if (input[i] > ymax) ymax = input[i];
            ymean += input[i];
            double dy_this = input[i] - input[i - 1];
            this->derivative[i - 1] = dy_this;
            if (dy_this < dy_min)
                dy_min = dy_this;
            output[i] = input[i];
        }

        ymean /= (double)n;
        double yrng = ymax - ymin;
        double dy_min_threshold = -yrng * 0.001;

        if (dy_min >= 0.0)
        {
            // Check to see if we are done.
            if (dy_min > 0.0
                || dy_min == 0.0 && type == monotonic_regression_type::non_decreasing)
            {
                // Let's say it was 0 trips through a loop.
                return 0;
            }
        }
        else
        {
            // Exhaustively annihilate decreasing energy in an unbiased way.  This is the best approach to achieving
            // a non-decreasing function, but there is no guarantee that it will work.  Even if it doesn't work, it
            // typically gets us most of the way home.
            num_trips_through_loop = annihilate_decreasing_energy(
                    stop_token,
                    dy_min,
                    dy_min_threshold,
                    limit_num_trips_through_loop,
                    output,
                    this->derivative);

            // This will ensure that 'output' is non-decreasing, and some other nice properties.
            eliminate_decreasing_energy(
                    output,
                    this->derivative,
                    ymean,
                    ymin,
                    ymax);
        }

        // Throw if stop requested.
        err::throw_if_stopped(stop_token);

        if (type == monotonic_regression_type::blended)
        {
            // Save the non-decreasing values.  We will need them later.
            this->non_decreasing_values.assign(output.begin(), output.end());
        }

        // Change a non-decreasing function into a monotonic function.
        if (type == monotonic_regression_type::increasing || type == monotonic_regression_type::blended)
        {
            //    Compute min, max, and recompute derivative.
            dy_min = update_derivative_and_find_min(output, this->derivative);

            convert_non_decreasing_to_increasing(
                    output,
                    this->derivative,
                    ymean,
                    ymin,
                    ymax,
                    dy_min);

            // Throw if stop requested.
            err::throw_if_stopped(stop_token);

            if (type == monotonic_regression_type::blended)
            {
                // Blend the non-decreasing and increasing functions together.

                // The proportion of weight given to the increasing function.
                double pw_increasing = calculate_blend_factor(this->non_decreasing_values, output.size());
                ThrowIf(pw_increasing < 0.0 || pw_increasing > 1.0);

                // The proportion of weight given to the non-decreasing function.
                double pw_non_decreasing = 1.0 - pw_increasing;

                // Blend
                for (i = 0; i < n; i++)
                {
                    output[i] = pw_increasing * output[i] + pw_non_decreasing * non_decreasing_values[i];
                }
            }
        }

        // Defend against a degenerate case.
        if (std::isnan(output[0]))
        {
            // In the degenerate case where the output should be flat, we might have some NaN and really it should just be flat.
            for (i = 0; i < n; i++)
                output[i] = ymean;
        }

        return num_trips_through_loop;
    }

    // --------------------------------
    // private static functions
    // --------------------------------

    int32_t monotonic_regressor::annihilate_decreasing_energy(
            _In_    std::stop_token stop_token,
            _Inout_ double& dy_min,
            _In_    double dy_min_threshold,
            _In_    int32_t limit_num_trips_through_loop,
            _Inout_ std::vector<double>& y,
            _Inout_ std::vector<double>& dy)
    {
        // The number of trips through the loop.
        int32_t output = 0;

        int32_t nm1 = y.size() - 1;

        while (dy_min < dy_min_threshold)
        {
            int32_t i = 0;      // Point to beginning of dy
            int32_t j = y.size() - 2;  // Point to end of dy.

            //    Go through all of dy
            while (j > 0)
            {
                // Check for cancellation.
                err::throw_if_stopped(stop_token);

                // Make an adjustment at point 'i'.
                repel_decreasing_energy_from_index(i, y, dy);

                // Make an adjustment at point 'j'.
                if (j != i)
                {
                    repel_decreasing_energy_from_index(j, y, dy);
                }

                j--;
                i++;
            }

            // Update minimum of dy
            dy_min = dy[0];
            for (i = 1; i < nm1; i++)
                if (dy[i] < dy_min)
                    dy_min = dy[i];

            // Increment number of trips.
            output++;
            if (output >= limit_num_trips_through_loop)
                break;
        }

        return output;
    }

    double monotonic_regressor::calculate_blend_factor(
            _In_    const std::vector<double>& non_decreasing_values,
            _In_    int32_t length)
    {
        // The blend factor is a number in the range (0, 1) representing the proportion of weight given to the "increasing" function.
        double output = 1.0;

        // Measure the size of the largest non-decreasing span.
        int run_len = 0;
        int max_run_len = 0;

        for (int i = 1; i < length; i++)
        {
            // Measure the size of the largest non-decreasing span.
            if (non_decreasing_values[i] != non_decreasing_values[i - 1])
                run_len = 0;
            else
                run_len++;
            if (run_len > max_run_len)
                max_run_len = run_len;
        }

        // Compute the blend factor
        output = 1.0 / (1.0 + std::pow((double)(4 * max_run_len) / (double)(length - 1), 3.0));

        return output;
    }

    void monotonic_regressor::convert_non_decreasing_to_increasing(
            _Inout_ std::vector<double>& y,
            _Inout_ std::vector<double>& dy,
            _In_    double ymean,
            _In_    double ymin,
            _In_    double ymax,
            _In_    double dymin)
    {
        int n = y.size();
        int nm1 = n-1;

        double zero_floor = (y[nm1] - y[0]) * 0.00001 / n;

        //    Only continue if the minimum derivative is zero (or tiny) and if the function is not totally flat.
        if (y[nm1] > y[0] && dymin <= zero_floor)
        {
            //----------------------------------------------------------------------------------
            // First we alter the derivative to be positive everywhere.
            //----------------------------------------------------------------------------------

            // Init indexes at a location where [iLow, iMid] straddles the first flat region.
            int i_low, i_mid = 0, i_high;
            //    iMid:  The index of a non-flat spot
            //    iLow:  The preceding non-flat index before iMid
            //    iHigh: The subsequent non-flat index following iMid
            i_low = -1;
            while (i_mid < nm1 && dy[i_mid] <= zero_floor) i_mid++;


            // Advance from iMid's initial value up through iMid = nm-1.
            while (i_mid < nm1)
            {
                // Advance index of iHigh to the subsequent non-flat index following iMid
                i_high = i_mid + 1;
                while (i_high < nm1 && dy[i_high] <= zero_floor) i_high++;

                if (i_high - i_mid > 1)
                {
                    if (i_mid - i_low > 1)
                    {
                        //---------------------------
                        //    Distribute mass to both sides
                        //---------------------------
                        //    The mass to each index (center index counts twice and receives double mass).
                        double dy_this = dy[i_mid] / (i_high - i_low);
                        //    Mass to iMid (double mass)
                        dy[i_mid] = dy_this * 2.0;
                        //    Mass below iMid
                        for (int i = i_low + 1; i < i_mid; i++)
                            dy[i] += dy_this;
                        //    Mass above iMid
                        for (int i = i_mid + 1; i < i_high; i++)
                            dy[i] += dy_this;
                        //---------------------------
                    }
                    else
                    {
                        //---------------------------
                        //    Distribute mass to high side
                        //---------------------------
                        //    The mass to each index (center index will get single mass).
                        double dy_this = dy[i_mid] / (i_high - i_mid);
                        dy[i_mid] = dy_this;
                        //    Mass from iMid to iHigh.
                        for (int i = i_mid + 1; i < i_high; i++)
                            dy[i] += dy_this;
                        //---------------------------
                    }

                }
                else if (i_mid - i_low > 1)
                {
                    //---------------------------
                    //    Distribute mass to low side
                    //---------------------------
                    //    The mass to each index (center index will get single mass).
                    double dy_this = dy[i_mid] / (i_mid - i_low);
                    //    Mass to iMid
                    dy[i_mid] = dy_this;
                    //    Mass from iLow to iMid (must be added to mass that may already exist there).
                    for (int i = i_low + 1; i < i_mid; i++)
                        dy[i] += dy_this;
                    //---------------------------
                }

                //    Advance indexes for next trip
                i_low = i_mid; // The "subseqent flat region" on the next trip is set to the "preceding flat region" from the last trip.
                i_mid = i_high; // We center on last trip's high index.
            }

            //----------------------------------------------------------------------------------
            // Next we integrate the derivative while ensuring the mean of the original input.
            //----------------------------------------------------------------------------------
            double energyBelow = ymean - y[0];
            double energyAbove = 0.0;
            for (int i = 0; i < nm1;)
            {
                y[i + 1] = y[i] + dy[i];
                if (y[++i] > ymean)
                    energyAbove += y[i] - ymean;
                else
                    energyBelow += ymean - y[i];
            }

            // Prepare to scale the values above and below differently.
            double middleEnergy = (energyAbove + energyBelow) / 2.0;
            double scalarBelow = middleEnergy / energyBelow;
            double scalarAbove = middleEnergy / energyAbove;

            // Limit the size of the scalars to ensure that the span of y-values is correct.
            double carry = 1.0;
            if (scalarBelow > 1.0)
                carry = (ymean - ymin) / (ymean - y[0]) / scalarBelow;
            else if (scalarAbove > 1.0)
                carry = (ymax - ymean) / (y[nm1] - ymean) / scalarAbove;
            if (carry < 1.0)
            {
                scalarBelow *= carry;
                scalarAbove *= carry;
            }

            // Adjust elements on big and small sides so that mean is preserved.
            for (int i = 0; i < n; i++)
                if (y[i] < ymean)
                    y[i] = (y[i] - ymean) * scalarBelow + ymean;
                else if (y[i] > ymean)
                    y[i] = (y[i] - ymean) * scalarAbove + ymean;
        }
    }

    void monotonic_regressor::eliminate_decreasing_energy(
            _Inout_ std::vector<double>& y,
            _Inout_ std::vector<double>& dy,
            _In_    double ymean,
            _In_    double ymin,
            _In_    double ymax)
    {
        int n = y.size();
        int nm1 = n - 1;

        // Fix small dy values to be non-negative (so that "y" is non-decreasing).
        // This operation will bias 'y' in an undesirable way:  The function will increase too rapidly, so there is a "clean-up" procedure after.
        double a_carry = 0.0;
        for (int i = 0; i < nm1;)
        {
            if (dy[i++] < 0.0)
                a_carry -= dy[i - 1];

            // Set 'y'
            y[i] += a_carry;

            // Watch out for tiny errors with floating point precision.  Ensure it is really non-decreasing.
            if (y[i] < y[i - 1])
            {
                a_carry += y[i - 1] - y[i];
                y[i] = y[i - 1];
            }
        }

        // The sum of deviations (from the mean) for values above the mean.
        double sum_above_mean = 0.0;

        // The sum of deviations (from the mean) for values below the mean.
        double sum_below_mean = 0.0;

        // The above operation may have (1) expanded the range of our function and (2) increased the mean of the function.
        // Let's measure this so we can fix it.

        double y0 = y[0];
        double yrng = y[nm1] - y0;

        // We can multiply each y-value with this scalar to eliminate the range expansion.
        double scalar = (yrng - a_carry) / yrng;

        for (int i = 0; i < n; i++)
        {
            // Counter the expansion of our function's range
            double yi = y[i] = (y[i] - y0) * scalar + y0;

            // Track what happened to the mean of the function
            if (yi < ymean)
                sum_below_mean += ymean - yi;
            else
                sum_above_mean += yi - ymean;
        }

        // See if we changed the mean.
        if (sum_above_mean != sum_below_mean)
        {
            // Yes it changed, we must fix it.
            if (sum_below_mean == 0.0 || sum_above_mean == 0.0)
            {
                // This is a degenerate situation, like if the function is flat.
                if (sum_below_mean > 0.0)
                {
                    double adj = sum_below_mean / (2 * n); // The extent to which the actual mean is below the desired mean.
                    for (int i = 0; i < n; i++)
                        y[i] += adj;
                }
                else
                {
                    double adj = sum_above_mean / (2 * n); // The extent to which the actual mean is above the desired mean.
                    for (int i = 0; i < n; i++)
                        y[i] -= adj;
                }
            }
            else
            {
                // The first constraint:  The distance of all values from the mean will be scaled so that sum_below_mean == sum_above_mean
                // (i.e. so that the mean is unaltered).
                double target = (sum_above_mean + sum_below_mean) / 2.0;
                double adj_below = target / sum_below_mean;
                double adj_above = target / sum_above_mean;

                // The second constraint:  Ensure that the range does not exceed the original range.
                double adj = 1.0;
                if (adj_below > 1.0)
                    adj = (ymean - ymin) / (ymean - y[0]) / adj_below;
                else if (adj_above > 1.0)
                    adj = (ymax - ymean) / (y[nm1] - ymean) / adj_above;
                if (adj < 1.0)
                {
                    adj_below *= adj;
                    adj_above *= adj;
                }

                for (int i = 0; i < n; i++)
                {
                    double yi = y[i];

                    if (yi < ymean)
                        y[i] = (yi - ymean) * adj_below + ymean;
                    else if (yi > ymean)
                        y[i] = (yi - ymean) * adj_above + ymean;
                }
            }
        }
    }

    void monotonic_regressor::repel_decreasing_energy_from_index(
            _In_    int32_t idx,
            _Inout_ std::vector<double>& y,
            _Inout_ std::vector<double>& dy)
    {
        // Perform this operation on point 'i'.
        double half_dy = dy[idx];
        if (half_dy < 0.0)
        {
            // Make it half.
            half_dy /= 2.0;

            //    Element y(i+1) must be incremented and y(i) must be decremented.
            y[idx + 1] -= half_dy; // increment y(i+1)
            y[idx] += half_dy; // decrement y(i)

            //    Now, 'dy' must be adjusted to account for this new change.
            if (idx > 0)
                dy[idx - 1] += half_dy; // The difference between y(i-1) and y(i) goes down.

            if (idx < y.size() - 2)
                dy[idx + 1] += half_dy; // The difference between y(i+1) and y(i+2) goes down.

            dy[idx] = 0.0; // The difference between y(i) and y(i+1) is now zero.
        }
    }

    double monotonic_regressor::update_derivative_and_find_min(
            _In_    std::vector<double>& y,
            _Inout_ std::vector<double>& dy)
    {
        double output = std::numeric_limits<double>::max();

        for (int i = 1; i < y.size(); i++)
        {
            double diff = y[i] = - y[i-1];
            if (diff < output)
            {
                output = diff;
            }

            dy[i-1] = diff;
        }

        return output;
    }
}}}
