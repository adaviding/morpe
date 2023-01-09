#include "morpe.h"

namespace morpe { namespace numerics { namespace D1
{
    /// PRIVATE STATIC
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

    /// PRIVATE STATIC
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

    /// PRIVATE STATIC
    /// If decreasing energy is present at the given index, then it is moved to the adjacent indices (idx-1, idx+1).
    /// If either of the adjacent indices has increasing energy, it will annihilate the decreasing energy in equal
    /// measure.
    ///
    /// @idx The zero-based index into 'y' where the adjustment is to be made.
    /// @param y The tabulated values of the function 'y'.
    /// @param dy The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].
    ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
    ///     can be ignored.
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
    int32_t monotonic_regressor::run(
            _In_    std::stop_token stop_token,
            _In_    monotonic_regression_type type,
            _In_    std::vector<double>& input,
            _Inout_ std::vector<double>& output)
    {
        ThrowIf(input.size() < 2);

        // Count the number of trips through an inner loop.
        int32_t num_trips_through_loop = 0;
        int limit_num_trips_through_loop = 100 + (int)(30.0 * std::log((double)input.size()));

        // Some shortcut varaibles.
        int n = input.size();
        int nm1 = input.size() - 1;

        // Ensure the output is the same size as the input.
        output.resize(n);

        // Prevent concurrent execution.  We are using member variables to store intermediate answers,
        // and concurrency would create catastrophic interference.
        std::unique_lock<std::recursive_mutex> lock(this->mutex);

        // Ensure that intermediate memory is allocated.
        if (this->dy.size() < n)
        {
            // Allocate a bit extra to reduce unnecessary allocations in the future.
            this->dy.resize((size_t)(n * 1.2 + 0.5));
            this->ynd.resize((size_t)(n * 1.2 + 0.5));
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
            this->dy[i - 1] = dy_this;
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
                    this->dy);

            // This will ensure that 'output' is non-decreasing, and some other nice properties.
            eliminate_decreasing_energy(
                    output,
                    this->dy,
                    ymean,
                    ymin,
                    ymax);
        }

        // The blend factor is a number in the range (0, 1) representing the proportion of weight given to the "increasing" function.
        double proportionBlendIncreasing = 1.0;
        if (type == monotonic_regression_type::blended)
        {
            // Measure the size of the largest non-decreasing span.
            int runLen = 0;
            int maxRunLen = 0;

            this->ynd[0] = output[0];     // Save the non-decreasing values
            for (i = 1; i < n; i++)
            {
                this->ynd[i] = output[i]; // Save the non-decreasing values

                // Measure the size of the largest non-decreasing span.
                if (output[i] != output[i - 1])
                    runLen = 0;
                else
                    runLen++;
                if (runLen > maxRunLen)
                    maxRunLen = runLen;
            }

            // Compute the blend factor
            proportionBlendIncreasing = 1.0 / (1.0 + std::pow((double)(4 * maxRunLen) / (double)(n - 1), 3.0));
        }

        // Change a non-decreasing function into a monotonic function.
        if (type == monotonic_regression_type::increasing || type == monotonic_regression_type::blended)
        {
            //    Compute min, max, and recompute derivative.
            dy_min = output[1] - output[0];
            for (i = 1; i < n; i++)
            {
                dy_this = output[i] - output[i - 1];
                this->dy[i - 1] = dy_this;
                if (dy_this < dy_min)
                    dy_min = dy_this;
            }

            sum_small = (output[nm1] - output[0]) * 0.00001 / (double)n;

            //    Only continue if the minimum derivative is zero (or tiny) and if the function is not totally flat.
            if (output[nm1] > output[0] && dy_min <= sum_small)
            {
                //----------------------------------------------------------------------------------
                //    Alter the derivative to be positive everywhere, then integrate the altered derivative.
                //----------------------------------------------------------------------------------

                //---------------------------
                //    Init indexes at a location where [iLow, iMid] straddles the first flat region.
                //---------------------------
                int iLow, iMid = 0, iHigh;
                //    iMid:  The index of a non-flat spot
                //    iLow:  The preceding non-flat index before iMid
                //    iHigh: The subsequent non-flat index following iMid
                iLow = -1;
                while (iMid < nm1 && this->dy[iMid] <= sum_small) iMid++;
                //---------------------------

                //    Advance from iMid's initial value up through iMid = nm-1.
                while (iMid < nm1)
                {
                    //    Advance index of iHigh to the subsequent non-flat index following iMid
                    iHigh = iMid + 1;
                    while (iHigh < nm1 && this->dy[iHigh] <= sum_small) iHigh++;

                    if (iHigh - iMid > 1)
                    {
                        if (iMid - iLow > 1)
                        {
                            //---------------------------
                            //    Distribute mass to both sides
                            //---------------------------
                            //    The mass to each index (center index counts twice and receives double mass).
                            dy_this = this->dy[iMid] / (double)(iHigh - iLow);
                            //    Mass to iMid (double mass)
                            this->dy[iMid] = dy_this * 2.0;
                            //    Mass below iMid
                            for (i = iLow + 1; i < iMid; i++)
                                this->dy[i] += dy_this;
                            //    Mass above iMid
                            for (i = iMid + 1; i < iHigh; i++)
                                this->dy[i] += dy_this;
                            //---------------------------
                        }
                        else
                        {
                            //---------------------------
                            //    Distribute mass to high side
                            //---------------------------
                            //    The mass to each index (center index will get single mass).
                            dy_this = this->dy[iMid] / (double)(iHigh - iMid);
                            this->dy[iMid] = dy_this;
                            //    Mass from iMid to iHigh.
                            for (i = iMid + 1; i < iHigh; i++)
                                this->dy[i] += dy_this;
                            //---------------------------
                        }

                    }
                    else if (iMid - iLow > 1)
                    {
                        //---------------------------
                        //    Distribute mass to low side
                        //---------------------------
                        //    The mass to each index (center index will get single mass).
                        dy_this = this->dy[iMid] / (iMid - iLow);
                        //    Mass to iMid
                        this->dy[iMid] = dy_this;
                        //    Mass from iLow to iMid (must be added to mass that may already exist there).
                        for (i = iLow + 1; i < iMid; i++)
                            this->dy[i] += dy_this;
                        //---------------------------
                    }

                    //    Advance indexes for next trip
                    iLow = iMid; // The "subseqent flat region" on the next trip is set to the "preceding flat region" from the last trip.
                    iMid = iHigh; // We center on last trip's high index.
                }
                //----------------------------------------------------------------------------------

                //----------------------------------------------------------------------------------
                //    Integrate derivative. Ensure the mean of the original input.
                //----------------------------------------------------------------------------------
                sum_small = ymean - output[0];
                sum_big = 0.0;
                for (i = 0; i < nm1;)
                {
                    output[i + 1] = output[i] + this->dy[i];
                    if (output[++i] > ymean)
                        sum_big += output[i] - ymean;
                    else
                        sum_small += ymean - output[i];
                }

                //    All values will be shifted by a factor of 0.5*(sBig-sSmall)/(sBig+sSmall)
                sTarget = (sum_big + sum_small) / 2.0;
                sum_small = sTarget / sum_small;
                sum_big = sTarget / sum_big;
                //    Ensure that the range is not violated
                a_carry = 1.0;
                if (sum_small > 1.0)
                    a_carry = (ymean - ymin) / (ymean - output[0]) / sum_small;
                else if (sum_big > 1.0)
                    a_carry = (ymax - ymean) / (output[nm1] - ymean) / sum_big;
                if (a_carry < 1.0)
                {
                    sum_small *= a_carry;
                    sum_big *= a_carry;
                }

                //    Adjust elements on big and small sides so that mean is preserved.
                for (i = 0; i < n; i++)
                    if (output[i] < ymean)
                        output[i] = (output[i] - ymean) * sum_small + ymean;
                    else if (output[i] > ymean)
                        output[i] = (output[i] - ymean) * sum_big + ymean;
                //----------------------------------------------------------------------------------
            }
        }

        //    If the output is NaN (which is theoretically possible, but generally should not happen), change it all to be flat (yMean)
        if (std::isnan(output[0]))
        {
            for (i = 0; i < n; i++)
                output[i] = ymean;
        }
        else if (type == monotonic_regression_type::blended)
        {
            //    Blend the non-decreasing and monotonic functions according to blendFactor (a 1.0-->0.0 monotonic sigmoid transform of the length of the longest flat portion).
            a_carry = 1.0 - proportionBlendIncreasing;
            for (i = 0; i < n; i++)
            {
                output[i] = proportionBlendIncreasing * output[i] + a_carry * ynd[i];
            }
        }

        return num_trips_through_loop;
    }
}}}
