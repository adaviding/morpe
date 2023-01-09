using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading;

namespace Morpe.Numerics.D1
{
    /// <summary>
    /// An object that performs monotonic regression.
    /// </summary>
    public class MonotonicRegressor
    {
        /// <summary>
        /// Performs a monotonic regression of a tabulated function.  This method calculates a monotonic function that is approximately equal
        /// to the tabulated function provided.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="type">The type of monotonic regression.</param>
        /// <param name="input">The tabulated function.</param>
        /// <param name="output">The monotonic function.  This should be supplied as vector of length identical to input.
        /// When the function has finished executing, the following are guaranteed to be true:
        ///        mean(output) == mean(input)
        ///        min(output) &gt;= min(input)
        ///        max(output) &lt;= max(input)
        /// </param>
        /// <returns>The number of repetitions through a loop.</returns>
        public int Run(
            CancellationToken cancellationToken,
            MonotonicRegressionType type,
            [NotNull] double[] input,
            [NotNull] double[] output)
        {
            if (input == null || output == null || output.Length < input.Length)
                throw new ArgumentException("The arguments cannot be null, and the length of the output must be at least the length of the input.");

            // Count the number of trips through the loop.
            int numTripsThroughLoop = 0;

            // Prevent concurrent execution.  We are using member variables to store intermediate answers,
            // and concurrency would create catastrophic interference.
            lock (this.mutex)
            {
                if (this.dy == null || this.dy.Length < input.Length)
                {
                    this.dy = new double[input.Length];
                    this.ynd = new double[input.Length];
                }

                int limitNumTripsThroughLoop = 100 + (int)(30.0 * Math.Log((double)input.Length));
                int n = input.Length;
                int nm1 = input.Length - 1;
                int i, j;
                //    Compute min, max, mean, etc.
                double yMin = input[0];
                double yMax = input[0];
                double yMean = input[0];
                double dyThis;
                double dyMin = input[1] - input[0];
                output[0] = input[0];
                for (i = 1; i < n; i++)
                {
                    if (input[i] < yMin) yMin = input[i];
                    if (input[i] > yMax) yMax = input[i];
                    yMean += input[i];
                    dyThis = input[i] - input[i - 1];
                    this.dy[i - 1] = dyThis;
                    if (dyThis < dyMin)
                        dyMin = dyThis;
                    output[i] = input[i];
                }

                yMean /= n;
                double rngY = yMax - yMin;
                double dyMinThreshold = -rngY * 0.001;

                if (dyMin >= 0.0)
                {
                    if (dyMin > 0.0
                        || dyMin == 0.0 && type == MonotonicRegressionType.NonDecreasing)
                    {
                        // Let's say we took 0 trips through a loop.
                        return 0;
                    }
                }
                else
                {
                    // Exhaustively annihilate decreasing energy in an unbiased way.  This is the best approach to achieving
                    // a non-decreasing function, but there is no guarantee that it will work.  Even if it doesn't work, it
                    // typically gets us most of the way home.
                    numTripsThroughLoop = AnnihilateDecreasingEnergy(
                        cancellationToken,
                        ref dyMin,
                        dyMinThreshold,
                        limitNumTripsThroughLoop,
                        y: output,
                        dy: this.dy);

                    // This will ensure that 'output' is non-decreasing, and some other nice properties.
                    EliminateDecreasingEnergy(
                        y: output,
                        dy: this.dy,
                        yMean: yMean,
                        yMin: yMin,
                        yMax: yMax);
                }

                // Prepare for blending.
                // The blend factor is a number in the range (0, 1) representing the proportion of weight given to the "increasing" function.
                double proportionBlendIncreasing = 1.0;
                if (type == MonotonicRegressionType.Blended)
                {
                    // We prepare for blending in 2 ways:
                    // 1.  Compute the blend factor (proportionBlendIncreasing)
                    // 2.  Save the non-decreasing values so that later they be blended with the increasing values.

                    // Measure the size of the largest non-decreasing span.
                    int runLen = 0;
                    int maxRunLen = 0;

                    this.ynd[0] = output[0];  // Save the non-decreasing values
                    for (i = 1; i < n; i++)
                    {
                        this.ynd[i] = output[i];  // Save the non-decreasing values

                        // Measure the size of the largest non-decreasing span.
                        if (output[i] != output[i - 1])
                            runLen = 0;
                        else
                            runLen++;
                        if (runLen > maxRunLen)
                            maxRunLen = runLen;
                    }

                    // Compute the blend factor
                    proportionBlendIncreasing = 1.0 / (1.0 + Math.Pow((double)(4 * maxRunLen) / (double)(n - 1), 3.0));
                }

                //    Change a non-decreasing function into a monotonic function.
                if (type == MonotonicRegressionType.Increasing || type == MonotonicRegressionType.Blended)
                {
                    //    Compute min, max, and recompute derivative.
                    dyMin = output[1] - output[0];
                    for (i = 1; i < n; i++)
                    {
                        dyThis = output[i] - output[i - 1];
                        this.dy[i - 1] = dyThis;
                        if (dyThis < dyMin)
                            dyMin = dyThis;
                    }

                    sumBelow = (output[nm1] - output[0]) * 0.00001 / (double)n;

                    //    Only continue if the minimum derivative is zero (or tiny) and if the function is not totally flat.
                    if (output[nm1] > output[0] && dyMin <= sumBelow)
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
                        while (iMid < nm1 && this.dy[iMid] <= sumBelow) iMid++;
                        //---------------------------

                        //    Advance from iMid's initial value up through iMid = nm-1.
                        while (iMid < nm1)
                        {
                            //    Advance index of iHigh to the subsequent non-flat index following iMid
                            iHigh = iMid + 1;
                            while (iHigh < nm1 && this.dy[iHigh] <= sumBelow) iHigh++;

                            if (iHigh - iMid > 1)
                            {
                                if (iMid - iLow > 1)
                                {
                                    //---------------------------
                                    //    Distribute mass to both sides
                                    //---------------------------
                                    //    The mass to each index (center index counts twice and receives double mass).
                                    dyThis = this.dy[iMid] / (double)(iHigh - iLow);
                                    //    Mass to iMid (double mass)
                                    this.dy[iMid] = dyThis * 2.0;
                                    //    Mass below iMid
                                    for (i = iLow + 1; i < iMid; i++)
                                        this.dy[i] += dyThis;
                                    //    Mass above iMid
                                    for (i = iMid + 1; i < iHigh; i++)
                                        this.dy[i] += dyThis;
                                    //---------------------------
                                }
                                else
                                {
                                    //---------------------------
                                    //    Distribute mass to high side
                                    //---------------------------
                                    //    The mass to each index (center index will get single mass).
                                    dyThis = this.dy[iMid] / (double)(iHigh - iMid);
                                    this.dy[iMid] = dyThis;
                                    //    Mass from iMid to iHigh.
                                    for (i = iMid + 1; i < iHigh; i++)
                                        this.dy[i] += dyThis;
                                    //---------------------------
                                }

                            }
                            else if (iMid - iLow > 1)
                            {
                                //---------------------------
                                //    Distribute mass to low side
                                //---------------------------
                                //    The mass to each index (center index will get single mass).
                                dyThis = this.dy[iMid] / (iMid - iLow);
                                //    Mass to iMid
                                this.dy[iMid] = dyThis;
                                //    Mass from iLow to iMid (must be added to mass that may already exist there).
                                for (i = iLow + 1; i < iMid; i++)
                                    this.dy[i] += dyThis;
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
                        sumBelow = yMean - output[0];
                        sumAbove = 0.0;
                        for (i = 0; i < nm1;)
                        {
                            output[i + 1] = output[i] + this.dy[i];
                            if (output[++i] > yMean)
                                sumAbove += output[i] - yMean;
                            else
                                sumBelow += yMean - output[i];
                        }

                        //    All values will be shifted by a factor of 0.5*(sBig-sSmall)/(sBig+sSmall)
                        sTarget = (sumAbove + sumBelow) / 2.0;
                        sumBelow = sTarget / sumBelow;
                        sumAbove = sTarget / sumAbove;
                        //    Ensure that the range is not violated
                        aCarry = 1.0;
                        if (sumBelow > 1.0)
                            aCarry = (yMean - yMin) / (yMean - output[0]) / sumBelow;
                        else if (sumAbove > 1.0)
                            aCarry = (yMax - yMean) / (output[nm1] - yMean) / sumAbove;
                        if (aCarry < 1.0)
                        {
                            sumBelow *= aCarry;
                            sumAbove *= aCarry;
                        }

                        //    Adjust elements on big and small sides so that mean is preserved.
                        for (i = 0; i < n; i++)
                            if (output[i] < yMean)
                                output[i] = (output[i] - yMean) * sumBelow + yMean;
                            else if (output[i] > yMean)
                                output[i] = (output[i] - yMean) * sumAbove + yMean;
                        //----------------------------------------------------------------------------------
                    }
                }

                //    If the output is NaN (which is theoretically possible, but generally should not happen), change it all to be flat (yMean)
                if (Double.IsNaN(output[0]))
                {
                    for (i = 0; i < n; i++)
                        output[i] = yMean;
                }
                else if (type == MonotonicRegressionType.Blended)
                {
                    //    Blend the non-decreasing and monotonic functions according to blendFactor (a 1.0-->0.0 monotonic sigmoid transform of the length of the longest flat portion).
                    aCarry = 1.0 - proportionBlendIncreasing;
                    for (i = 0; i < n; i++)
                    {
                        output[i] = proportionBlendIncreasing * output[i] + aCarry * this.ynd[i];
                    }
                }
            }

            return numTripsThroughLoop;
        }

        /// <summary>
        /// Intermediate storage for the derivative and altered derivative.
        /// </summary>
        private double[] dy;

        /// <summary>
        /// A mutex ensures that the <see cref="Run"/> method cannot run more than once.
        /// </summary>
        private readonly object mutex = new object();

        /// <summary>
        /// Intermediate storage for an intermediate result (the non-decreasing function).
        /// </summary>
        private double[] ynd;

        /// <summary>
        /// Repeatedly spread decreasing energy to immediate neighbors until a threshold is met, or until we have gone through the
        /// loop too many times.  Decreasing energy is annihilated by increasing energy.  Thus, it tends to be annihilated by this
        /// spreading action.
        ///
        /// This is the best approach to annihilating the decreasing energy, but this method may not eliminate all of it.  Even so,
        /// it typically gets us most of the way home.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="dyMin">The minimum value of 'dy'.  This value will be updated with each trip through the loop.</param>
        /// <param name="dyMinThreshold">The loop may exit when 'dyMin' exceeds this value.  This should be a small negative number.</param>
        /// <param name="limitNumTripsThroughLoop">The maximum number of trips through the loop.</param>
        /// <param name="y">The tabulated values of the function 'y'.  These values will be modified.</param>
        /// <param name="dy">The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].  These values will be modified.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.</param>
        /// <returns>The number of trips taken through the loop.</returns>
        private static int AnnihilateDecreasingEnergy(
            CancellationToken cancellationToken,
            ref double dyMin,
            double dyMinThreshold,
            int limitNumTripsThroughLoop,
            double[] y,
            double[] dy)
        {
            // The number of trips through the loop.
            int output = 0;

            int nm1 = y.Length - 1;

            while (dyMin < dyMinThreshold)
            {
                int i = 0; // Point to beginning of dy
                int j = y.Length - 2; // Point to end of dy.
                //    Go through all of dy
                while (j > 0)
                {
                    // Check for cancellation.
                    cancellationToken.ThrowIfCancellationRequested();

                    RepelDecreasingEnergyFromIndex(
                        idx: i,
                        y: y,
                        dy: dy);

                    if (j != i)
                    {
                        RepelDecreasingEnergyFromIndex(
                            idx: j,
                            y: y,
                            dy: dy);
                    }

                    j--;
                    i++;
                }

                //    Update minimum of dy
                dyMin = dy[0];
                for (i = 1; i < nm1; i++)
                    if (dy[i] < dyMin)
                        dyMin = dy[i];

                //    Increment number of trips.
                output++;
                if (output >= limitNumTripsThroughLoop)
                    break;
            }

            return output;
        }

        private static void ConvertNonDecreasingToIncreasing()
        {

        }

        /// <summary>
        /// This is executed after <see cref="AnnihilateDecreasingEnergy"/> to ensure that all decreasing energy has been finally eliminated.  This method is
        /// not the best way to annihilate the decreasing energy, but it does guarantee the following:
        /// * 'y' will be non-decreasing.
        /// * 'dy' will be non-negative.
        /// * The mean of 'y' will be 'yMean'.
        /// * The min of 'y' shall not be below 'yMin'.
        /// * The max of 'y' shall not be above 'yMax'.
        /// </summary>
        /// <param name="y">The tabulated values of the function 'y'.  These values will be modified.</param>
        /// <param name="dy">The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].  These values will be modified.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.</param>
        /// <param name="yMean">The mean of 'y'.  Although we will change the 'y' values, we will ensure that it still has this mean.</param>
        /// <param name="yMin">The minimum allowed value of 'y'.  No value of 'y' will be below this minimum.</param>
        /// <param name="yMax">The maximum allowed value of 'y'.  No value of 'y' will be above this maximum.</param>
        private static void EliminateDecreasingEnergy(
            double[] y,
            double[] dy,
            double yMean,
            double yMin,
            double yMax)
        {
            int n = y.Length;
            int nm1 = n - 1;

            // Fix small dy values to be non-negative (so that 'y' is non-decreasing).
            // This operation will bias 'y' in an undesirable way:  The function will increase too rapidly, so there is a "clean-up" procedure after.
            double aCarry = 0.0;
            for (int i = 0; i < nm1;)
            {
                if (dy[i++] < 0.0)
                    aCarry -= dy[i - 1];

                // Set 'y'
                y[i] += aCarry;

                // Watch out for tiny errors with floating point precision.  Ensure it is really non-decreasing.
                if (y[i] < y[i - 1])
                {
                    aCarry += y[i - 1] - y[i];
                    y[i] = y[i - 1];
                }
            }

            // The sum of deviations (from the mean) for values above the mean.
            double sumAboveMean = 0.0;

            // The sum of deviations (from the mean) for values below the mean.
            double sumBelowMean = 0.0;

            // The above operation may have (1) expanded the range of our function and (2) increased the mean of the function.
            // Let's measure this so we can fix it.
            double y0 = y[0];
            double rngY = y[nm1] - y0;

            // We can multiply each y-value with this scalar to eliminate the range expansion.
            double scalar = (rngY - aCarry) / rngY;

            for (int i = 0; i < n; i++)
            {
                // Counter the expansion of our function's range
                double yi = y[i] = (y[i] - y0) * scalar + y0;

                // Now measure what happened to the mean of the function
                if (yi < yMean)
                    sumBelowMean += yMean - yi;
                else
                    sumAboveMean += yi - yMean;
            }

            // See if we changed the mean.
            if (sumAboveMean != sumBelowMean)
            {
                // Yes it changed, we must fix it.
                if (sumBelowMean == 0.0 || sumAboveMean == 0.0)
                {
                    // This is a degenerate situation, like if the function is flat.
                    if (sumBelowMean > 0.0)
                    {
                        double adj = sumBelowMean / (2 * n); // The extent to which the actual mean is below the desired mean.
                        for (int i = 0; i < n; i++)
                            y[i] += adj;
                    }
                    else
                    {
                        double adj = sumAboveMean / (2 * n); // The extent to which the actual mean is above the desired mean.
                        for (int i = 0; i < n; i++)
                            y[i] -= adj;
                    }
                }
                else
                {
                    // The first constraint:  The distance of all values from the mean will be scaled so that sumBelowMean==sumAboveMean
                    // (i.e. so that the mean unaltered).
                    double target = (sumAboveMean + sumBelowMean) / 2.0;
                    double adjBelow = target / sumBelowMean;
                    double adjAbove = target / sumAboveMean;

                    // The second constraint:  Ensure that the range does not exceed the original range.
                    double adj = 1.0;
                    if (adjBelow > 1.0)
                        adj = (yMean - yMin) / (yMean - y[0]) / adjBelow;
                    else if (adjAbove > 1.0)
                        adj = (yMax - yMean) / (y[nm1] - yMean) / adjAbove;
                    if (adj < 1.0)
                    {
                        adjBelow *= adj;
                        adjAbove *= adj;
                    }

                    // Now apply both constraints simultaneously.
                    for (int i = 0; i < n; i++)
                    {
                        double yi = y[i];

                        if (yi < yMean)
                            y[i] = (yi - yMean) * adjBelow + yMean;
                        else if (yi > yMean)
                            y[i] = (yi - yMean) * adjAbove + yMean;
                    }
                }
            }
        }

        /// <summary>
        /// If decreasing energy is present at the given index, then it is moved to the adjacent indices (idx-1, idx+1).
        /// If either of the adjacent indices has increasing energy, it will annihilate the decreasing energy in equal
        /// measure.
        /// </summary>
        /// <param name="idx">The zero-based index into 'y' where the adjustment is to be made.</param>
        /// <param name="y">The tabulated values of the function 'y'.  These values will be modified.</param>
        /// <param name="dy">The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].  These values will be modified.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.</param>
        private static void RepelDecreasingEnergyFromIndex(
            int idx,
            double[] y,
            double[] dy)
        {
            double halfDy = dy[idx];
            if (halfDy < 0.0)
            {
                // Make it half.
                halfDy /= 2.0;

                //    Element y(i+1) must be incremented and y(i) must be decremented.
                y[idx + 1] -= halfDy; // increment y(i+1)
                y[idx] += halfDy; // decrement y(i)

                //    Now, 'dy' must be adjusted to account for this new change.
                if (idx > 0)
                    dy[idx - 1] += halfDy; // The difference between y(i-1) and y(i) goes down.
                if (idx < y.Length - 2)
                    dy[idx + 1] += halfDy; // The difference between y(i+1) and y(i+2) goes down.
                dy[idx] = 0.0; // The difference between y(i) and y(i+1) goes to zero.
            }
        }
    }
}
