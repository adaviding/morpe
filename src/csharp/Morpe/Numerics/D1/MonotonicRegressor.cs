using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading;
using Morpe.Validation;

namespace Morpe.Numerics.D1
{
    /// <summary>
    /// An object that performs monotonic regression.
    ///
    /// This just exposes a single <see cref="Run"/> method, but we use a class to recycle some heap resources (for efficiency) in case the caller wants to
    /// repeatedly call <see cref="Run"/>.
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
                if (this.derivative == null || this.derivative.Length < input.Length)
                {
                    this.derivative = new double[input.Length];
                    this.nonDecreasingValues = new double[input.Length];
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
                    this.derivative[i - 1] = dyThis;
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
                        dy: this.derivative);

                    // This will ensure that 'output' is non-decreasing, and some other nice properties.
                    EliminateDecreasingEnergy(
                        y: output,
                        dy: this.derivative,
                        yMean: yMean,
                        yMin: yMin,
                        yMax: yMax);
                }

                // Throw if cancellation requested
                cancellationToken.ThrowIfCancellationRequested();

                if (type == MonotonicRegressionType.Blended)
                {
                    // Save the non-decreasing values.  We will need them later.
                    Array.Copy(
                        sourceArray: output,
                        destinationArray: this.nonDecreasingValues,
                        length: output.Length);
                }

                //    Change a non-decreasing function into a monotonic function.
                if (type == MonotonicRegressionType.Increasing || type == MonotonicRegressionType.Blended)
                {
                    // Update derivative and its minimum value.
                    dyMin = UpdateDerivativeAndFindMin(y: output, dy: this.derivative);

                    // Throw if cancellation requested
                    cancellationToken.ThrowIfCancellationRequested();

                    // Convert a non-decreasing function to an increasing function (if necessary).
                    ConvertNonDecreasingToIncreasing(
                        y: output,
                        dy: this.derivative,
                        yMean: yMean,
                        yMin: yMin,
                        yMax: yMax,
                        dyMin: dyMin);

                    // Blending?
                    if (type == MonotonicRegressionType.Blended)
                    {
                        // Blend the non-decreasing and increasing functions together.

                        // The proportion of weight given to the increasing function.
                        double pwIncreasing = CalculateBlendFactor(
                            nonDecreasingValues: this.nonDecreasingValues,
                            length: output.Length);

                        Chk.True(pwIncreasing >= 0.0 && pwIncreasing <= 1.0, "{0} = {1}, expected to be in the range [0, 1].", nameof(pwIncreasing), pwIncreasing);

                        // The proportion of weight given to the non-decreasing function.
                        double pwNonDecreasing = 1.0 - pwIncreasing;

                        // Blend
                        for (i = 0; i < n; i++)
                        {
                            output[i] = pwIncreasing * output[i] + pwNonDecreasing * this.nonDecreasingValues[i];
                        }
                    }
                }

                // Defend against a degenerate case.
                if (Double.IsNaN(output[0]))
                {
                    // In the degenerate case where the output should be flat, we might have some NaN and really it should just be flat.
                    for (i = 0; i < n; i++)
                        output[i] = yMean;
                }
            }

            return numTripsThroughLoop;
        }

        /// <summary>
        /// Intermediate storage for the derivative of 'y'.  This changes as 'y' changes.
        /// </summary>
        private double[] derivative;

        /// <summary>
        /// A mutex ensures that the <see cref="Run"/> method cannot run more than once.
        /// </summary>
        private readonly object mutex = new object();

        /// <summary>
        /// Intermediate storage for an intermediate result:  The non-decreasing function.
        /// </summary>
        private double[] nonDecreasingValues;

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

        /// <summary>
        /// This is needed in the context of calculating a <see cref="MonotonicRegressionType.Blended"/> monotonic regression,
        /// an operation basically blends the outputs of the other operation types.
        /// * <see cref="MonotonicRegressionType.NonDecreasing"/>
        /// * <see cref="MonotonicRegressionType.Increasing"/>
        ///
        /// This method calculates a "blend factor", or basically a weighting used to mix the two other output types.
        /// </summary>
        /// <param name="nonDecreasingValues">The monotonic regression output for <see cref="MonotonicRegressionType.NonDecreasing"/>.</param>
        /// <param name="length">The actual number of non-decreasing values, since the vector passed in may have extra padding.</param>
        /// <returns>The proportion of weight to the <see cref="MonotonicRegressionType.Increasing"/> values.</returns>
        private static double CalculateBlendFactor(
            [NotNull] double[] nonDecreasingValues,
            int length)
        {
            Chk.LessOrEqual(length, nonDecreasingValues.Length, "The length is {0}, but the vector length is {1}.", length, nonDecreasingValues.Length);

            // Prepare for blending.
            // The blend factor is a number in the range (0, 1) representing the proportion of weight given to the "increasing" function.
            double output = 1.0;

            // We prepare for blending in 2 ways:
            // 1.  Compute the blend factor (proportionBlendIncreasing)
            // 2.  Save the non-decreasing values so that later they be blended with the increasing values.

            // Measure the size of the largest non-decreasing span.
            int runLen = 0;
            int maxRunLen = 0;

            for (int i = 1; i < length; i++)
            {
                // Measure the size of the largest non-decreasing span.
                if (nonDecreasingValues[i] != nonDecreasingValues[i - 1])
                    runLen = 0;
                else
                    runLen++;
                if (runLen > maxRunLen)
                    maxRunLen = runLen;
            }

            // Compute the blend factor
            output = 1.0 / (1.0 + Math.Pow((double)(4 * maxRunLen) / (double)(length - 1), 3.0));
            return output;
        }

        /// <summary>
        /// Modifies the non-decreasing values of 'y' to be strictly increasing, so that each value of 'y' is greater than (not equal to) the prior value.
        /// </summary>
        /// <param name="y">The tabulated values of the function 'y'.  On input, we know that these values are non-decreasing, but not necessarily
        /// increasing.  On output, these values will be modified and will contain the increasing values.</param>
        /// <param name="dy">The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].  These values will be modified.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        ///
        ///     On exit, these values will be meaningless:  They are not kept in sync with 'y'.</param>
        /// <param name="yMean">The intended mean of 'y'.  This value will be preserved.</param>
        /// <param name="yMin">The intended min of 'y'.  All values of 'y' will be greater than or equal to this value, but no guarantee that 'y' will actually
        /// have this min.</param>
        /// <param name="yMax">The intended max of 'y'.  All values of 'y' will be less than or equal to this value, but no guarantee that 'y' will actually
        /// have this max.</param>
        /// <param name="dyMin">The minimum value of 'dy' (pre-computed for us).</param>
        private static void ConvertNonDecreasingToIncreasing(
            double[] y,
            double[] dy,
            double yMean,
            double yMin,
            double yMax,
            double dyMin)
        {
            int n = y.Length;
            int nm1 = n - 1;

            double zeroFloor = (y[nm1] - y[0]) * 0.00001 / n;

            //    Only continue if the minimum derivative is zero (or tiny) and if the function is not totally flat.
            if (y[nm1] > y[0] && dyMin <= zeroFloor)
            {
                //----------------------------------------------------------------------------------
                // First we alter the derivative to be positive everywhere.
                //----------------------------------------------------------------------------------

                // Init indexes at a location where [iLow, iMid] straddles the first flat region.
                int iLow, iMid = 0, iHigh;
                //    iMid:  The index of a non-flat spot
                //    iLow:  The preceding non-flat index before iMid
                //    iHigh: The subsequent non-flat index following iMid
                iLow = -1;
                while (iMid < nm1 && dy[iMid] <= zeroFloor) iMid++;


                // Advance from iMid's initial value up through iMid = nm-1.
                while (iMid < nm1)
                {
                    // Advance index of iHigh to the subsequent non-flat index following iMid
                    iHigh = iMid + 1;
                    while (iHigh < nm1 && dy[iHigh] <= zeroFloor) iHigh++;

                    if (iHigh - iMid > 1)
                    {
                        if (iMid - iLow > 1)
                        {
                            //---------------------------
                            //    Distribute mass to both sides
                            //---------------------------
                            //    The mass to each index (center index counts twice and receives double mass).
                            double dyThis = dy[iMid] / (iHigh - iLow);
                            //    Mass to iMid (double mass)
                            dy[iMid] = dyThis * 2.0;
                            //    Mass below iMid
                            for (int i = iLow + 1; i < iMid; i++)
                                dy[i] += dyThis;
                            //    Mass above iMid
                            for (int i = iMid + 1; i < iHigh; i++)
                                dy[i] += dyThis;
                            //---------------------------
                        }
                        else
                        {
                            //---------------------------
                            //    Distribute mass to high side
                            //---------------------------
                            //    The mass to each index (center index will get single mass).
                            double dyThis = dy[iMid] / (iHigh - iMid);
                            dy[iMid] = dyThis;
                            //    Mass from iMid to iHigh.
                            for (int i = iMid + 1; i < iHigh; i++)
                                dy[i] += dyThis;
                            //---------------------------
                        }

                    }
                    else if (iMid - iLow > 1)
                    {
                        //---------------------------
                        //    Distribute mass to low side
                        //---------------------------
                        //    The mass to each index (center index will get single mass).
                        double dyThis = dy[iMid] / (iMid - iLow);
                        //    Mass to iMid
                        dy[iMid] = dyThis;
                        //    Mass from iLow to iMid (must be added to mass that may already exist there).
                        for (int i = iLow + 1; i < iMid; i++)
                            dy[i] += dyThis;
                        //---------------------------
                    }

                    //    Advance indexes for next trip
                    iLow = iMid; // The "subseqent flat region" on the next trip is set to the "preceding flat region" from the last trip.
                    iMid = iHigh; // We center on last trip's high index.
                }

                //----------------------------------------------------------------------------------
                // Next we integrate the derivative while ensuring the mean of the original input.
                //----------------------------------------------------------------------------------
                double energyBelow = yMean - y[0];
                double energyAbove = 0.0;
                for (int i = 0; i < nm1;)
                {
                    y[i + 1] = y[i] + dy[i];
                    if (y[++i] > yMean)
                        energyAbove += y[i] - yMean;
                    else
                        energyBelow += yMean - y[i];
                }

                // Prepare to scale the values above and below differently.
                double middleEnergy = (energyAbove + energyBelow) / 2.0;
                double scalarBelow = middleEnergy / energyBelow;
                double scalarAbove = middleEnergy / energyAbove;

                // Limit the size of the scalars to ensure that the span of y-values is correct.
                double carry = 1.0;
                if (scalarBelow > 1.0)
                    carry = (yMean - yMin) / (yMean - y[0]) / scalarBelow;
                else if (scalarAbove > 1.0)
                    carry = (yMax - yMean) / (y[nm1] - yMean) / scalarAbove;
                if (carry < 1.0)
                {
                    scalarBelow *= carry;
                    scalarAbove *= carry;
                }

                // Adjust elements on big and small sides so that mean is preserved.
                for (int i = 0; i < n; i++)
                    if (y[i] < yMean)
                        y[i] = (y[i] - yMean) * scalarBelow + yMean;
                    else if (y[i] > yMean)
                        y[i] = (y[i] - yMean) * scalarAbove + yMean;
            }
        }

        /// <summary>
        /// This is executed after <see cref="AnnihilateDecreasingEnergy"/> to ensure that all decreasing energy has been finally eliminated.  This method is
        /// not the best way to annihilate the decreasing energy, but it does guarantee the following:
        /// * 'y' will be non-decreasing.
        /// * 'dy' will be non-negative.
        /// * The mean of 'y' will be 'yMean'.
        /// * The min of 'y' shall not be below 'yMin'.
        /// * The max of 'y' shall not be above 'yMax'.
        ///
        /// WARNING:  When this function exits, the values of 'dy' will be invalid.
        /// </summary>
        /// <param name="y">The tabulated values of the function 'y'.  These values will be modified.</param>
        /// <param name="dy">The diff of 'y' values.  For all y:  dy[i] = y[i+1] - y[i].  These values are NOT modified.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.
        ///
        ///     When this function returns, the values of this variable will be invalid.  The values of 'y' have changed
        ///     but the values of 'dy' have not.</param>
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
            double carry = 0.0;
            for (int i = 0; i < nm1;)
            {
                if (dy[i++] < 0.0)
                    carry -= dy[i - 1];

                // Set 'y'
                y[i] += carry;

                // Watch out for tiny errors with floating point precision.  Ensure it is really non-decreasing.
                if (y[i] < y[i - 1])
                {
                    carry += y[i - 1] - y[i];
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
            double scalar = (rngY - carry) / rngY;

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

        /// <summary>
        /// Updates the values of 'dy' given the current values of 'y' and also finds the minimum of 'dy'.
        /// </summary>
        /// <param name="y">The current values of 'y'.</param>
        /// <param name="dy">A container which holds the output values of 'dy'.
        ///     This vector may have extra length allocated (just for efficient heap utilization).  The extra elements
        ///     can be ignored.</param>
        /// <returns>The minimum value of 'dy'.</returns>
        private static double UpdateDerivativeAndFindMin(
            double[] y,
            double[] dy)
        {
            double output = double.MaxValue;

            for (int i = 1; i < y.Length; i++)
            {
                double diff = y[i] - y[i - 1];
                if (diff < output)
                {
                    output = diff;
                }

                dy[i - 1] = diff;
            }

            return output;
        }
    }
}
