namespace Morpe.Numerics.D1
{
    public static class Util
    {
        /// <summary>
        /// Performs a binary search on a tabulated function.
        /// </summary>
        /// <param name="xTabulated">The tabulated function.  The values must be non-decreasing.</param>
        /// <param name="xTarget">The target value.</param>
        /// <returns>Returns the zero-based index of the target value within the range of the tabulated function.
        /// It can be a continuous (non-integer) value, in the range [0.0, xTabulated.Length], representing an exact
        /// placement between two tabulated values (as a linear-interpolant between those values).  It can also be
        /// -1.0 if below the tabulated range, or xTabulated.Length if above the range.</returns>
        public static double BinarySearchOfNonDecreasing(double[] xTabulated, double xTarget)
        {
            int iMin = 0;
            if (xTarget < xTabulated[iMin]) return -1.0;

            int iMax = xTabulated.Length - 1;
            if (xTarget > xTabulated[iMax]) return (double)xTabulated.Length;

            int iHigh = iMax;
            int iLow = iMin;
            int iMid;
            while (true)
            {
                iMid = (iHigh + iLow) / 2;

                //-----------------------------------------------------------------------------
                //    Finish
                //-----------------------------------------------------------------------------
                if (xTabulated[iMid] == xTarget)
                {
                    //    Put both indexes in range and rely on execution of next "if" statement to finish.
                    iLow = iHigh = iMid;
                }
                if (iHigh - iLow < 5)
                {
                    //    Get iLow to be the lowest index of xTabulated possibly equal to the xTarget.
                    //    Otherwise it should be just below the xTarget.
                    while (iLow > iMin && xTabulated[iLow - 1] >= xTarget) iLow--;
                    while (iLow < iMax && xTabulated[iLow + 1] < xTarget) iLow++;
                    if (iLow < iMax && xTabulated[iLow + 1] == xTarget) iLow++;
                    //    Get iHigh to be the highest index of xTabulated possibly equal to the xTarget.
                    //    Otherwise it should be just above the xTarget.
                    while (iHigh < iMax && xTabulated[iHigh + 1] <= xTarget) iHigh++;
                    while (iHigh > iMin && xTabulated[iHigh - 1] > xTarget) iHigh--;
                    if (iHigh > iMin && xTabulated[iHigh - 1] == xTarget) iHigh--;

                    if (iHigh == iLow) // One unique xTarget value.  Return its index.
                        return (double)iLow;
                    if (xTabulated[iLow] == xTarget) // Many unique xTarget values exist.  Return the center of the range.
                        return 0.5 * (double)(iLow + iHigh);
                    if (iHigh - iLow != 1)
                        //    Return linear interpolant.
                        return (double)iLow + (xTarget - xTabulated[iLow]) / (xTabulated[iHigh] - xTabulated[iLow]);
                }

                //-----------------------------------------------------------------------------
                //    Continue
                //-----------------------------------------------------------------------------
                //  Reduce possible range by 1/2
                if (xTabulated[iMid] > xTarget)
                    iLow = iMid;
                else // Must be < xTarget
                    iHigh = iMid;
            }
        }
        
        /// <summary>
        /// Performs linear interpolation of a tabulated function y -> f(x).
        /// </summary>
        /// <param name="xTabulated">Non-decreasing tabulated values of an independent variable "x".</param>
        /// <param name="yTabulated">Tabulated values of a dependent variable "y".</param>
        /// <param name="xTarget">The value of x for which the matching value of y is linearly interpolated.</param>
        /// <returns>The interpolated value of the function y -> f(x) where x=xTarget.</returns>
        public static double Linterp(double[] xTabulated, double[] yTabulated, double xTarget)
        {
            double i = BinarySearchOfNonDecreasing(xTabulated, xTarget);
            if (i < 0.0 || i >= (double)(xTabulated.Length - 1))
                return double.NaN;
            double iMod = i - (double)(int)i;
            if (iMod < 1e-10 || iMod > 0.9999999999)
                return yTabulated[(int)(i + 0.5)];
            else
                return (1.0 - iMod) * yTabulated[(int)i] + iMod * yTabulated[1 + (int)i];
        }
        
        /// <summary>
        /// Quicksorts an indexed list by changing the order of the indices and not the list itself.
        /// </summary>
        /// <param name="idx">Zero-based indices into the list x.  On output, these indices are reordered so that x[idx] is sorted.</param>
        /// <param name="x">The list to be sorted.  The elements are not changed on output.</param>
        /// <param name="left">A zero-based index identifying the lowest element of idx where sorting is conducted.</param>
        /// <param name="right">A zero-based index identifying the highest element of idx where sorting is conducted.</param>
        public static void QuickSortIndex(int[] idx, double[] x, int left, int right)
        {
            int i = left, j = right;
            int tmp;
            double pivot = x[idx[(left + right) / 2]];

            //    Partition
            while (i <= j)
            {
                while (x[idx[i]] < pivot) i++;
                while (x[idx[j]] > pivot) j--;
                if (i <= j)
                {
                    tmp = idx[i];
                    idx[i++] = idx[j];
                    idx[j--] = tmp;
                }
            }

            //    Recursion
            if (left < j)
                QuickSortIndex(idx, x, left, j);
            if (i < right)
                QuickSortIndex(idx, x, i, right);
        }
    }
}