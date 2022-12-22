using System;

namespace Morpe.Numerics.F1
{
    public static class Util
    {
        /// <summary>
        /// Gets a subset of the input array.
        /// </summary>
        /// <param name="input">The input array.</param>
        /// <param name="iFirst">The 0-based index of the first element to include in the output.</param>
        /// <param name="iLast">The 0-based inde of the last element to include in the output.</param>
        /// <returns>A newly created array with elements equal to a subset of the input array.</returns>
        public static float[] GetSubarray(float[] input, int iFirst, int iLast)
        {
            int len = iLast - iFirst + 1;
            if (len < 0 || iFirst<0 || iLast >= input.Length)
                throw new ArgumentException("The indices specified do not properly define a section of the input array.");

            float[] output = new float[len];
            for (int i = 0; i < len; i++)
                output[i] = input[i + iFirst];
            return output;
        }
    }
}
