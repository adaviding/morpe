using System;

namespace Morpe.Numerics.I
{
    public static class Util
    {
        /// <summary>
        /// Returns the entry (a,b) of the Pascal matrix.  The parameters are interchangeable because the matrix is symmetric.
        /// The value returned is equal to:  (a+b)! / a! / b!
        /// </summary>
        /// <param name="a">The row of the Pascal matrix.</param>
        /// <param name="b">The column of the Pascal matrix.</param>
        public static int Pascal(int a, int b)
        {
            int aa = Math.Max(a, b);
            int bb = Math.Min(a, b);
            int c = a + b;
            int output = 1;
            int i;
            for (i = aa + 1; i <= c; i++)
                output *= i;
            for (i = 2; i <= bb; i++)
                output /= i;
            return output;
        }
    }
}
