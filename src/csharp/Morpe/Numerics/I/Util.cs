using System;
using System.Diagnostics.CodeAnalysis;
using Morpe.Validation;

namespace Morpe.Numerics.I
{
    public static class Util
    {
        /// <summary>
        /// Fill the vector 'x' with a series of numbers like 0, 1, 2, ..., x.Length-1.
        /// </summary>
        /// <param name="x">The vector to be filled.</param>
        public static void FillSeries(
            [NotNull] int[] x)
        {
            Chk.NotNull(x, nameof(x));

            for (int i = 0; i < x.Length; i++)
            {
                x[i] = i;
            }
        }
        
        /// <summary>
        /// Allocates a jagged 2D array.
        /// </summary>
        /// <param name="height">Length of first index.</param>
        /// <param name="width">Length of second index.</param>
        /// <returns>The jagged 2D array.</returns>
        [return: NotNull]
        public static int[][] JaggedSquare(int height, int width)
        {
            int[][] output = new int[height][];

            for (int i = 0; i < height; i++)
            {
                output[i] = new int[width];
            }

            return output;
        }
        
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
