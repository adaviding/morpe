namespace Morpe.Numerics.I1
{
    public static class Util
    {
        /// <summary>
        /// Fills the vector with the series {0, 1, ..., toFill.Length-2, toFill.Length-1}.
        /// </summary>
        /// <param name="toFill">The vector to be filled.</param>
        public static void FillSeries(int[] toFill)
        {
            if (toFill == null)
                return;
            for (int i = 0; i < toFill.Length; i++)
                toFill[i] = i;
        }
        
        /// <summary>
        /// Generates a random ordering of elements.  Useful for shuffling.
        /// </summary>
        /// <param name="n">The length of the random ordering.</param>
        /// <returns>An array of length n containing a random ordering of the set of integers {0, ..., n-1}.</returns>
        public static int[] RandomOrder(int n)
        {
            int temp;
            int[] output = new int[n];
            FillSeries(output);
            for (int i = 0; i < n; i++)
            {
                int j = Morpe.Util.Rand.Next(n);
                if (j != i)
                {
                    temp = output[i];
                    output[i] = output[j];
                    output[j] = temp;
                }
            }
            return output;
        }
    }
}