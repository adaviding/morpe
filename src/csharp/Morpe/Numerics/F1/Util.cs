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
        
        /// <summary>
		/// Quicksorts an indexed list by changing the order of the indices and not the list itself.
		/// </summary>
		/// <param name="idx">Zero-based indices into the list x.  On output, these indices are reordered so that x[idx] is sorted.</param>
		/// <param name="x">The list to be sorted.  The elements are not changed on output.</param>
		/// <param name="left">A zero-based index identifying the lowest element of idx where sorting is conducted.</param>
		/// <param name="right">A zero-based index identifying the highest element of idx where sorting is conducted.</param>
		public static void QuickSortIndex(int[] idx, float[] x, int left, int right)
		{
			int i = left, j = right;
			int tmp;
			double pivot = x[idx[(left + right) / 2]];

			//	Partition
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

			//	Recursion
			if (left < j)
				QuickSortIndex(idx, x, left, j);
			if (i < right)
				QuickSortIndex(idx, x, i, right);
		}
		
		/// <summary>
		/// Quicksorts an indexed table by an individual column by changing the order of the indices and not the table itself.
		/// </summary>
		/// <param name="idx">Zero-based indices into the table x.  On output, these indices are reordered so that x[idx][iCol] is sorted.</param>
		/// <param name="x">The table to be sorted.  The elements are not changed on output.</param>
		/// <param name="iCol">The zero-based column index which determines the rank order of rows within the table.</param>
		/// <param name="left">A zero-based index identifying the lowest row of idx where sorting is conducted.</param>
		/// <param name="right">A zero-based index identifying the highest row of idx where sorting is conducted.</param>
		public static void QuickSortIndex(int[] idx, float[][] x, int iCol, int left, int right)
		{
			int i = left, j = right;
			int tmp;
			double pivot = x[idx[(left + right) / 2]][iCol];

			//	Partition
			while (i <= j)
			{
				while (x[idx[i]][iCol] < pivot) i++;
				while (x[idx[j]][iCol] > pivot) j--;
				if (i <= j)
				{
					tmp = idx[i];
					idx[i++] = idx[j];
					idx[j--] = tmp;
				}
			}

			//	Recursion
			if (left < j)
				QuickSortIndex(idx, x, iCol, left, j);
			if (i < right)
				QuickSortIndex(idx, x, iCol, i, right);
		}
    }
}