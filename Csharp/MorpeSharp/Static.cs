using System;

namespace Morpe
{
	/// <summary>
	/// Static methods of general purpose.
	/// </summary>
	public class Static
	{
		/// <summary>
		/// Copies the values from one grid of values to another.
		/// </summary>
		/// <param name="src">The source grid.</param>
		/// <param name="dest">The destination grid.</param>
		public static void Copy(float[][] src, float[][] dest)
		{
			if (src == null || dest == null)
				return;
			if (src.Length > dest.Length)
				throw new ArgumentException("The destination length must be at least equal to the source length");

			for (int iRow = 0; iRow < src.Length; iRow++)
			{
				float[] x = src[iRow];
				if(x!=null)
				{
					if (dest[iRow] == null || dest[iRow].Length < x.Length)
						dest[iRow] = new float[x.Length];
					Array.Copy(x, dest[iRow], x.Length);
				}
			}
		}
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
		/// Initializes an array of arrays simulating a 2D arrow with nRows rows and nCols columns.
		/// </summary>
		/// <typeparam name="T">The array element type.</typeparam>
		/// <param name="nRows">The number of rows.</param>
		/// <param name="nCols">The number of columns.</param>
		/// <returns>An array of arrays.</returns>
		public static T[][] NewArrays<T>(int nRows, int nCols)
		{
			T[][] output = new T[nRows][];
			for (int iRow = 0; iRow < nRows; iRow++)
				output[iRow] = new T[nCols];
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
		/// Returns the entry (a,b) of the Pascal matrix (which is symmetric).  (a+b)! / a / b
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
		/// <summary>
		/// Calculates the sum of a vector of integers.
		/// </summary>
		/// <param name="vals">The vector of integers.</param>
		public static int Sum(int[] vals)
		{
			int output = 0;
			foreach (int v in vals)
				output += v;
			return output;
		}
		/// <summary>
		/// Calculates the sum of a vector of double.
		/// </summary>
		/// <param name="vals">The vector of double.</param>
		public static double Sum(double[] vals)
		{
			double output = 0.0;
			foreach (double v in vals)
				output += v;
			return output;
		}
	}
}