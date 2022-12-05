using System;
using System.Runtime.CompilerServices;

namespace Morpe
{
	/// <summary>
	/// Static methods of general purpose.
	/// </summary>
	public static class Util
	{
		/// <summary>
		/// Copies an array of arrays.
		/// </summary>
		/// <typeparam name="T">The array element type.</typeparam>
		/// <param name="src">An array of arrays.</param>
		/// <returns>An array of arrays with values equal to the input values.</returns>
		public static T[][] Copy<T>(T[][] src)
		{
			if (src == null)
				return null;
			T[][] output = new T[src.Length][];
			for (int iRow = 0; iRow < src.Length; iRow++)
			{
				T[] row = src[iRow];
				if (row != null)
				{
					output[iRow] = new T[row.Length];
					for (int iCol = 0; iCol < row.Length; iCol++)
						output[iRow][iCol] = row[iCol];
				}
			}
			return output;
		}
		
		/// <summary>
		/// Copies the values from one grid of values to another.
		/// </summary>
		/// <param name="src">The source grid.</param>
		/// <param name="dest">The destination grid.</param>
		public static void Copy<T>(T[][] src, T[][] dest)
		{
			if (src == null || dest == null)
				return;
			if (src.Length > dest.Length)
				throw new ArgumentException("The destination length must be at least equal to the source length");

			for (int iRow = 0; iRow < src.Length; iRow++)
			{
				T[] x = src[iRow];
				if (x != null)
				{
					if (dest[iRow] == null || dest[iRow].Length < x.Length)
						dest[iRow] = new T[x.Length];
					Array.Copy(x, dest[iRow], x.Length);
				}
			}
		}

		/// <summary>
		/// Initializes an array of arrays simulating a 2D array with nRows rows and nCols columns.
		/// </summary>
		/// <typeparam name="T">The array element type.</typeparam>
		/// <param name="rows">The number of rows.</param>
		/// <param name="cols">The number of columns.</param>
		/// <returns>An array of arrays.</returns>
		public static T[][] NewArrays<T>(int rows, int cols)
		{
			T[][] output = new T[rows][];
			for (int iRow = 0; iRow < rows; iRow++)
				output[iRow] = new T[cols];
			return output;
		}
		/// <summary>
		/// Initializes an array of arrays simulating a 3D array with nPages pages, nRows rows, and nCols columns.
		/// </summary>
		/// <typeparam name="T">The array element type.</typeparam>
		/// <param name="pages">The number of pages.</param>
		/// <param name="rows">The number of rows.</param>
		/// <param name="cols">The number of columns.</param>
		/// <returns>An array of arrays of arrays.</returns>
		public static T[][][] NewArrays<T>(int pages, int rows, int cols)
		{
			T[][][] output = new T[pages][][];
			for (int iPage = 0; iPage < pages; iPage++)
				output[iPage] = NewArrays<T>(rows, cols);
			return output;
		}
		
		/// <summary>
		/// Gets the random number generator for the current thread.
		/// </summary>
		public static Random Rand
		{
			get
			{
				if (rand == null)
				{
					// In .NET Framework, the clock has a limited resolution.  To avoid having the same seed on multiple
					// threads we subtract the managed thread id.
					rand = new Random((int)DateTime.UtcNow.Ticks - Environment.CurrentManagedThreadId);
				}

				return rand;
			}
		}
		[ThreadStatic]
		private static Random rand;
		
		/// <summary>
		/// Sets all elements of "output" equal to a value.
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="output">The array whose elements are set.</param>
		/// <param name="value">The value.</param>
		public static void SetValues<T>(T[][] output, T value)
		{
			if (output == null)
				return;
			foreach(T[] row in output)
			{
				if(row!=null)
				{
					for (int i = 0; i < output.Length; i++)
						row[i] = value;
				}
			}
		}
		
		/// <summary>
		/// Shuffles the pages in a random order.
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="x">Pages of sheets.</param>
		public static void Shuffle<T>(T[][][] x)
		{
			T[][] sheet;
			for(int i=0; i<x.Length; i++)
			{
				int j = Rand.Next(x.Length);
				if(j!=i)
				{
					sheet = x[i];
					x[i] = x[j];
					x[j] = sheet;
				}
			}
		}
	}
}