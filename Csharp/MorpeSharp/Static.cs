using System;

namespace Morpe
{
	/// <summary>
	/// Static methods of general purpose.
	/// </summary>
	public class Static
	{
		/// <summary>
		/// Performs a binary search on a tabulated function.
		/// </summary>
		/// <param name="xIncreasing">The tabulated function.  The values must be non-decreasing.</param>
		/// <param name="xTarget">The target value.</param>
		/// <returns>Returns the zero-based index of the target value within the range of the tabulated function.
		/// It can be a continuous (non-integer) value, in the range [0.0, xTabulated.Length], representing an exact
		/// placement between two tabulated values (as a linear-interpolant between those values).  It can also be
		/// -1.0 if below the tabulated range, or xTabulated.Length if above the range.</returns>
		public static double BinarySearch_NonDecreasing(double[] xTabulated, double xTarget)
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
				//	Finish
				//-----------------------------------------------------------------------------
				if (xTabulated[iMid] == xTarget)
				{
					//	Put both indexes in range and rely on execution of next "if" statement to finish.
					iLow = iHigh = iMid;
				}
				if (iHigh - iLow < 5)
				{
					//	Get iLow to be the lowest index of xTabulated possibly equal to the xTarget.
					//	Otherwise it should be just below the xTarget.
					while (iLow > iMin && xTabulated[iLow - 1] >= xTarget) iLow--;
					while (iLow < iMax && xTabulated[iLow + 1] < xTarget) iLow++;
					if (iLow < iMax && xTabulated[iLow + 1] == xTarget) iLow++;
					//	Get iHigh to be the highest index of xTabulated possibly equal to the xTarget.
					//	Otherwise it should be just above the xTarget.
					while (iHigh < iMax && xTabulated[iHigh + 1] <= xTarget) iHigh++;
					while (iHigh > iMin && xTabulated[iHigh - 1] > xTarget) iHigh--;
					if (iHigh > iMin && xTabulated[iHigh - 1] == xTarget) iHigh--;

					if (iHigh == iLow) // One unique xTarget value.  Return its index.
						return (double)iLow;
					if (xTabulated[iLow] == xTarget) // Many unique xTarget values exist.  Return the center of the range.
						return 0.5 * (double)(iLow + iHigh);
					if (iHigh - iLow != 1)
						//	Return linear interpolant.
						return (double)iLow + (xTarget - xTabulated[iLow]) / (xTabulated[iHigh] - xTabulated[iLow]);
				}

				//-----------------------------------------------------------------------------
				//	Continue
				//-----------------------------------------------------------------------------
				//  Reduce possible range by 1/2
				if (xTabulated[iMid] > xTarget)
					iLow = iMid;
				else // Must be < xTarget
					iHigh = iMid;
			}
		}
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
		/// Performs linear interpolation of a tabulated function y -> f(x).
		/// </summary>
		/// <param name="xTabulated">Non-decreasing tabulated values of an independent variable "x".</param>
		/// <param name="yTabulated">Tabulated values of a dependent variable "y".</param>
		/// <param name="xTarget">The value of x for which the matching value of y is linearly interpolated.</param>
		/// <returns>The interpolated value of the function y -> f(x) where x=xTarget.</returns>
		public static double Linterp(double[] xTabulated, double[] yTabulated, double xTarget)
		{
			double i = BinarySearch_NonDecreasing(xTabulated, xTarget);
			if (i < 0.0 || i >= (double)(xTabulated.Length - 1))
				return double.NaN;
			double iMod = i - (double)(int)i;
			if (iMod < 1e-10 || iMod > 0.9999999999)
				return yTabulated[(int)(i + 0.5)];
			else
				return (1.0 - iMod) * yTabulated[(int)i] + iMod * yTabulated[1 + (int)i];
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
		public static readonly Random Rand = new Random();
		/// <summary>
		/// Generates a random n * n rotation matrix.  This is useful for generating a random orthonormal basis.
		/// </summary>
		/// <param name="n">The size of the matrix.</param>
		/// <returns>The random rotation matrix.</returns>
		public static double[,] RandomRotationMatrix(int n)
		{
			int i, j, k;
			double theta, c, s, z;

			//	Initialize the identity matrix.
			double[,] output = new double[n, n];
			for (i = 0; i < n; i++)
				output[i, i] = 1.0;
			
			//	Special case, n==1.
			if (n == 1)
			{
				if (Rand.NextDouble() < 0.5)
					output[0, 0] = -1.0;
				return output;
			}

			//	For each pair (i,j) of rows in the rotation matrix
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					if (i != j)
					{
						//	Random rotation angle.
						theta = Rand.NextDouble() * D.TwoPi; // planar rotation counterclockwise by theta
						c = Math.Cos(theta);	// R(i,i) and  R(j,j)
						s = Math.Sin(theta);	// R(j,i) and -R(i,j)

						//	For each column R(:,i) and R(:,j), multiply each row M(i,:) and M(j,:).  In-place storage (in M) is possible.
						for (k = 0; k < n; k++)
						{
							//	X(i,k)	=	R(i,i)*M(i,k)	+	R(i,j)*M(j,k)
							z = c * output[i,k] - s * output[j,k];
							//	X(j,k)	=	R(j,i)*M(i,k)	+	R(j,j)*M(j,k)
							output[j,k] = s * output[i,k] + c * output[j,k];
							//	In place storage.
							output[i,k] = z;
						}
					}
				}
			}
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