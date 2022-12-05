using System;

namespace Morpe.Numerics.F1
{
	public class F
	{
		/// <summary>
		/// Adds the elements of "input" to the elements of "output".
		/// </summary>
		/// <param name="output">The matrix where summation is compiled.</param>
		/// <param name="scalar">The addends.  The values added onto the output.</param>
		public static void Add(float[][] output, float[][] input)
		{
			if (input == null)
				return;
			if (output == null)
				throw new ArgumentException("Output cannot be null.");
			if (output.Length != input.Length)
				throw new ArgumentException("The arrays must be of equal size.");
			for (int iRow = 0; iRow < output.Length; iRow++)
			{
				float[] outRow = output[iRow];
				float[] inRow = input[iRow];
				if (outRow != null && inRow != null)
				{
					if (outRow.Length != inRow.Length)
						throw new ArgumentException("The arrays must be of equal size.");
					for (int iCol = 0; iCol < outRow.Length; iCol++)
						outRow[iCol] += inRow[iCol];
				}
			}
		}

		/// <summary>
		/// Multiplies the elements of a matrix "output" by a scalar.
		/// </summary>
		/// <param name="output">The matrix.</param>
		/// <param name="scalar">The scalar</param>
		public static void Multiply(float[][] output, float scalar)
		{
			if (output == null)
				return;
			for (int iRow = 0; iRow < output.Length; iRow++)
			{
				float[] row = output[iRow];
				if (row != null)
				{
					for (int iCol = 0; iCol < row.Length; iCol++)
						row[iCol] *= scalar;
				}
			}
		}
	}
}