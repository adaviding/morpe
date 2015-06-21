using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe
{
	public class MonotonicRegression
	{
		/// <summary>
		/// Performs a monotonic regression of a tabulated function.  This method calculates a monotonic function that is approximately equal
		/// to the tabulated function provided.
		/// </summary>
		/// <param name="output">The monotonic function.  This should be supplied as vector of length identical to input.
		/// When the function has finished executing, the following are guaranteed to be true:
		///		mean(output) == mean(input)
		///		min(output) >= min(input)
		///		max(output) <= max(input)
		/// </param>
		/// <param name="input">The tabulated function.</param>
		/// <returns>The number of repetitions through a loop.</returns>
		public int Run(float[] output, float[] input)
		{
			if (input == null || output == null || output.Length < input.Length)
				return 0;

			throw new NotImplementedException();
		}
	}
	public enum MonotonicRegressionType
	{
		/// <summary>
		/// A monotonic regression where the output is likely to have "flat spots".  Thus, the output is not truly monotonic.
		/// </summary>
		NonDecreasing,
		/// <summary>
		/// A monotonic regression where the flat spots of a non-decreasing function are linearly interpolated.  The output is monotonic.
		/// </summary>
		Increasing,
		/// <summary>
		/// A blending of two other methods:  NonDecreasing and Increasing.
		/// </summary>
		Blended
	}
}
