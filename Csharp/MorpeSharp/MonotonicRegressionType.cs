using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe
{
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
