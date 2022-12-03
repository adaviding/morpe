using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe
{
	public enum MonotonicRegressionType
	{
		/// <summary>
		/// A monotonic regression where the output probably has "flat spots".  Thus, the output is not truly monotonic.
		/// </summary>
		NonDecreasing,
		/// <summary>
		/// A monotonic regression where the flat spots of a non-decreasing function are linearly interpolated.  The output is monotonic.
		/// </summary>
		Increasing,
		/// <summary>
		/// This is the recommended setting.  It is a blending of two other methods:  NonDecreasing and Increasing.  The output is monotonic.
		/// </summary>
		Blended
	}
}
