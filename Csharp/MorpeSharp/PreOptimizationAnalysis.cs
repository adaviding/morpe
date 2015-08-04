using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe
{
	/// <summary>
	/// This encapsulates measurements of the conditioned and expanded training data made prior to optimiztion.  The purpose of these measurements is to
	/// provide good initial estimates of the parameters to be optimized.  Also, some of these measurements are used to constrain the optimization algorithm.
	/// </summary>
	public class PreOptimizationAnalysis
	{
		/// <summary>
		/// The space conditioner.
		/// </summary>
		public SpatialConditioner Conditioner;
		/// <summary>
		/// Measures the spatial condition.  Contains extra information.
		/// </summary>
		public SpatialConditionMeasurer ConditionMeasurer;
		/// <summary>
		/// Unidimensional accuracy-maximizing criteria.
		/// </summary>
		public UniCrit[][] Crits;
		/// <summary>
		/// The initial values of the polynomial coefficients.  Indexed as [iPoly][iCoeff].
		/// </summary>
		public float[][] ParamInit;
		/// <summary>
		/// The scale of each polynomial coefficient.  The scales are the same across all polynomials.
		/// </summary>
		public float[] ParamScale;
		/// <summary>
		/// The scale of the training data (for each column of data).
		/// </summary>
		public float[] Xscale;
	}
}
