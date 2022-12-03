using System;

namespace Morpe
{
	/// <summary>
	/// Represents Morpe's options, and options that adjust the solver's behavior.
	/// </summary>
	public class SolverOptions
	{
		/// <summary>
		/// This is called to construct a new instance of the <see cref="Morpe.Options"/> class.
		/// </summary>
		public SolverOptions() { }
		/// <summary>
		/// Optimization halts when the change in conditional entropy is less than this amount.  This number must be positive.  Default = 0.005.
		/// </summary>
		public float EntropyTol=0.005f;
		/// <summary>
		/// If true, the solver will pick initial values for the classifier's parameters.
		/// If false, the classifier's parameters will be used as initial parameters.
		/// </summary>
		public bool InitializeParams = true;
		/// <summary>
		/// The optimization routine will run this many times.  Each time, it will approach the solution along a different path.
		/// </summary>
		public int NumberOfApproaches=5;
		/// <summary>
		///	The rate at which parameters shrink during optimization.  Must be greater than 1.0.
		/// </summary>
		public float ParamShrinkFactor=2.0f;
		/// <summary>
		/// Must be greater than 0.  Optimization cannot change parameters less than the following magnitude for any [p,i]-th parameter.
		//		[MAGNITUDE] = ParamDiffMax * Morpe.Trainer.ParamScale[p][i]
		/// </summary>
		public float ParamDiffMax = 2.0f;
		/// <summary>
		/// Must be greater than 0.  Optimization cannot change parameters less than the following magnitude for any [p,i]-th parameter.
		//		[MAGNITUDE] = ParamDiffMin * Morpe.Trainer.ParamScale[p][i]
		/// </summary>
		public float ParamDiffMin = 0.005f;
		/// <summary>
		/// Determines how the optimization routine weights data from each category.
		/// </summary>
		public WeightingRule WeightingRule = WeightingRule.EqualPriors;
		public SolverOptions Copy()
		{
			SolverOptions output = new SolverOptions();
			output.EntropyTol = this.EntropyTol;
			output.InitializeParams = this.InitializeParams;
			output.NumberOfApproaches = this.NumberOfApproaches;
			output.ParamShrinkFactor = this.ParamShrinkFactor;
			output.ParamDiffMax = this.ParamDiffMax;
			output.ParamDiffMin = this.ParamDiffMin;
			output.WeightingRule = this.WeightingRule;
			return output;
		}
	}
}