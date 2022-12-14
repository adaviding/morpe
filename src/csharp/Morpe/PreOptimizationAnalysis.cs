using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe
{
	/// <summary>
	/// This encapsulates measurements of the conditioned and expanded training data made prior to optimization.  The purpose of these measurements is to
	/// provide good initial estimates of the parameters to be optimized.  Also, some of these measurements are used to constrain the optimization algorithm.
	/// </summary>
	public class PreOptimizationAnalysis
	{
		/// <summary>
		/// The spatial conditioner.
		/// </summary>
		public SpatialConditioner Conditioner;
		
		/// <summary>
		/// Measures the spatial conditioner based on training data.  Contains statistical information.
		/// </summary>
		public SpatialConditionMeasurer ConditionMeasurer;
		
		/// <summary>
		/// Unidimensional accuracy-maximizing criteria for each conditioned-expanded dimension.
		/// Indexed as [iCat][iCoeff]
		/// </summary>
		public UniCrit[][] Crits;
		
		/// <summary>
		/// The initial values of the polynomial coefficients.  Indexed as [iPoly][iCoeff].
		/// </summary>
		public float[][] ParamInit;
		
		/// <summary>
		/// The scale of data for each conditioned-expanded dimension.  These scales are the same
		/// across all polynomials.
		/// </summary>
		public float[] ParamScale;
		
		/// <summary>
		/// Deep copy.
		/// </summary>
		/// <returns>A deep copy.</returns>
		public PreOptimizationAnalysis Copy()
		{
			PreOptimizationAnalysis output = new PreOptimizationAnalysis();
			
			if(this.Conditioner!=null)
				output.Conditioner = this.Conditioner.Clone();

			if(this.ConditionMeasurer!=null)
				output.ConditionMeasurer = this.ConditionMeasurer.Copy();
			
			if(this.Crits!=null)
			{
				output.Crits = Util.Copy<UniCrit>(this.Crits);
				for(int iCat=0; iCat<output.Crits.Length; iCat++)
				{
					UniCrit[] crits = output.Crits[iCat];
					for(int iCoeff=0; iCoeff<crits.Length; iCoeff++)
						crits[iCoeff] = crits[iCoeff].Clone();
				}
			}
			
			if(this.ParamInit!=null)
				output.ParamInit = Util.Copy<float>(this.ParamInit);
			
			if(this.ParamScale!=null)
				output.ParamScale = (float[])this.ParamScale.Clone();
			
			return output;
		}
		
		/// <summary>
		/// Creates a deep copy of this instance which is cast for a 2 category sub-problem where the target category
		/// is cast as category 0 and the remaining samples are cast as category 1.
		/// </summary>
		/// <param name="iCat">The target category.</param>
		/// <returns>A copy of the instance, cast for the dual problem.</returns>
		public PreOptimizationAnalysis CopyAsDual(int iCat)
		{
			if(iCat<0)
				return null;

			PreOptimizationAnalysis output = new PreOptimizationAnalysis();
			
			if (this.Conditioner != null)
				output.Conditioner = this.Conditioner.Clone();

			if (this.Crits != null && this.Crits.Length > iCat)
			{
				output.Crits = new UniCrit[1][];
				UniCrit[] crits = output.Crits[iCat];
				for (int iCoeff = 0; iCoeff < crits.Length; iCoeff++)
					crits[iCoeff] = crits[iCoeff].Clone();
				output.Crits[0] = crits;
			}
			
			if (this.ParamInit != null && this.ParamInit.Length > iCat)
			{
				output.ParamInit = new float[1][];
				output.ParamInit[0] = (float[])this.ParamInit[iCat].Clone();
			}
			
			if (this.ParamScale != null)
				output.ParamScale = (float[])this.ParamScale.Clone();
			
			return output;
		}
	}
}
