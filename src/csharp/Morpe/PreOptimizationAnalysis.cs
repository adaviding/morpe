using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Morpe.Validation;

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
        /// Indexed as [iPoly][iCoeff]
        /// </summary>
        public UniCrit[][] Crits;
        
        /// <summary>
        /// The initial values of the polynomial coefficients.  Indexed as [iRank][iPoly][iCoeff].
        /// </summary>
        public float[][][] ParamInit;
        
        /// <summary>
        /// The scale of each polynomial coefficient.  This is inversely proportional to the spread of the data for
        /// each polynomial coefficient.
        ///
        /// The length of this vector is equal to the number of polynomial coefficients for the rank <see cref="Rank"/>.
        /// </summary>
        public float[] ParamScale;

        /// <summary>
        /// The norm of <see cref="ParamScale"/> for each iRank rank, where iRank is rank - 1.
        /// </summary>
        public float[] ParamScaleNorm;

        /// <summary>
        /// The maximum polynomial rank which is covered by this analysis.
        /// </summary>
        public int Rank;
        
        /// <summary>
        /// Creates a deep copy.
        /// </summary>
        /// <returns>The deep copy.</returns>
        public PreOptimizationAnalysis Clone()
        {
            PreOptimizationAnalysis output = (PreOptimizationAnalysis)this.MemberwiseClone();
            
            output.Conditioner = this.Conditioner?.Clone();
            output.ConditionMeasurer = this.ConditionMeasurer?.Clone();
            output.Crits = Util.Clone(this.Crits);
            output.ParamInit = Util.Clone(this.ParamInit);
            output.ParamScale = (float[])this.ParamScale?.Clone();
            output.ParamScaleNorm = (float[])this.ParamScaleNorm?.Clone();
            
            return output;
        }
    }
}
