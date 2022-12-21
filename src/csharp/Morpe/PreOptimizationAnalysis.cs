using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Morpe.Validation;

using F = Morpe.Numerics.F;

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
        public UniCrit[,] Crits;
        
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
        
        /// <summary>
        /// Gets a <see cref="PreOptimizationAnalysis"/> for the given subspace.
        /// </summary>
        /// <param name="fullPoly">The polygon that describes the full space (the space of this analysis).</param>
        /// <param name="subPoly">The polygon that describes the subspace.</param>
        /// <param name="subDims">For each dimension of the subspace polygon, this gives the corresponding index into
        /// the fullspace polygon.  The values must be increasing.</param>
        /// <returns>The analysis for the given subspace.</returns>
        [return: NotNull]
        public PreOptimizationAnalysis Subspace(
            [NotNull] Poly fullPoly,
            [NotNull] Poly subPoly,
            [NotNull] int[] subDims)
        {
            Chk.Less(subPoly.Rank, this.Rank,
                "The rank of the sub-polygon should be less than this instance.");
            Chk.Less(subPoly.Rank, fullPoly.Rank,
                "The rank of the subspace polygon should be less than the fullspace polygon.");
            Chk.Equal(subDims.Length, subPoly.NumDims,
                "{0}.{1} != {2}.{3}",
                nameof(subDims), nameof(subDims.Length),
                nameof(subPoly), nameof(subPoly.NumDims));
            Chk.Increasing(subDims, "The vector of subspace dimension indices must be increasing.");
            Chk.LessOrEqual(0, subDims.First(), "The vector of subspace dimension indices is out of range (too low).");
            Chk.Less(subDims.Last(), fullPoly.NumDims, "The vector of subspace dimension indices is out of range (too high).");

            PreOptimizationAnalysis output = this.Clone();

            int numPoly = this.Crits.Length;
            output.Rank = subPoly.Rank;
            output.Crits = new UniCrit[numPoly, subPoly.NumCoeffs];
            output.ParamScale = new float[subPoly.NumCoeffs];
            output.ParamScaleNorm = new float[subPoly.Rank];

            // Allocate space for initial parameter values.
            output.ParamInit = Util.NewArrays<float[]>(subPoly.Rank, numPoly);
            for (int iRank = 0; iRank < subPoly.Rank; iRank++)
                for (int iPoly = 0; iPoly < numPoly; iPoly++)
                    output.ParamInit[iRank][iPoly] = new float[subPoly.NumCoeffsForRank[iRank]];

            // Map subspace coefficients to fullspace coefficients.
            int[] mapSubToFull = Poly.SubspaceToFullspaceCoefficientMapping(fullPoly, subPoly, subDims);

            // For each subspace coefficient
            for (int iCoeff = 0; iCoeff < subPoly.NumCoeffs; iCoeff++)
            {
                int jCoeff = mapSubToFull[iCoeff];
                
                output.ParamScale[iCoeff] = this.ParamScale[jCoeff];

                for (int iPoly = 0; iPoly < numPoly; iPoly++)
                    output.Crits[iPoly, iCoeff] = this.Crits[iPoly, jCoeff];
                
                for(int iRank=subPoly.Coeffs[iCoeff].Length-1; iRank<output.Rank; iRank++)
                    for (int iPoly = 0; iPoly < numPoly; iPoly++)
                        output.ParamInit[iRank][iPoly][iCoeff] = this.ParamInit[iRank][iPoly][jCoeff];
            }
            
            for(int iRank=0; iRank<output.Rank; iRank++)
                output.ParamScaleNorm[iRank] = (float)F.Util.NormL2(output.ParamInit[iRank]);

            return output;
        }
    }
}
