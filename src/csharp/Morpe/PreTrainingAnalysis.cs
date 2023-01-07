using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using Morpe.Validation;

using F = Morpe.Numerics.F;
using F1 = Morpe.Numerics.F1;
using I = Morpe.Numerics.I;

namespace Morpe
{
    /// <summary>
    /// This encapsulates measurements of the conditioned and expanded training data made prior to optimization.  The purpose of these measurements is to
    /// provide good initial estimates of the parameters to be optimized.  Also, some of these measurements are used to constrain the optimization algorithm.
    /// </summary>
    public class PreTrainingAnalysis
    {
        /// <summary>
        /// Performs the pre-optimization analysis based on conditioned and expanded data.
        /// </summary>
        /// <param name="data">Training data that has already been conditioned and expanded, where the <see cref="CategorizedData.State"/> encapsulates
        /// those operations.</param>
        /// <returns>The analysis.</returns>
        public static PreTrainingAnalysis BuildFromExpandedData(
            CancellationToken cancellationToken,
            [NotNull] CategorizedData data)
        {
            Chk.NotNull(data, nameof(data));
            Chk.NotNull(data.State, "{0}.{1}", nameof(data), nameof(data.State));
            
            PreTrainingAnalysis output = new PreTrainingAnalysis();

            CategoryWeights weights = CategoryWeights.Measure(
                numEach: data.NumEach,
                rule: CategoryWeightingRule.EqualPriors,
                targetCategory: null);
            
            CategoryWeights[] dualWeights = data.NumCats <= 2
                ? null
                : CategoryWeights.MeasureAllDuals(
                    numEach: data.NumEach,
                    rule: CategoryWeightingRule.EqualPriors);

            int numPolys = data.NumCats == 2
                ? 1
                : data.NumCats;

            output.Rank = data.State.Polynomial.Rank;
            output.ParamScale = new float[data.State.Polynomial.NumCoeffs];
            output.Crits = new UniCrit[numPolys, data.State.Polynomial.NumCoeffs];

            output.ParamInit = new float[output.Rank][][];
            for (int i = 0; i < output.Rank; i++)
                output.ParamInit[i] = Util.NewArrays<float>(numPolys, data.State.Polynomial.NumCoeffsForRank[i]);
            
            int i20 = (int)(0.5 + 0.20 * (double)(data.NumTotal - 1));
            int i80 = (int)(0.5 + 0.80 * (double)(data.NumTotal - 1));
            
            float[] xVec = new float[data.NumTotal];   // A single coordinate for each datum.
            int[] iVec = new int[data.NumTotal];       // An index for each datum.
            byte[] catVec = data.GetCategoryVector();  // The category label for each datum.

            // For each polynomial coefficient (i.e. each column of the expanded data)
            for (int iCoeff = 0; iCoeff < data.State.Polynomial.NumCoeffs; iCoeff++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                
                // Fill the iVec.  We will be sorting 'x' through this index so that we can measure quantiles.
                I.Util.FillSeries(iVec);
                
                // Fill the xVec.
                int iDatum = 0;
                for (int iCat = 0; iCat < data.NumCats; iCat++)
                {
                    int catLen = data.NumEach[iCat];
                    for (int jDatum = 0; jDatum < catLen; jDatum++)
                        xVec[iDatum++] = data.X[iCat][jDatum][iCoeff];
                }
                
                // Sort the iVec.
                F.Util.QuickSortIndex(iVec, xVec, 0, xVec.Length-1);
                
                // Measure the spread which gives us the parameter scale.
                float spread = xVec[iVec[i80]] - xVec[iVec[i20]];
                output.ParamScale[iCoeff] = 1f / spread;
                
                // Measure norms for the parameter scale (for each polynomial rank).
                output.ParamScaleNorm = new float[output.Rank];
                for (int iRank = 0; iRank < output.Rank; iRank++)
                {
                    output.ParamScaleNorm[iRank] = (float)Norm(output.ParamScale, data.State.Polynomial.NumCoeffsForRank[iRank]);
                }
                
                // Get a univariate classification criteria.
                for (int iPoly = 0; iPoly < numPolys; iPoly++)
                {
                    output.Crits[iPoly, iCoeff] = UniCrit.MaximumAccuracy(
                        targetCat: (byte)iPoly,
                        cat: catVec,
                        x: xVec,
                        idx: iVec,
                        catWeight: numPolys == 1
                            ? weights
                            : dualWeights[iPoly]);
                }
            }
            
            // At this point, we have done everything besides starting parameters.  In order to calculate these, we will
            // use the univariate criteria for each dimension, and weight each parameter explosively by classification accuracy.
            double chanceAcc = 0.5;
            double invDataLength = 1.0 / data.NumTotal;
            
            // For each polynomial coefficient (i.e. each column of the expanded data)
            for (int iCoeff = 0; iCoeff < data.State.Polynomial.NumCoeffs; iCoeff++)
            {
                for (int iPoly = 0; iPoly < numPolys; iPoly++)
                {
                    UniCrit crit = output.Crits[iPoly, iCoeff];
                    double accAboveChance = Math.Max(0.0, crit.Accuracy - chanceAcc);
                    double weight = accAboveChance / (1.0 - chanceAcc - invDataLength);
                    double x = weight * output.ParamScale[iCoeff]; 

                    if (!crit.TargetUpper)
                    {
                        x = -x;
                    }

                    // Fill all the ranks
                    int iRank = output.Rank;
                    while (--iRank >= 0 && iCoeff < data.State.Polynomial.NumCoeffsForRank[iRank])
                    {
                        output.ParamInit[iRank][iPoly][iCoeff] = (float)x;
                    }
                }
            }
            
            // A final adjustment to parameter scale.
            for (int iRank = 0; iRank < output.Rank; iRank++)
            {
                // Basically the norm of the parameters should match the norm of the param scale vector.
                F.Util.Multiply(
                    output: output.ParamInit[iRank],
                    scalar: (float)(output.ParamScaleNorm[iRank] / F.Util.NormL2(output.ParamInit[iRank])));
            }

            return output;
        }
        
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
        public PreTrainingAnalysis Clone()
        {
            PreTrainingAnalysis output = (PreTrainingAnalysis)this.MemberwiseClone();
            
            output.Conditioner = this.Conditioner?.Clone();
            output.ConditionMeasurer = this.ConditionMeasurer?.Clone();
            output.Crits = Util.Clone(this.Crits);
            output.ParamInit = Util.Clone(this.ParamInit);
            output.ParamScale = (float[])this.ParamScale?.Clone();
            output.ParamScaleNorm = (float[])this.ParamScaleNorm?.Clone();
            
            return output;
        }
        
        /// <summary>
        /// Gets a <see cref="PreTrainingAnalysis"/> for the given subspace.
        /// </summary>
        /// <param name="fullPolynomial">The polynomial that describes the full space (the space of this analysis).</param>
        /// <param name="subPolynomial">The polynomial that describes the subspace.</param>
        /// <param name="subDims">For each dimension of the subspace polynomial, this gives the corresponding index into
        /// the fullspace polynomial.  The values must be increasing.</param>
        /// <returns>The analysis for the given subspace.</returns>
        [return: NotNull]
        public PreTrainingAnalysis Subspace(
            [NotNull] Polynomial fullPolynomial,
            [NotNull] Polynomial subPolynomial,
            [NotNull] int[] subDims)
        {
            Chk.Less(subPolynomial.Rank, this.Rank,
                "The rank of the sub-polynomial should be less than this instance.");
            Chk.Less(subPolynomial.Rank, fullPolynomial.Rank,
                "The rank of the subspace polynomial should be less than the fullspace polynomial.");
            Chk.Equal(subDims.Length, subPolynomial.NumDims,
                "{0}.{1} != {2}.{3}",
                nameof(subDims), nameof(subDims.Length),
                nameof(subPolynomial), nameof(subPolynomial.NumDims));
            Chk.Increasing(subDims, "The vector of subspace dimension indices must be increasing.");
            Chk.LessOrEqual(0, subDims.First(), "The vector of subspace dimension indices is out of range (too low).");
            Chk.Less(subDims.Last(), fullPolynomial.NumDims, "The vector of subspace dimension indices is out of range (too high).");

            PreTrainingAnalysis output = this.Clone();

            int numPoly = this.Crits.Length;
            output.Rank = subPolynomial.Rank;
            output.Crits = new UniCrit[numPoly, subPolynomial.NumCoeffs];
            output.ParamScale = new float[subPolynomial.NumCoeffs];
            output.ParamScaleNorm = new float[subPolynomial.Rank];

            // Allocate space for initial parameter values.
            output.ParamInit = Util.NewArrays<float[]>(subPolynomial.Rank, numPoly);
            for (int iRank = 0; iRank < subPolynomial.Rank; iRank++)
                for (int iPoly = 0; iPoly < numPoly; iPoly++)
                    output.ParamInit[iRank][iPoly] = new float[subPolynomial.NumCoeffsForRank[iRank]];

            // Map subspace coefficients to fullspace coefficients.
            int[] mapSubToFull = Polynomial.SubspaceToFullspaceCoefficientMapping(fullPolynomial, subPolynomial, subDims);

            // For each subspace coefficient
            for (int iCoeff = 0; iCoeff < subPolynomial.NumCoeffs; iCoeff++)
            {
                int jCoeff = mapSubToFull[iCoeff];
                
                output.ParamScale[iCoeff] = this.ParamScale[jCoeff];

                for (int iPoly = 0; iPoly < numPoly; iPoly++)
                    output.Crits[iPoly, iCoeff] = this.Crits[iPoly, jCoeff];
                
                for(int iRank=subPolynomial.Coeffs[iCoeff].Length-1; iRank<output.Rank; iRank++)
                    for (int iPoly = 0; iPoly < numPoly; iPoly++)
                        output.ParamInit[iRank][iPoly][iCoeff] = this.ParamInit[iRank][iPoly][jCoeff];
            }
            
            for(int iRank=0; iRank<output.Rank; iRank++)
                output.ParamScaleNorm[iRank] = (float)F.Util.NormL2(output.ParamInit[iRank]);

            return output;
        }

        /// <summary>
        /// Calculates the L2 norm of the first 'len' elements of 'x'.
        /// </summary>
        /// <param name="x">The vector.</param>
        /// <param name="len">The number of elements that contribute to the norm.</param>
        /// <returns>The L2 norm of the first 'len' elements of 'x'.</returns>
        private static double Norm([NotNull] float[] x, int len)
        {
            double output = 0.0;

            for (int i = 0; i < len; i++)
            {
                double d = x[i];
                output += d * d;
            }

            output = Math.Sqrt(output);
            
            return output;
        }
    }
}
