using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Morpe.Validation;

using D = Morpe.Numerics.D;

namespace Morpe
{
    /// <summary>
    /// Basic Gaussian statistics per category, basic descriptions of <see cref="CategorizedData"/>.
    /// </summary>
    public class GaussianStats
    {
        /// <summary>
        /// Measures Gaussian stats for each category.
        /// </summary>
        /// <param name="data"></param>
        /// <returns></returns>
        [return: NotNull]
        public static GaussianStats Measure([NotNull] CategorizedData data)
        {
            GaussianStats output = new GaussianStats();
            output.Means = new double[data.NumCats][];
            output.Covs = new double[data.NumCats][,];

            for (int i = 0; i < data.NumCats; i++)
            {
                output.Means[i] = MeasureMean(data.X[i]);
                output.Covs[i] = MeasureCov(data.X[i], output.Means[i]);
            }

            return output;
        }
        
        /// <summary>
        /// The covariance matrix for each category.
        /// </summary>
        public double[][,] Covs;

        /// <summary>
        /// The number of categories.
        /// </summary>
        public int NumCats
        {
            get
            {
                return this.Means?.Length ?? -1;
            }
        }

        /// <summary>
        /// The number of spatial dimensions;
        /// </summary>
        public int NumDims
        {
            get
            {
                return this.Means?[0]?.Length ?? -1;
            }
        }

        /// <summary>
        /// The mean of each category;
        /// </summary>
        public double[][] Means;
        
        /// <summary>
        /// In a scenario where <see cref="NumCats"/> is greater than 2, we sometimes create a "dual" view of the data
        /// categories where 'targetCat' is category 0, and all other categories are combined into category 1.
        /// </summary>
        /// <param name="targetCat">The target category.</param>
        /// <param name="weights">The relative weighting for each category.  If null, each category is given equal
        /// weighting (like with <see cref="CategoryWeightingRule.EqualPriors"/>).</param>
        /// <returns>The <see cref="GaussianStats"/> for the "dual" view of the target category.</returns>
        public GaussianStats ForDual(int targetCat, [MaybeNull] double[] weights)
        {
            // Cache some values for computational efficiency.
            int numCats = this.NumCats;
            int numDims = this.NumDims;
            
            Chk.Less(1, numCats, "There must be at least 2 categories.");
            Chk.Less(0, numDims, "There must be at least 1 spatial dimension.");
            Chk.Less(targetCat, numCats,
                "The target category {0} is greater or equal to the number of categories {1}.", targetCat, numCats);
            Chk.LessOrEqual(0, targetCat, "The target category {0} cannot be negative.", targetCat);
            
            // Weights
            if (weights == null)
            {
                // Assume equal priors.
                weights = new double[numCats];
                for (int i = 0; i < numCats; i++)
                {
                    weights[i] = 1.0;
                }
            }
            else
            {
                // Validation on the weights that got passed in.
                Chk.Equal(numCats, weights.Length, "There must be 1 weight per category.");
                Chk.True(weights.All(a => a >= 0.0), "At least one negative weight was supplied.");
            }
            double weightOfNonTargets = weights.Sum() - weights[targetCat];
            Chk.Less(0, weightOfNonTargets, "The total weight of non-target categories must be greater than zero.");
            
            GaussianStats output = new GaussianStats
            {
                Means = new double[2][],
                Covs = new double[2][,]
            };
            
            // The target category is easy.
            output.Means[0] = this.Means[targetCat].Clone() as double[];
            output.Covs[0] = this.Covs[targetCat].Clone() as double[,];

            if (numCats == 2)
            {
                // This is a degenerate case.  Just return the non-target as category "1".
                output.Means[1] = this.Means[1 - targetCat].Clone() as double[];
                output.Covs[1] = this.Covs[1 - targetCat].Clone() as double[,];

                return output;
            }

            // Calculate the mean.
            output.Means[1] = new double[numDims];
            for (int iCat = 0; iCat < numCats; iCat++)
            {
                if (iCat == targetCat)
                    continue;
                
                double relativeWeight = weights[iCat] / weightOfNonTargets;
                double[] mean = this.Means[iCat];
                for (int iDim = 0; iDim < numDims; iDim++)
                {
                    output.Means[1][iDim] += mean[iDim] * relativeWeight;
                }
            }
            
            // Calculate the difference between the category means from the grand mean.
            double[][] meanDiffs = new double[numCats][];
            for (int iCat = 0; iCat < numCats; iCat++)
            {
                if (iCat == targetCat)
                    continue;
                
                meanDiffs[iCat] = D.Util.Subtract(this.Means[iCat], output.Means[1]);
            }

            // Calculate the covariance matrix.
            double[,] covCombined = new double[numDims, numDims];
            output.Covs[1] = covCombined;

            for (int iCat = 0; iCat < numCats; iCat++)
            {
                if (iCat == targetCat)
                    continue;
                
                double relativeWeight = weights[iCat] / weightOfNonTargets;
                double[] meanDiff = meanDiffs[iCat];
                double[,] cov = this.Covs[iCat];
                
                for (int iDim = 0; iDim < numDims; iDim++)
                {
                    for (int jDim = iDim; jDim < numDims; jDim++)
                    {
                        // Fill in the diagonal and upper triangle.
                        covCombined[iDim, jDim] += relativeWeight * (cov[iDim, jDim] + meanDiff[iDim] * meanDiff[jDim]);
                    }
                }
            }
            
            for (int iDim = 0; iDim < numDims; iDim++)
            {
                for (int jDim = iDim + 1; jDim < numDims; jDim++)
                {
                    // Reflect the upper triangle onto the lower triangle.
                    covCombined[jDim, iDim] = covCombined[iDim, jDim];
                }
            }

            return output;
        }

        /// <summary>
        /// Gets the quadratic parameters which are optimal for a 2-category problem.
        ///
        /// This returns all first-order (linear) and second-order (quadratic) terms of the polynomial equation.  That
        /// is 1 less parameter than the Fisher Quadratic Discriminant which contains an additional scalar parameter.
        /// MoRPE does not need that scalar parameter so we don't compute it.
        ///
        /// This calculation is related to Appendix B in the MoRPE paper:  In this case we are computing the parameters
        /// for '-y' (negative "y") because it is larger for category 0 and smaller for category 1.
        /// </summary>
        /// <returns>The parameters of a quadratic function which are optimal for a Gaussian assumption.  This is
        /// indexed in a way that matches the definition of <see cref="Poly.Coeffs"/>.
        /// </returns>
        /// <remarks>This is typically used to calculate starting parameters for a MoRPE classifier.</remarks>
        [return: NotNull]
        public double[] QuadraticParams([NotNull] D.MatrixInvertor matrixInvertor)
        {
            Chk.Equal(2, this.NumCats, "This method can only be called for the 2-category problem.");
            Chk.Equal(this.NumDims, matrixInvertor.Rank,
                "The rank of the matrix invertor ({0}) must be equal to the number of spatial dimensions ({1}).",
                matrixInvertor.Rank,
                this.NumDims);

            // Cache this value for computational efficiency.
            int numDims = this.NumDims;
            
            // Calculate the number of output parameters and initialize the output.
            int numParams = numDims + (numDims * (numDims + 1)) / 2;
            double[] output = new double[numParams];

            // Invert the covariance matrices.
            double[][,] invCovs = new double[2][,]
            {
                matrixInvertor.Invert(this.Covs[0]),
                matrixInvertor.Invert(this.Covs[1])
            };
            
            // Compute 'b' the vector of 1-order params.
            double[] b = D.Util.Subtract(
                D.Util.Product(this.Means[1], invCovs[1]),
                D.Util.Product(this.Means[0], invCovs[0]));
            
            // Compute 'C' the matrix of 2-order params
            //double[,] c = 0.5 * D.Util.Subtract(invCovs[0], invCovs[1]);
            
            // Fill the output in a manner consistent with Poly.Coeffs.
            int iOutput = 0;
            
            // The linear parameters
            for (int i = 0; i < numDims; i++)
            {
                output[iOutput++] = -b[i];
            }
            
            // The quadratic params
            for (int i = 0; i < numDims; i++)
            {
                for (int j = i; j < numDims; j++)
                {
                    if (i == j)
                    {
                        output[iOutput++] = 0.5 * (invCovs[1][i, j] - invCovs[0][i, j]);
                    }
                    else
                    {
                        output[iOutput++] =
                            0.5 * (invCovs[1][i, j] - invCovs[0][i, j] + invCovs[1][j, i] - invCovs[0][j, i]);
                    }
                }
            }

            return output;
        }

        /// <summary>
        /// Measure the covariance matrix.
        /// </summary>
        /// <param name="x">Rows of data.</param>
        /// <param name="mean">The mean of the rows data.</param>
        /// <returns>The covariance matrix.</returns>
        private static double[,] MeasureCov([NotNull] float[][] x, [NotNull] double[] mean)
        {
            Chk.Less(1, x.Length, "There is not enough data.");
            int numDims = x[0].Length;
            Chk.Less(0, numDims, "The spatial dimensionality is 0.");
            
            double[,] output = new double[numDims, numDims];
            foreach (float[] row in x)
            {
                for (int iDim = 0; iDim < numDims; iDim++)
                {
                    double xi = row[iDim] - mean[iDim];
                    output[iDim, iDim] += xi * xi;

                    for (int jDim = iDim + 1; jDim < numDims; jDim++)
                    {
                        double xj = row[jDim] - mean[jDim];
                        output[iDim, jDim] += xi * xj;
                    }
                }
            }

            double scalar = 1.0 / (x.Length - 1);
            for (int iDim = 0; iDim < numDims; iDim++)
            {
                output[iDim, iDim] *= scalar;
                for (int jDim = iDim + 1; jDim < numDims; jDim++)
                {
                    output[iDim, jDim] *= scalar;
                    output[jDim, iDim] = output[iDim, jDim];
                }
            }

            return output;
        }

        /// <summary>
        /// Measure the mean.
        /// </summary>
        /// <param name="x">Rows of data;</param>
        /// <returns>The mean.</returns>
        private static double[] MeasureMean([NotNull] float[][] x)
        {
            Chk.Less(0, x.Length, "There is no data.");
            int numDims = x[0].Length;
            Chk.Less(0, numDims, "The spatial dimensionality is 0.");

            double[] output = new double[numDims];
            foreach (float[] row in x)
            {
                for (int iDim = 0; iDim < numDims; iDim++)
                {
                    output[iDim] += row[iDim];
                }
            }
            
            for (int iDim = 0; iDim < numDims; iDim++)
            {
                output[iDim] /= x.Length;
            }

            return output;
        }
    }
}
