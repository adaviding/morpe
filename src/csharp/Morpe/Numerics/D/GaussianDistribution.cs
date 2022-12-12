using System.Diagnostics.CodeAnalysis;

namespace Morpe.Numerics.D
{
    public static class GaussianDistribution
    {
        /// <summary>
        /// Calculates the probability density of the coordinate 'x' with respect to a Gaussian distribution having the
        /// specified properties. 
        /// </summary>
        /// <param name="x">The spatial coordinate for which the z-score is calculated.</param>
        /// <param name="mean">The mean of the Gaussian distribution.</param>
        /// <param name="invChol">The inverse of the Cholesky factor of the covariance matrix for the Gaussian
        /// distribution.</param>
        /// <returns>The probability density.</returns>
        public static double Density(
            [NotNull] double[] x,
            [NotNull] double[] mean,
            [NotNull] double[,] invChol)
        {
            double z = Zscore(x, mean, invChol);
            double output = D1.GaussianDistribution.Pdf(z);
            return output;
        }

        /// <summary>
        /// Calculate the z-score of the coordinate 'x' with respect to a Gaussian distribution having the specified
        /// properties.
        /// </summary>
        /// <param name="x">The spatial coordinate for which the z-score is calculated.</param>
        /// <param name="mean">The mean of the Gaussian distribution.</param>
        /// <param name="invChol">The inverse of the Cholesky factor of the covariance matrix for the Gaussian
        /// distribution.</param>
        /// <returns>The z-score.</returns>
        public static double Zscore(
            [NotNull] double[] x,
            [NotNull] double[] mean,
            [NotNull] double[,] invChol)
        {
            // Subtract the mean
            double[] xc = D.Util.AddScaled(-1.0, mean, x);
            
            // Divide by the standard deviation.
            double[] z = D.Util.Product(xc, invChol);
            
            // Length of 'z'
            double output = D.Util.NormL2(z);

            return output;
        }
    }
}
