using System;
using System.Diagnostics.CodeAnalysis;
using Morpe;
using Morpe.Validation;
using NUnit.Framework.Constraints;
using D = Morpe.Numerics.D;
using D1 = Morpe.Numerics.D1;
using F = Morpe.Numerics.F;

namespace Test.Morpe
{
    /// <summary>
    /// A mock data set based on Gaussian populations for which ground truth is known.
    /// </summary>
    public class MockGaussianDataSet
    {
        /// <summary>
        /// Creates some spherical distributions having pairwise d' ("d prime" from signal detection theory) as specified. 
        /// </summary>
        /// <param name="numEach">The number of samples in each category.</param>
        /// <param name="numDims">The number of spatial dimensions.  This must be at least numEach.Length because of the pairwise d' constraint.</param>
        /// <param name="dPrime">The value of pairwise d' ("d prime" from signal detection theory).  So any two categories (selected at random) will be separated
        /// by this much d'.</param>
        /// <returns>The mock Gaussian data set.</returns>
        public static MockGaussianDataSet CreateSpherical(
            [NotNull] int[] numEach,
            int numDims,
            double dPrime)
        {
            Chk.NotNull(numEach, nameof(numEach));
            Chk.Less(1, numEach.Length, "There must be at least 2 categories.");
            Chk.Less(0, numDims, "There must be at least 1 spatial dimension.");
            Chk.LessOrEqual(0.0, dPrime, "{0} should be non-negative.", nameof(dPrime));
            Chk.LessOrEqual(numEach.Length, numDims, "The number of spatial dimensions cannot be less than the number of categories because of the pairwise d' constraint.");

            MockGaussianDataSet output = new MockGaussianDataSet();
            output.NumCats = numEach.Length;
            output.NumDims = numDims;
            
            output.Means = new double[output.NumCats][];
            output.Covs = new double[output.NumCats][,];
            output.InvChols = new double[output.NumCats][,];

            double[,] randRot = D.Util.RandomRotationMatrix(numDims);

            double meanScalar = dPrime / Math.Sqrt(numDims);
            
            for (int c = 0; c < output.NumCats; c++)
            {
                output.Means[c] = new double[numDims];

                for (int i = 0; i < output.NumDims; i++)
                {
                    output.Means[c][i] = randRot[c, i] * meanScalar;
                    output.Covs[c] = D.Util.IdentityMatrix(numDims);
                    output.InvChols[c] = D.Util.IdentityMatrix(numDims);
                }
            }

            output.Data = output.CreateRandomSample(numEach);

            return output;
        }
        
        /// <summary>
        /// The number of categories.
        /// </summary>
        public int NumCats;
        
        /// <summary>
        /// The number of spatial dimensions.
        /// </summary>
        public int NumDims;

        /// <summary>
        /// The mean for each category.  Indexed as [c][i].
        /// c = Category
        /// i = Spatial dimension
        /// </summary>
        public double[][] Means;

        /// <summary>
        /// The covariance matrix for each category.  Indexed as [c][i,j].
        /// c = Category
        /// i = Spatial dimension
        /// j = Spatial dimension
        /// </summary>
        public double[][,] Covs;
        
        /// <summary>
        /// The inverse of the cholesky factor of the covariance matrix for each category.  Indexed as [c][i,j].
        /// c = Category
        /// i = Spatial dimension
        /// j = Spatial dimension
        ///
        /// These are lower-triangular matrices, as described in <see cref="D.Util.CholeskyFactor"/>
        /// </summary>
        public double[][,] InvChols;

        /// <summary>
        /// The categorized data.  This is typically set as the output of <see cref="CreateRandomSample"/>.
        /// </summary>
        public CategorizedData Data;

        public MockGaussianDataSet()
        {
        }

        /// <summary>
        /// Creates a random sample of data.
        /// </summary>
        /// <param name="numEach">The number of samples in each category.  The length of this vector must be equal
        /// to <see cref="NumCats"/>.</param>
        /// <returns>The random sample of data.</returns>
        [return: NotNull]
        public CategorizedData CreateRandomSample(
            [NotNull] int[] numEach)
        {
            Chk.Equal(this.NumCats, numEach.Length, "The length of the input must be equal to the number of categories.");

            CategorizedData output = new CategorizedData(
                numEach: numEach,
                numDims: this.NumDims);

            for (int iCat = 0; iCat < numEach.Length; iCat++)
            {
                double[,] chol = D.Util.CholeskyFactor(this.Covs[iCat]);
                for (int iRow = 0; iRow < numEach[iCat]; iRow++)
                {
                    double[] z = D1.GaussianDistribution.Rand(this.NumDims);
                    double[] x = D.Util.Add(D.Util.Product(chol, z), this.Means[iCat]);
                    output.X[iCat][iRow] = F.Util.Convert(x);
                }
            }

            return output;
        }

        /// <summary>
        /// Measures the optimal fit for the <see cref="Data"/>.  This can be measured because we know the properties of the Gaussian populations.
        /// </summary>
        /// <param name="weights">The category weights.</param>
        /// <returns>The optimal fit (accuracy and entropy).</returns>
        public (double accuracy, double entropy) MeasureOptimalFit([NotNull] CategoryWeights weights)
        {
            return this.MeasureOptimalFit(weights, this.Data);
        }

        /// <summary>
        /// Measures the optimal fit for the given data.  This can be measured because we know the properties of the Gaussian populations.
        /// </summary>
        /// <param name="weights">The category weights.</param>
        /// <param name="data">The given data.</param>
        /// <returns>The optimal fit (accuracy and entropy).</returns>
        public (double accuracy, double entropy) MeasureOptimalFit(
            [NotNull] CategoryWeights weights,
            [NotNull] CategorizedData data)
        {
            Chk.NotNull(weights, nameof(weights));
            Chk.NotNull(data, nameof(data));

            (double accuracy, double entropy) output = (0.0, 0.0);

            double totalWeight = 0.0;
            
            for (int iCat = 0; iCat < data.NumCats; iCat++)
            {
                double weight = weights.Weights[iCat]; 
                    
                for (int iDatum = 0; iDatum < data.X[iCat].Length; iDatum++)
                {
                    int catMax = -1;
                    double pdMax = double.MinValue;
                    double pdCat = 0.0;
                    double pdNotCat = 0.0;
                    double pd;

                    float[] x = data.X[iCat][iDatum];
                    
                    for (int jCat = 0; jCat < data.NumCats; jCat++)
                    {
                        // Measure the probability density for the current category.
                        pd = D.GaussianDistribution.Density(
                            x: D.Util.Convert(x),
                            mean: this.Means[jCat],
                            invChol: this.InvChols[jCat]);

                        if (pd > pdMax)
                        {
                            pdMax = pd;
                            catMax = jCat;
                        }

                        if (jCat == iCat)
                        {
                            pdCat += pd;
                        }
                        else
                        {
                            pdNotCat += pd;
                        }
                    }

                    if (iCat == catMax)
                    {
                        output.accuracy += weight;
                    }

                    pd = pdCat / (pdCat + pdNotCat);
                    output.entropy += -weight * Math.Log(pd, data.NumCats);
                    
                    totalWeight += weight;
                }
            }

            output.accuracy /= totalWeight;
            output.entropy /= totalWeight;

            return output;
        }
    }
}
