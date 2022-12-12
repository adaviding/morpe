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
        /// The categorized data.
        /// </summary>
        public CategorizedData Data;

        private MockGaussianDataSet()
        {
        }

        /// <summary>
        /// Creates a mock Gaussian data set with the given number of categories and dimensions.  The category means
        /// are the columns of a random rotation matrix.  The covariance matrices are identity.
        /// </summary>
        /// <param name="numCats"><see cref="NumCats"/></param>
        /// <param name="numDims"><see cref="NumDims"/></param>
        public MockGaussianDataSet(int numCats, int numDims)
        {
            this.NumCats = numCats;
            this.NumDims = numDims;
            
            Chk.Less(1, numCats, "There must be at least 2 categories.");
            Chk.Less(0, numCats, "There must be at least 1 spatial dimension.");

            this.Means = new double[numCats][];
            this.Covs = new double[numCats][,];
            this.InvChols = new double[numCats][,];

            double[,] randRot = D.Util.RandomRotationMatrix(Math.Max(numCats, numDims));
            
            for (int c = 0; c < this.NumCats; c++)
            {
                this.Means[c] = new double[numDims];

                for (int i = 0; i < this.NumDims; i++)
                {
                    this.Means[c][i] = randRot[c, i];
                    this.Covs[c] = D.Util.IdentityMatrix(numDims);
                    this.InvChols[c] = D.Util.IdentityMatrix(numDims);
                }
            }
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
                    output.X[iCat][iRow] = F.Util.From(x);
                }
            }

            return output;
        }
    }
}