using System.Collections.Generic;
using System.Text;
using Morpe;
using NUnit.Framework;

namespace Test.Morpe
{
    public class PolynomialTests
    {
        [TestCase(1, 1, 1, 1)]
        [TestCase(2, 1, 2, 2)]
        [TestCase(2, 2, 5, 3)]
        [TestCase(3, 1, 3, 3)]
        [TestCase(3, 2, 9, 6)]
        [TestCase(3, 3, 19, 10)]
        [TestCase(3, 4, 34, 15)]
        [TestCase(4, 1, 4, 4)]
        [TestCase(4, 2, 14, 10)]
        public void TestNumCoeffs(int numDims, int rank, int numCoeffs, int numCoeffsHomo)
        {
            Assert.AreEqual(numCoeffs,     Poly.NumCoeff(numDims, rank),     "Number of inhomogeneous coefficients is wrong.");
            Assert.AreEqual(numCoeffsHomo, Poly.NumCoeffHomo(numDims, rank), "Number of homogeneous coefficients is wrong.");
        }

        [TestCase(1,1)]
        [TestCase(2,1)]
        [TestCase(3,1)]
        [TestCase(4,1)]
        [TestCase(2,2)]
        [TestCase(3,2)]
        [TestCase(4,2)]
        [TestCase(2,3)]
        [TestCase(3,3)]
        [TestCase(4,3)]
        [TestCase(2,4)]
        [TestCase(3,4)]
        [TestCase(4,4)]
        public void TestPolyCoeffs(int numDims, int rank)
        {
            Poly poly = new Poly(numDims, rank);
            
            Assert.AreEqual(numDims, poly.NumDims, nameof(poly.NumDims));
            Assert.AreEqual(rank, poly.Rank, nameof(poly.Rank));

            int numCoeff = Poly.NumCoeff(numDims, rank);
            Assert.AreEqual(numCoeff, poly.Coeffs.Length, "Wrong number of coefficients.");

            // We are going to build a unique string to represent each coefficient, to ensure there are no duplicates.
            StringBuilder uniqueCoeffStringBuilder = new StringBuilder();
            HashSet<string> uniqueCoffStrings = new HashSet<string>();

            for (int i = 0; i < poly.NumCoeffs; i++)
            {
                int len = poly.Coeffs[i].Length;
                Assert.LessOrEqual(len, rank, "The length of the coefficient vector cannot be greater than the polynomial rank.");
                
                // The prior value of 'current' (defined below).
                int prior = -1; // initialize with a dummy value.
                
                for (int j = 0; j < len; j++)
                {
                    int current = poly.Coeffs[i][j];

                    Assert.GreaterOrEqual(current, 0, "Poly coeffs must be non-negative.");
                    Assert.Less(current, numDims, "Poly coeffs must be less than the number of spatial dimensions.");
                    
                    if (i > 0)
                    {
                        Assert.LessOrEqual(prior, current, "Poly coeffs should be non-decreasing.");
                        uniqueCoeffStringBuilder.Append("_");
                    }

                    uniqueCoeffStringBuilder.Append(current);

                    prior = current;
                }

                // Make sure that each coefficient is unique.
                string uniqueCoeffString = uniqueCoeffStringBuilder.ToString();
                uniqueCoeffStringBuilder.Clear();
                Assert.False(uniqueCoffStrings.Contains(uniqueCoeffString),
                    "The coefficient {0} has at least 1 duplicate.  Coefficients must be unique.", uniqueCoeffString);
                uniqueCoffStrings.Add(uniqueCoeffString);
            }
            
        }
    }
}
