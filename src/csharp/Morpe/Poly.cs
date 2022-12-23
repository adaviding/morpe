using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Runtime.InteropServices;
using Morpe.Validation;
using I = Morpe.Numerics.I;

namespace Morpe
{
    /// <summary>
    /// Represents a multivariate inhomogeneous polynomial.
    /// <see cref="http://mathworld.wolfram.com/MultivariatePolynomial.html"/>
    /// </summary>
    public class Poly
    {
        /// <summary>
        /// Computes a mapping from subspace polynomial coefficients up to fullspace polynomial coefficients.
        /// </summary>
        /// <param name="fullPoly">The fullspace polynomial.</param>
        /// <param name="subPoly">The subspace polynomial.</param>
        /// <param name="subDims">For each spatial dimension of the subspace polynomial, this gives the index of the
        /// corresponding spatial dimension in the fullspace polynomial.</param>
        /// <returns>For each coefficient of the subspace polynomial, this gives the corresponding index of the coefficient
        /// in the fullspace polynomial.</returns>
        [return: NotNull]
        public static int[] SubspaceToFullspaceCoefficientMapping(
            [NotNull] Poly fullPoly,
            [NotNull] Poly subPoly,
            [NotNull] int[] subDims)
        {
            Chk.Less(subPoly.Rank, fullPoly.Rank,
                "The rank of the sub-polynomial should be less than this instance.");
            Chk.Less(subPoly.Rank, fullPoly.Rank,
                "The rank of the subspace polynomial should be less than the fullspace polynomial.");
            Chk.Equal(subDims.Length, subPoly.NumDims,
                "{0}.{1} != {2}.{3}",
                nameof(subDims), nameof(subDims.Length),
                nameof(subPoly), nameof(subPoly.NumDims));
            Chk.Increasing(subDims, "The vector of subspace dimension indices must be increasing.");
            Chk.LessOrEqual(0, subDims.First(), "The vector of subspace dimension indices is out of range (too low).");
            Chk.Less(subDims.Last(), fullPoly.NumDims, "The vector of subspace dimension indices is out of range (too high).");
            
            int[] output = new int[subPoly.NumCoeffs];

            int jCoeff = 0;
            for(int iCoeff=0; iCoeff<subPoly.NumCoeffs; iCoeff++)
            {
                int[] sub = subPoly.Coeffs[iCoeff];
                int[] mapping = new int[sub.Length];
                for(int iRank=0; iRank < sub.Length; iRank++)
                    mapping[iRank] = subDims[sub[iRank]];
                while (!Test.EqualListing(fullPoly.Coeffs[jCoeff], mapping))
                    jCoeff++;
                output[iCoeff] = jCoeff;
                jCoeff++;
            }

            return output;
        }
        
        /// <summary>
        /// Returns the number of inhomogeneous polynomial coefficients for a given dimensionality and rank, including
        /// all coefficients of lesser rank.  However, it does not include the 0-th order term which MoRPE does not use.
        /// </summary>
        /// <param name="numDims"><see cref="NumDims"/>.</param>
        /// <param name="rank"><see cref="Rank"/>.</param>
        public static int NumCoeff(int numDims, int rank)
        {
            // This is the number of coefficients including the 0-th order term.
            int output = I.Util.Pascal(numDims, rank);

            //                0    1    2    3    4    5    6    7    8    9
            //                            Rank of Tensor
            //    0            1    1    1    1    1    1    1    1    1    1
            //    1            1    2    3    4    5    6    7    8    9    10
            //    2    Ndims   1    3    6    10    15    21    28    36    45    55
            //    3    of      1    4    10   20    35    56    84    120   165   220
            //    4    Space   1    5    15   35    70    126   210   330   495   715
            //    5            1    6    21   56    126   252   462   924   1716
            
            // We subtract 1 to remove the 0-th order term.  (MoRPE does not use this term.)
            output--;

            return output;
        }
        
        /// <summary>
        /// Returns the number of homogeneous polynomial coefficients for a given dimensionality and rank.
        ///
        /// This does NOT including coefficients of lesser rank (as per the term "homogeneous").
        /// </summary>
        /// <param name="nDims"><see cref="NumDims"/>.</param>
        /// <param name="rank"><see cref="Rank"/>.</param>
        public static int NumCoeffHomo(int nDims, int rank)
        {
            return I.Util.Pascal(nDims-1, rank);
        }
        
        /// <summary>
        /// The inhomogeneous polynomial coefficients.  Each row is non-decreasing and defines a coefficient.
        ///
        ///    For example...
        ///        A row of {0, 0, 2, 5} would correspond to the following polynomial term.
        ///            x[0] * x[0] * x[2] * x[5]
        ///        A row of {0, 1, 2} would correspond to the following polynomial term.
        ///            x[0] * x[1] * x[2]
        /// 
        /// The number of rows is <see cref="NumCoeffs"/>.
        /// The number of values in each row is equal to the rank of that coefficient which is no greater than <see cref="Rank"/>. 
        /// </summary>
        [NotNull]
        public int [][] Coeffs { get; private set; }
        
        /// <summary>
        /// The number of inhomogeneous polynomial coefficients.  Also, the number of rows of <see cref="Poly.Coeffs"/> .
        /// </summary>
        public readonly int NumCoeffs;
        
        /// <summary>
        /// The number of inhomogeneous polynomial coefficients at the given rank.
        ///
        /// [rank-1] The index of zero is actually a rank of 1.  Subtract 1 from the rank to get the index.
        /// </summary>
        [NotNull]
        public int[] NumCoeffsForRank { get; private set; }
        
        /// <summary>
        /// The number of homogeneous polynomial coefficients at the given rank.
        ///
        /// [rank-1] The index of zero is actually a rank of 1.  Subtract 1 from the rank to get the index.
        /// </summary>
        [NotNull]
        public int[] NumCoeffsForRankHomo { get; private set; }
        
        /// <summary>
        /// The number of spatial dimensions over which the polynomial is defined.  (i.e. The number of coordinate axes.)
        /// </summary>
        public readonly int NumDims;
        
        /// <summary>
        ///    The rank of the polynomial.  (i.e. The maximum power of a term.)
        ///
        ///    For example...
        ///        A linear polynomial would be Rank=1.
        ///        A quadratic polynomial would be Rank=2.
        ///        A cubic polynomial would be Rank=3.
        /// </summary>
        public readonly int Rank;
        
        /// <summary>
        /// Constructs a multivariate polynomial of a given rank and spatial dimensionality.
        /// </summary>
        /// <param name="numDims"><see cref="NumDims"/>.</param>
        /// <param name="rank"><see cref="Poly.Rank"/>.</param>
        public Poly(int numDims, int rank)
        {
            if (numDims < 1)
                throw new ArgumentException("The number of spatial dimensions must be at least 1.");
            if (rank < 1)
                throw new ArgumentException("The polynomial rank must be at least 1.");

            this.NumDims = numDims;
            this.Rank = rank;

            // How many polynomial coefficients total?
            this.NumCoeffs = NumCoeff(numDims, rank);
            
            // Calculate the number of polynomial coefficients for all ranks up to the given rank, and for the 
            // homogeneous and inhomogeneous cases.
            this.NumCoeffsForRank = new int[rank];      // inhomogeneous
            this.NumCoeffsForRankHomo = new int[rank];  // homogeneous

            int numPrior = 0;
            int numCurrent=0;
            for(int iRank=0; iRank < rank; iRank++)
            {
                numPrior = numCurrent;
                numCurrent = Poly.NumCoeff(numDims, iRank + 1);
                this.NumCoeffsForRank[iRank] = numCurrent;
                this.NumCoeffsForRankHomo[iRank] = numCurrent - numPrior;
            }

            //    Allocate coefficients.
            this.Coeffs = new int[this.NumCoeffs][];

            //    One coefficient (or row) at a time.
            List<int> row = new List<int>(rank);
            row.Add(-1);
            
            //int[] row = new int[rank];
            //for(int i=0; i<row.Length; i++)
            //    row[i] = -1;

            //    For each coefficient.
            for(int iCoeff=0; iCoeff < this.NumCoeffs; iCoeff++)
            {
                int i, j;
                
                // Increment the smallest digit and begin to handle "carry-over" arithmetic.
                for (i = 0; i < row.Count; i++)
                {
                    int val = ++row[i];
                    if (val < numDims)
                    {
                        break;
                    }
                    else if (i == row.Count - 1)
                    {
                        i++;
                        row.Add(0);
                    }
                }

                // Finish the "carry-over" by ensuring that any leftward maxed digits have been reset.
                for( j=i-1; j>=0; j-- )
                {
                    if( row[j]==numDims )
                        row[j]=row[j+1];
                }

                this.Coeffs[iCoeff] = ((IEnumerable<int>)row).Reverse().ToArray();
            }
        }

        /// <summary>
        /// Creates a deep copy of the instance.
        /// </summary>
        /// <returns>The deep copy.</returns>
        public Poly Clone()
        {
            Poly output = this.MemberwiseClone() as Poly;

            output.Coeffs = Util.Clone(this.Coeffs);
            output.NumCoeffsForRank = (int[])this.NumCoeffsForRank.Clone();
            output.NumCoeffsForRankHomo = (int[])this.NumCoeffsForRankHomo.Clone();

            return output;
        }
        
        /// <summary>
        /// Computes 'output' as the polynomial expansion of 'x'.
        /// </summary>
        /// <param name="output">OUTPUT:  The expanded feature vector.  (Memory is already pre-allocated.)</param>
        /// <param name="x">INPUT:  The feature vector.</param>
        public void Expand(float[] output, float[] x)
        {
            Chk.Equal(x.Length, this.NumDims, "The length of the feature vector {0} must be equal to the spatial dimensionality {1}.", x.Length, this.NumDims);
            
            if (output.Length < this.NumCoeffs)
                throw new ArgumentException ("The argument y cannot have fewer elements than the number of polynomial coefficients.");
            
            float yy;
            for(int iCoeff=0; iCoeff<this.NumCoeffs; iCoeff++)
            {
                yy = 1.0f;
                int[] thisCoeff = this.Coeffs[iCoeff];
                
                for (int iTerm=0; iTerm<thisCoeff.Length; iTerm++)
                {
                    int ix = this.Coeffs[iCoeff][iTerm];
                    yy *= x[ix];
                }
                output[iCoeff] = yy;
            }
        }
        
        /// <summary>
        /// Computes 'output' as the polynomial expansion of 'x'.
        /// </summary>
        /// <param name="x">The feature vector.</param>
        /// <returns>The expanded feature vector.</returns>
        public float[] Expand(float[] x)
        {
            float[] output = new float[this.NumCoeffs];
            this.Expand(output, x);
            return output;
        }
    }
}
