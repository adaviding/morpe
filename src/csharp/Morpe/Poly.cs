using System;
using System.Diagnostics.CodeAnalysis;
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
		/// Returns the number of inhomogeneous polynomial coefficients for a given dimensionality and rank, including
		/// all coefficients of lesser rank.
		/// </summary>
		/// <param name="numDims"><see cref="NumDims"/>.</param>
		/// <param name="rank"><see cref="Rank"/>.</param>
		public static int NumCoeff(int numDims, int rank)
		{
			return I.Util.Pascal(numDims, rank);

			//				0	1	2	3	4	5	6	7	8	9
			//							Rank of Tensor
			//	0			1	1	1	1	1	1	1	1	1	1
			//	1			1	2	3	4	5	6	7	8	9	10
			//	2	Ndims	1	3	6	10	15	21	28	36	45	55
			//	3	of		1	4	10	20	35	56	84	120	165	220
			//	4	Space	1	5	15	35	70	126	210	330	495	715
			//	5			1	6	21	56	126	252	462	924	1716
		}
		
		/// <summary>
		/// Returns the number of inhomogeneous polynomial coefficients for a given dimensionality and rank, NOT
		/// including coefficients of lesser rank.
		/// </summary>
		/// <param name="nDims"><see cref="NumDims"/>.</param>
		/// <param name="rank"><see cref="Rank"/>.</param>
		public static int NumCoeffAtRank(int nDims, int rank)
		{
			return I.Util.Pascal(nDims-1, rank);
		}
		
		/// <summary>
		/// The polynomial coefficients.  A 2D array, right-padded with -1 values.  The -1 values signal "empty".  Each
		/// row is non-decreasing and defines a coefficient.
		///
		///	For example...
		///		A row of {0, 0, 2, 5, -1, -1, -1} would correspond to the following polynomial term.
		///			x[0] * x[0] * x[2] * x[5]
		///		A row of {0, 1, 2, -1, -1, -1, -1} would correspond to the following polynomial term.
		///			x[0] * x[1] * x[2]
		/// 
		/// The number of rows is <see cref="NumCoeffs"/>.
		/// The number of columns is <see cref="NumDims"/>.
		/// </summary>
		[NotNull]
		public int [,] Coeffs { get; private set; }
		
		/// <summary>
		/// The number of polynomial coefficients.  Also, the number of rows of <see cref="Poly.Coeffs"/> .
		/// </summary>
		public readonly int NumCoeffs;
		
		/// <summary>
		/// The number of spatial dimensions over which the polynomial is defined.  (i.e. The number of coordinate axes.)
		/// </summary>
		public readonly int NumDims;
		
		/// <summary>
		///	The rank of the polynomial.  (i.e. The maximum power of a term.)
		///
		///	For example...
		///		A linear polynomial would be Rank=1.
		///		A quadratic polynomial would be Rank=2.
		///		A cubic polynomial would be Rank=3.
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

			//	How many polynomial coefficients total?
			this.NumCoeffs = NumCoeff(numDims, rank);

			//	Allocate coefficients.
			this.Coeffs = new int[this.NumCoeffs,numDims];

			//	One coefficient (or row) at a time.
			int[] row = new int[numDims];
			for(int i=0; i<row.Length; i++)
				row[i] = -1;

			//	For each coefficient.
			for(int iCoeff=0; iCoeff < this.NumCoeffs; iCoeff++)
			{
				int i, j;
				//	Increment the smallest digit and begin to handle "carry-over" arithmetic.
				i=-1;
				while( ++row[++i]==numDims );	//	When a digit is maxed, keep incrementing rightward until we're not maxed out.

				//	Finish the "carry-over" by ensuring that any leftward maxed digits have been reset.
				for( j=i-1; j>=0; j-- )
				{
					if( row[j]==numDims )
						row[j]=row[j+1];
				}

				//	Copy last row to output
				for(j=0; j<numDims; j++)
					this.Coeffs[iCoeff,j] = row[j];
			}
		}

		/// <summary>
		/// Creates a deep copy of the instance.
		/// </summary>
		/// <returns>The deep copy.</returns>
		public Poly Clone()
		{
			Poly output = this.MemberwiseClone() as Poly;

			output.Coeffs = this.Coeffs.Clone() as int[,];

			return output;
		}
		
		/// <summary>
		/// Computes 'output' as the polynomial expansion of 'x'.
		/// </summary>
		/// <param name="output">OUTPUT:  The expanded feature vector.  (Memory is already pre-allocated.)</param>
		/// <param name="x">INPUT:  The feature vector.</param>
		public void Expand(float[] output, float[] x)
		{
			if(x.Length < this.NumDims)
				throw new ArgumentException("The argument x cannot have fewer elements than the number of spatial dimensions.");
			if (output.Length < this.NumCoeffs)
				throw new ArgumentException ("The argument y cannot have fewer elements than the number of polynomial coefficients.");
			int iExp, iRank;
			float yy;
			for(iExp=0; iExp<this.NumCoeffs; iExp++)
			{
				yy = 1.0f;
				for (iRank=0; iRank<this.Rank; iRank++)
				{
					int ix = this.Coeffs[iExp,iRank];
					if(ix<0) break;
					yy *= x[ix];
				}
				output[iExp] = yy;
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
