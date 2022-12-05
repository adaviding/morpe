using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using Morpe.Validation;

namespace Morpe.Numerics.D
{
	/// <summary>
	/// Static operations and static readonlyants for double-precision floating point numbers.
	/// </summary>
	public class Util
	{
		public static readonly double Pi_1_4   = 0.7853981633974483;
		public static readonly double Pi_1_2   = 1.570796326794897;
		public static readonly double Pi_3_4   = 2.356194490192345;
		public static readonly double Pi_5_4   = 3.926990816987241;
		public static readonly double Pi_3_2   = 4.712388980384689;
		public static readonly double Pi_7_4   = 5.497787143782138;
		public static readonly double Pi_2     = 6.283185307179586;
		public static readonly double Sqrt_1_2 = 0.7071067811865475;
		public static readonly double Sqrt_2   = 1.414213562373095;

		/// <summary>
		/// Computes the Cholesky factor of a symmetric, positive semidefinite matrix.
		/// </summary>
		/// <param name="matrix">The matrix to be factored: A symmetric, positive semidefinite matrix.  Only the
		/// upper-right triangle is used (including the diagonal).</param>
		/// <returns>If the matrix was positive semidefinite, then the Cholesky factor is returned as a lower-left
		/// triangular matrix, where the upper-right triangle is zero; otherwise null is returned.</returns>
		[return: MaybeNull]
		public static double[,] CholeskyFactor([NotNull] double[,] matrix)
		{
			Chk.NotNull(matrix, nameof(matrix));
			
			int n = matrix.GetLength(0);
			Chk.Equal(matrix.GetLength(1), n, "The matrix must be square.");
			
			double[,] output = new double[n, n];
			double sum;
			
			for (int i = 0; i < n ; i++)
			{
				for (int j = i ;j < n; j++)
				{
					sum = matrix[i,j];
					for (int k = i - 1;k >= 0;k--)
						sum -= matrix[i,k] * matrix[j,k];
					if (i == j)
					{
						if (sum <= 0.0)
							return null;  // Not positive semidefinite.

						output[i,i] = Math.Sqrt(sum);
					}
					else
						output[j,i] = sum / output[i,i];
				}
			}
			
			return output;
		}
		
		/// <summary>
		/// Generates a random n * n rotation matrix.  This is useful for generating a random orthonormal basis.
		/// </summary>
		/// <param name="n">The size of the matrix.</param>
		/// <returns>The random rotation matrix.</returns>
		public static double[,] RandomRotationMatrix(int n)
		{
			int i, j, k;
			double theta, c, s, z;

			//	Initialize the identity matrix.
			double[,] output = new double[n, n];
			for (i = 0; i < n; i++)
				output[i, i] = 1.0;
			
			//	Special case, n==1.
			if (n == 1)
			{
				if (Morpe.Util.Rand.NextDouble() < 0.5)
					output[0, 0] = -1.0;
				return output;
			}

			//	For each pair (i,j) of rows in the rotation matrix
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					if (i != j)
					{
						//	Random rotation angle.
						theta = Morpe.Util.Rand.NextDouble() * Pi_2; // planar rotation counterclockwise by theta
						c = Math.Cos(theta);	// R(i,i) and  R(j,j)
						s = Math.Sin(theta);	// R(j,i) and -R(i,j)

						//	For each column R(:,i) and R(:,j), multiply each row M(i,:) and M(j,:).  In-place storage (in M) is possible.
						for (k = 0; k < n; k++)
						{
							//	X(i,k)	=	R(i,i)*M(i,k)	+	R(i,j)*M(j,k)
							z = c * output[i,k] - s * output[j,k];
							//	X(j,k)	=	R(j,i)*M(i,k)	+	R(j,j)*M(j,k)
							output[j,k] = s * output[i,k] + c * output[j,k];
							//	In place storage.
							output[i,k] = z;
						}
					}
				}
			}
			return output;
		}
	}
}
