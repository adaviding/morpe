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
        /// Adds two vectors 'a' and 'b'.
        /// </summary>
        /// <param name="a">The vector 'a'.</param>
        /// <param name="b">The vector 'b'.</param>
        /// <returns>The vector whose elements ar the sum of 'a' and 'b'.</returns>
        public static double[] Add(
            [NotNull] double[] a,
            [NotNull] double[] b)
        {
            Chk.Equal(a.Length, b.Length, "Vectors must be the same length.");

            double[] output = new double[a.Length];

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = a[i] + b[i];
            }

            return output;
        }
        
        /// <summary>
        /// Adds a scaled version of the vector 'a' to the vector 'b' and returns the result.
        /// </summary>
        /// <param name="aScalar">The scalar of 'a'.</param>
        /// <param name="a">The vector 'a'.</param>
        /// <param name="b">The vector 'b'.</param>
        /// <returns>aScalar * a + b</returns>
        [return: NotNull]
        public static double[] AddScaled(
            double aScalar,
            [NotNull] double[] a,
            [NotNull] double[] b)
        {
            Chk.Equal(a.Length, b.Length, "The vectors must be the same length.");

            double[] output = new double[a.Length];

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = aScalar * a[i] + b[i];
            }

            return output;
        }

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
        /// Creates an identity matrix.
        /// </summary>
        /// <param name="n">The rank of the matrix.</param>
        /// <returns>The identity matrix.</returns>
        [return: NotNull]
        public static double[,] IdentityMatrix(int n)
        {
            double[,] output = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                output[i, i] = 1.0;
            }

            return output;
        }

        /// <summary>
        /// Calculates the Euclidian (L2) norm of the vector x.
        /// </summary>
        /// <param name="x">The input vector.</param>
        /// <returns>The norm.</returns>
        public static double NormL2(
            [NotNull] double[] x)
        {
            double output = x.Sum(a => a * a);
            return output;
        }

        /// <summary>
        /// Calculates the product of a vector and a matrix, vector * matrix.
        /// </summary>
        /// <param name="vector">The input vector.</param>
        /// <param name="matrix">The matrix.</param>
        /// <param name="output">The output vector.</param>
        public static void Product(
            [NotNull] double[] vector,
            [NotNull] double[,] matrix,
            [NotNull] double[] output)
        {
            Chk.Equal(vector.Length, matrix.GetLength(0), "The length of the input vector must match the number of rows of the matrix.");
            Chk.Equal(output.Length, matrix.GetLength(1), "The length of the output vector must match the number of columns of the matrix.");

            for (int iCol = 0; iCol < vector.Length; iCol++)
            {
                double x = 0.0;
                
                for (int iRow = 0; iRow < output.Length; iRow++)
                {
                    x += matrix[iRow, iCol] * vector[iRow];
                }
                
                output[iCol] = x;
            }
        }

        /// <summary>
        /// Calculates the product of a vector and a matrix, vector * matrix.
        /// </summary>
        /// <param name="vector">The input vector.</param>
        /// <param name="matrix">The matrix.</param>
        /// <returns>The output vector.
        [return: NotNull]
        public static double[] Product(
            [NotNull] double[] vector,
            [NotNull] double[,] matrix)
        {
            double[] output = new double[matrix.GetLength(1)];

            Product(vector, matrix, output);

            return output;
        }
        
        /// <summary>
        /// Calculates the product of a matrix and a vector, matrix * vector.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="vector">The input vector.</param>
        /// <param name="output">The output vector.</param>
        public static void Product(
            [NotNull] double[,] matrix,
            [NotNull] double[] vector,
            [NotNull] double[] output)
        {
            Chk.Equal(vector.Length, matrix.GetLength(1), "The length of the input vector must match the number of columns of the matrix.");
            Chk.Equal(output.Length, matrix.GetLength(0), "The length of the output vector must match the number of rows of the matrix.");

            for (int iRow = 0; iRow < output.Length; iRow++)
            {
                double x = 0.0;
                for (int iCol = 0; iCol < vector.Length; iCol++)
                {
                    x += matrix[iRow, iCol] * vector[iCol];
                }
                
                output[iRow] = x;
            }
        }

        /// <summary>
        /// Calculates the product of a matrix and a vector, matrix * vector.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="vector">The input vector.</param>
        [return: NotNull]
        public static double[] Product(
            [NotNull] double[,] matrix,
            [NotNull] double[] vector)
        {
            double[] output = new double[matrix.GetLength(0)];

            Product(matrix, vector, output);

            return output;
        }
        
        /// <summary>
        /// Calculates the product of two matrices, A * B.
        /// </summary>
        /// <param name="a">Matrix A.</param>
        /// <param name="b">Matrix B.</param>
        /// <param name="output">The output.  The initial values can be anything.</param>
        public static void Product(
            [NotNull] double[,] a,
            [NotNull] double[,] b,
            [NotNull] double[,] output)
        {
            int aRows = a.GetLength(0);
            int aCols = a.GetLength(1);
            int bRows = b.GetLength(0);
            int bCols = b.GetLength(1);
            
            Chk.True(aCols == bRows, "The number of columns for matrix 'A' must match the number of rows for matrix 'B'.");
            Chk.Equal(aRows, output.GetLength(0), "The output must have the same number of rows as matrix 'A'.");
            Chk.Equal(bCols, output.GetLength(1), "The output must have the same number of columns as matrix 'B'.");

            for (int iRow = 0; iRow < aRows; iRow++)
            {
                for (int iCol = 0; iCol < bCols; iCol++)
                {
                    double x = 0;

                    for (int iSlide = 0; iSlide < aCols; iSlide++)
                    {
                        x += a[iRow, iSlide] * b[iSlide, iCol];
                    }

                    output[iRow, iCol] = x;
                }
            }
        }

        /// <summary>
        /// Calculates the product of two matrices, A * B.
        /// </summary>
        /// <param name="a">Matrix A.</param>
        /// <param name="b">Matrix B.</param>
        [return: NotNull]
        public static double[,] Product(
            [NotNull] double[,] a,
            [NotNull] double[,] b)
        {
            double[,] output = new double[a.GetLength(0), b.GetLength(1)];
            
            Product(a, b, output);
            
            return output;
        }
        
        /// <summary>
        /// Generates a random n * n rotation matrix.  This is useful for generating a random orthonormal basis.
        /// </summary>
        /// <param name="n">The size of the matrix.</param>
        /// <returns>The random rotation matrix.</returns>
        [return: NotNull]
        public static double[,] RandomRotationMatrix(int n)
        {
            int i, j, k;
            double theta, c, s, z;

            //    Initialize the identity matrix.
            double[,] output = new double[n, n];
            for (i = 0; i < n; i++)
                output[i, i] = 1.0;
            
            //    Special case, n==1.
            if (n == 1)
            {
                if (Morpe.Util.Rand.NextDouble() < 0.5)
                    output[0, 0] = -1.0;
                return output;
            }

            //    For each pair (i,j) of rows in the rotation matrix
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i != j)
                    {
                        //    Random rotation angle.
                        theta = Morpe.Util.Rand.NextDouble() * Pi_2; // planar rotation counterclockwise by theta
                        c = Math.Cos(theta);    // R(i,i) and  R(j,j)
                        s = Math.Sin(theta);    // R(j,i) and -R(i,j)

                        //    For each column R(:,i) and R(:,j), multiply each row M(i,:) and M(j,:).  In-place storage (in M) is possible.
                        for (k = 0; k < n; k++)
                        {
                            //    X(i,k)    =    R(i,i)*M(i,k)    +    R(i,j)*M(j,k)
                            z = c * output[i,k] - s * output[j,k];
                            //    X(j,k)    =    R(j,i)*M(i,k)    +    R(j,j)*M(j,k)
                            output[j,k] = s * output[i,k] + c * output[j,k];
                            //    In place storage.
                            output[i,k] = z;
                        }
                    }
                }
            }
            return output;
        }

        /// <summary>
        /// Subtract one vector from another.
        /// </summary>
        /// <param name="a">A vector.</param>
        /// <param name="b">A vector, same size as the other.</param>
        /// <returns>a - b ('a' minus 'b') as a vector, same size as the inputs.</returns>
        [return: NotNull]
        public static double[] Subtract(
            [NotNull] double[] a,
            [NotNull] double[] b)
        {
            Chk.Equal(a.Length, b.Length, "The vectors must be the same length.");

            double[] output = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                output[i] = a[i] - b[i];
            }

            return output;
        }

        /// <summary>
        /// Subtract one matrix from another.
        /// </summary>
        /// <param name="a">A matrix.</param>
        /// <param name="b">A matrix, same size as the other.</param>
        /// <returns>a - b ('a' minus 'b') as a matrix, same size as the inputs.</returns>
        [return: NotNull]
        public static double[,] Subtract(
            [NotNull] double[,] a,
            [NotNull] double[,] b)
        {
            int numRows = a.GetLength(0);
            int numCols = a.GetLength(1);
            Chk.Equal(numRows, b.GetLength(0), "The matrices must be the same size.");
            Chk.Equal(numCols, b.GetLength(1), "The matrices must be the same size.");
            
            double[,] output = new double[numRows, numCols];
            for (int i = 0; i < numRows; i++)
            {
                for (int j = 0; j < numCols; j++)
                {
                    output[i, j] = a[i, j] - b[i, j];
                }
            }

            return output;
        }
    }
}
