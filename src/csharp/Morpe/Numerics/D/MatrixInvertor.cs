using System;
using System.Diagnostics.CodeAnalysis;
using Morpe.Validation;

namespace Morpe.Numerics.D
{
    /// <summary>
    /// Encapsulates matrix inversion.  This can also calculate the determinant of a matrix.
    ///
    /// This uses LU decomposition and backsubstitution.
    ///
    /// This class is designed to minimize heap allocations.  The heap allocations are done once (during construction
    /// and the first operation), and they are never done again.  For this reason, the class does not support parallel
    /// execution on multiple threads, but the caller does not need to worry:  An internal mutex is used to prevent
    /// multi-threaded execution.
    /// </summary>
    public class MatrixInvertor
    {
        /// <summary>
        /// The rank of the matrices that can be inverted by this class.
        /// </summary>
        public int Rank { get; private set; }

        /// <summary>
        /// Constructs a new instance.  Allocates most heap resources.
        /// </summary>
        /// <param name="rank">The rank of the square matrix.  This is the number of rows (or columns).</param>
        public MatrixInvertor(int rank)
        {
            this.Rank = rank;
            this.lhsMatrix = Morpe.Util.NewArrays<double>(rank, rank);
            this.rhsVector = new double[rank];
            this.idxRows = new int[rank];
            this.mutex = this.lhsMatrix;
        }

        /// <summary>
        /// Calculates the determinant of the given matrix.
        /// </summary>
        /// <param name="matrix">The given matrix.</param>
        /// <returns>The determinant.</returns>
        public double Determinant(
            [NotNull] double[,] matrix)
        {
            Chk.NotNull(matrix, "The matrix cannot be null.");
            
            int len = matrix.GetLength(0);
            
            Chk.Equal(len, this.Rank, "The rank of the matrix is not compatible with this MatrixInvertor.");
            Chk.Equal(len, matrix.GetLength(1), "The matrix must be square.");

            double output = 1.0;
            
            lock (this.mutex)
            {
                for (int i = 0; i < len; i++)
                for (int j = 0;j < len;j++)
                    //	Initialize
                    this.lhsMatrix[i][j] = matrix[i,j];
            
                //	Inversion by LU decomposition and then back-substitution
                if(!this.LuDecomp(ref output))
                    return Double.NaN;

                for (int i = 0; i < len; i++)
                {
                    output *= this.lhsMatrix[i][i];
                }
            }

            return output;
        }
        
        /// <summary>
        /// Calculates the inverse of the given matrix.
        /// </summary>
        /// <param name="matrix">The given matrix.</param>
        /// <param name="output">Outputs the inverse of the given matrix.</param>
        /// <returns>True if the matrix could be inverted, false otherwise.  If false, the matrix is probably singular.
        /// </returns>
        public bool Invert(
            [NotNull] double[,] matrix,
            [NotNull] double[,] output)
        {
            Chk.NotNull(matrix, nameof(matrix));
            Chk.NotNull(output, nameof(output));
            
            int len = matrix.GetLength(0);
            
            Chk.Equal(len, this.Rank, "The rank of the matrix is not compatible with this MatrixInvertor.");
            Chk.Equal(len, matrix.GetLength(1), "The matrix must be square.");
            Chk.True(len == output.GetLength(0) && len == output.GetLength(1),
                "The output must be the same size as the input matrix.");
            
            double determinant = 1.0;

            lock (this.mutex)
            {
                int i, j;
                for (i = 0; i < len; i++)
                for (j = 0; j < len; j++)
                    //	Initialize
                    this.lhsMatrix[i][j] = matrix[i, j];

                //	Inversion by LU decompusition and then back-substitution
                if (!this.LuDecomp(ref determinant))
                    return false;
                else
                {
                    for (j = 0; j < len; j++)
                    {
                        for (i = 0; i < len; i++)
                            this.rhsVector[i] = 0.0;
                        this.rhsVector[j] = 1.0;
                        this.LuBacksub();
                        for (i = 0; i < len; i++)
                            output[i, j] = this.rhsVector[i];
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Calculates the inverse of the given matrix.
        /// </summary>
        /// <param name="matrix">The given matrix.</param>
        /// <returns>The inverse.  If null, the given matrix could not be inverted because it was probably singular.
        /// </returns>
        [return: MaybeNull]
        public double[,] Invert(
            [NotNull] double[,] matrix)
        {
            Chk.NotNull(matrix, nameof(matrix));

            double[,] output = new double[matrix.GetLength(0), matrix.GetLength(1)];
            if (!this.Invert(matrix, output))
            {
                output = null;
            }
            return output;
        }

        /// <summary>
        /// An arbitrary tiny value.
        /// </summary>
        private static readonly double Tiny = 1e-36;
        
        /// <summary>
        /// In regards to the linear system Ax = b, this is the matrix 'A' which initially stores the original matrix
        /// prior to LU decomposition (<see cref="LuDecomp"/>), and after it stores the LU decomposed matrix.
        /// </summary>
        private double[][] lhsMatrix;
        
        /// <summary>
        /// In regards to the linear system Ax = b, this is the vector 'b' used for LU backsubstitution
        /// (<see cref="LuBacksub"/>).
        /// </summary>
        private double[] rhsVector;
        
        /// <summary>
        /// Holds the indexes of the permuted rows of <see cref="lhsMatrix"/>.  This is set during
        /// <see cref="LuDecomp"/>.
        /// </summary>
        private int[] idxRows;
        
        /// <summary>
        /// An arbitrary object used to prevent multiple operations from occurring in parallel.
        /// </summary>
        private object mutex;
        
        /// <summary>
        /// A vector used as working space during LU decomposition.
        /// </summary>
        private double[] workingVec;

        /// <summary>
        /// Performs the LU decomposition.
        /// </summary>
        /// <param name="d">This is an output value that equals + or - 1.0 depending on whether or the number of row
        /// interchanges was even or odd respectively.</param>
        /// <returns>False if a is singular, true otherwise.</returns>
        private bool LuDecomp(ref double d)
        {
            int i, j, k;
            double big, dum, sum, temp;
            int imax = 0;
            
            if (this.workingVec == null || this.workingVec.Length < this.Rank)
                this.workingVec = new double[this.Rank];

            d = 1.0;
            for (i = 0;i < this.Rank;i++)
            {
                big = 0.0;
                for (j = 0;j < this.Rank;j++)
                    if ((temp = Math.Abs(this.lhsMatrix[i][j])) > big)
                        big = temp;
                 
                if (big == 0.0)
                    return false;
             
                this.workingVec[i] = 1.0 / big;
            }
            for (j = 0;j < this.Rank;j++)
            {
                for (i = 0;i < j;i++)
                {
                    sum = this.lhsMatrix[i][j];
                    for (k = 0;k < i;k++)
                        sum -= this.lhsMatrix[i][k] * this.lhsMatrix[k][j];
                    this.lhsMatrix[i][j] = sum;
                }
                big = 0.0;
                for (i = j;i < this.Rank;i++)
                {
                    sum = this.lhsMatrix[i][j];
                    for (k = 0;k < j;k++)
                        sum -= this.lhsMatrix[i][k] * this.lhsMatrix[k][j];
                    this.lhsMatrix[i][j] = sum;
                    if ((dum = this.workingVec[i] * Math.Abs(sum)) >= big)
                    {
                        big = dum;
                        imax = i;
                    }
                 
                }
                if (j != imax)
                {
                    for (k = 0;k < this.Rank;k++)
                    {
                        dum = this.lhsMatrix[imax][k];
                        this.lhsMatrix[imax][k] = this.lhsMatrix[j][k];
                        this.lhsMatrix[j][k] = dum;
                    }
                    d = -d;
                    this.workingVec[imax] = this.workingVec[j];
                }
             
                this.idxRows[j] = imax;
                if (this.lhsMatrix[j][j] == 0.0)
                    this.lhsMatrix[j][j] = Tiny;
             
                if (j != this.Rank - 1)
                {
                    dum = 1.0 / (this.lhsMatrix[j][j]);
                    for (i = j + 1;i < this.Rank;i++)
                        this.lhsMatrix[i][j] *= dum;
                }
            }
            return true;
        }

        /// <summary>
        /// Uses LU Back Substituion to solve the set of <see cref="Rank"/> linear equations (Ax=b).
        ///
        /// Several member v
        /// <see cref="lhsMatrix"/>  Inputs the LU Decomposition of the matrix 'A', as determined by <see cref="LuDecomp"/>
        /// this.indx Inputs the index into the row permutations of the LU Decomposition.
        /// <see cref="rhsVector"/> Inputs the variable 'b'.  Outputs the value of 'x'.
        /// </summary>
        private void LuBacksub()
        {
            int i, ii = 0, j;
            double sum;
            for (i = 0;i < this.Rank; i++)
            {
                int ip = this.idxRows[i];
                sum = this.rhsVector[ip];
                this.rhsVector[ip] = this.rhsVector[i];
                if (ii != 0)
                    for (j = ii - 1;j < i;j++)
                        sum -= this.lhsMatrix[i][j] * this.rhsVector[j];
                else if (sum != 0.0)
                    ii = i + 1;
              
                this.rhsVector[i] = sum;
            }
            for (i = this.Rank - 1;i >= 0;i--)
            {
                sum = this.rhsVector[i];
                for (j = i + 1; j < this.Rank; j++)
                    sum -= this.lhsMatrix[i][j] * this.rhsVector[j];
                this.rhsVector[i] = sum / this.lhsMatrix[i][i];
            }
        }
    }
}