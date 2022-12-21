using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Morpe.Validation;

namespace Morpe.Numerics.F
{
    public class Util
    {
        /// <summary>
        /// Adds the elements of "input" to the elements of "output".
        /// </summary>
        /// <param name="input">The matrix of values to be added.  This is not modified.</param>
        /// <param name="output">The matrix where summation is compiled.  This is modified.</param>
        public static void Add(
            [NotNull] float[][] input,
            [NotNull] float[][] output)
        {
            Chk.NotNull(input, nameof(input));
            Chk.NotNull(output, nameof(output));
            
            if (output.Length != input.Length)
                throw new ArgumentException("The arrays must be of equal size.");
            
            for (int iRow = 0; iRow < output.Length; iRow++)
            {
                float[] outRow = output[iRow];
                float[] inRow = input[iRow];
                
                if (inRow == null && outRow == null)
                {
                    continue;
                }
                else if (inRow == null || outRow == null)
                {
                    throw new ArgumentException($"A null row was detected in either the {nameof(input)} or {nameof(output)} but not both.");
                }
                else
                {
                    if (outRow.Length != inRow.Length)
                        throw new ArgumentException("The arrays must be of equal size.");
                    for (int iCol = 0; iCol < outRow.Length; iCol++)
                        outRow[iCol] += inRow[iCol];
                }
            }
        }
        
        /// <summary>
        /// Adds the elements of "input" to the elements of "output".
        /// </summary>
        /// <param name="inputScalar">The input scalar.</param>
        /// <param name="input">The matrix of values to be added.  This is not modified.</param>
        /// <param name="output">The matrix where summation is compiled.  This is modified.</param>
        public static void AddScaled(
            double inputScalar,
            [NotNull] float[][] input,
            [NotNull] float[][] output)
        {
            Chk.Finite(inputScalar, nameof(inputScalar));
            Chk.NotNull(input, nameof(input));
            Chk.NotNull(output, nameof(output));
            
            if (output.Length != input.Length)
                throw new ArgumentException("The arrays must be of equal size.");
            
            for (int iRow = 0; iRow < output.Length; iRow++)
            {
                float[] outRow = output[iRow];
                float[] inRow = input[iRow];

                if (inRow == null && outRow == null)
                {
                    continue;
                }
                else if (inRow == null || outRow == null)
                {
                    throw new ArgumentException($"A null row was detected in either the {nameof(input)} or {nameof(output)} but not both.");
                }
                else
                {
                    if (outRow.Length != inRow.Length)
                        throw new ArgumentException("The arrays must be of equal size.");
                    for (int iCol = 0; iCol < outRow.Length; iCol++)
                        outRow[iCol] += (float)(inRow[iCol] * inputScalar);
                }
            }
        }

        /// <summary>
        /// Convert a vector from double to float.
        /// </summary>
        /// <param name="input">The input values (double).</param>
        /// <returns>The output values (float).</returns>
        [return: NotNull]
        public static float[] From([NotNull] double[] input)
        {
            float[] output = new float[input.Length];

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = (float)input[i];
            }

            return output;
        }

        /// <summary>
        /// Allocates a jagged 3D array.
        /// </summary>
        /// <param name="depth">Length of first index.</param>
        /// <param name="height">Length of second index.</param>
        /// <param name="width">Length of third index.</param>
        /// <returns>The jagged 3D array.</returns>
        [return: NotNull]
        public static float[][][] JaggedCube(int depth, int height, int width)
        {
            float[][][] output = new float[depth][][];

            for (int i = 0; i < depth; i++)
            {
                output[i] = JaggedSquare(height, width);
            }

            return output;
        }
        
        /// <summary>
        /// Allocates a jagged 2D array.
        /// </summary>
        /// <param name="height">Length of first index.</param>
        /// <param name="width">Length of second index.</param>
        /// <returns>The jagged 2D array.</returns>
        [return: NotNull]
        public static float[][] JaggedSquare(int height, int width)
        {
            float[][] output = new float[height][];

            for (int i = 0; i < height; i++)
            {
                output[i] = new float[width];
            }

            return output;
        }

        /// <summary>
        /// Multiplies the elements of a matrix "output" by a scalar.
        /// </summary>
        /// <param name="output">The matrix.</param>
        /// <param name="scalar">The scalar</param>
        public static void Multiply(float[][] output, float scalar)
        {
            if (output == null)
                return;
            for (int iRow = 0; iRow < output.Length; iRow++)
            {
                float[] row = output[iRow];
                if (row != null)
                {
                    for (int iCol = 0; iCol < row.Length; iCol++)
                        row[iCol] *= scalar;
                }
            }
        }
        
        /// <summary>
        /// Calculates the Euclidian (L2) norm of the vector x.
        /// </summary>
        /// <param name="x">The input vector.</param>
        /// <returns>The norm.</returns>
        public static double NormL2(
            [NotNull] float[] x)
        {
            double output = Math.Sqrt(x.Sum(a => (double)a * a));
            return output;
        }
        
        /// <summary>
        /// Calculates the Euclidian (L2) norm of the elements of 'x'.  The input is treated as a single vector, despite
        /// being passed in as a 2D jagged array.
        /// </summary>
        /// <param name="x">The 2D jagged array which we treat like a single input vector.</param>
        /// <returns>The norm.</returns>
        public static double NormL2(
            [NotNull] float[][] x)
        {
            double output = 0.0;
            
            foreach (float[] xx in x)
            {
                if (xx == null)
                    continue;
                
                output += xx.Sum(a => (double)a * a);
            }
            output = Math.Sqrt(output);
            return output;
        }

        /// <summary>
        /// Scales each element of the given array by the given scalar. 
        /// </summary>
        /// <param name="array">The array.  The contents are modified by this method.</param>
        /// <param name="scalar">The scalar.</param>
        public static void Scale([MaybeNull] float[][] array, float scalar)
        {
            if (array != null)
            {
                for(int iRow=0; iRow<array.Length; iRow++)
                    for(int iCol=0; iCol<array[iRow].Length; iCol++)
                        array[iRow][iCol] *= scalar;
            }
        }
    }
}
