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

        [return: MaybeNull]
        public static float[] Convert([MaybeNull] double[] input)
        {
            if (input == null)
            {
                return null;
            }
            
            float[] output = new float[input.Length];
            for (int i = 0; i < input.Length; i++)
                output[i] = (float)input[i];

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
        /// Quicksorts an indexed list by changing the order of the indices and not the list itself.
        /// </summary>
        /// <param name="idx">Zero-based indices into the list x.  On output, these indices are reordered so that x[idx] is sorted.</param>
        /// <param name="x">The list to be sorted.  The elements are not changed on output.</param>
        /// <param name="left">A zero-based index identifying the lowest element of idx where sorting is conducted.</param>
        /// <param name="right">A zero-based index identifying the highest element of idx where sorting is conducted.</param>
        public static void QuickSortIndex(int[] idx, float[] x, int left, int right)
        {
            int i = left, j = right;
            int tmp;
            double pivot = x[idx[(left + right) / 2]];

            //    Partition
            while (i <= j)
            {
                while (x[idx[i]] < pivot) i++;
                while (x[idx[j]] > pivot) j--;
                if (i <= j)
                {
                    tmp = idx[i];
                    idx[i++] = idx[j];
                    idx[j--] = tmp;
                }
            }

            //    Recursion
            if (left < j)
                QuickSortIndex(idx, x, left, j);
            if (i < right)
                QuickSortIndex(idx, x, i, right);
        }
        
        /// <summary>
        /// Quicksorts an indexed table by an individual column by changing the order of the indices and not the table itself.
        /// </summary>
        /// <param name="idx">Zero-based indices into the table x.  On output, these indices are reordered so that x[idx][iCol] is sorted.</param>
        /// <param name="x">The table to be sorted.  The elements are not changed on output.</param>
        /// <param name="iCol">The zero-based column index which determines the rank order of rows within the table.</param>
        /// <param name="left">A zero-based index identifying the lowest row of idx where sorting is conducted.</param>
        /// <param name="right">A zero-based index identifying the highest row of idx where sorting is conducted.</param>
        public static void QuickSortIndex(int[] idx, float[][] x, int iCol, int left, int right)
        {
            int i = left, j = right;
            int tmp;
            double pivot = x[idx[(left + right) / 2]][iCol];

            //    Partition
            while (i <= j)
            {
                while (x[idx[i]][iCol] < pivot) i++;
                while (x[idx[j]][iCol] > pivot) j--;
                if (i <= j)
                {
                    tmp = idx[i];
                    idx[i++] = idx[j];
                    idx[j--] = tmp;
                }
            }

            //    Recursion
            if (left < j)
                QuickSortIndex(idx, x, iCol, left, j);
            if (i < right)
                QuickSortIndex(idx, x, iCol, i, right);
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
