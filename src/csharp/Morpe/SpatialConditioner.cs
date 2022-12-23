using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Morpe.Validation;

namespace Morpe
{
    /// <summary>
    /// This object is used to transform data back and forth between (1) the original feature space and (2) the conditioned feature space.
    /// The conditioned feature space has a mean of 0 and a spread of 1 in each direction.  This is useful for minimizing the effect
    /// of errors such as floating point truncation.
    /// </summary>
    public class SpatialConditioner
    {
        /// <summary>
        /// Represents the origin of the feature space (specified in original units).  This variable has this.Ndims columns.
        /// </summary>
        public float[] Origin;
        
        /// <summary>
        /// Represents the spread of data in the feature space (specified in original units).  This measure is similar to a standard
        /// deviation, but it measures deviation from the origin instead of the mean.
        /// </summary>
        public float[] Spread;
        
        /// <summary>
        /// The spatial dimensionality.
        /// </summary>
        public int NumDims;
        
        /// <summary>
        /// Allocates memory for a new Spatial conditioner with the given dimensionality.
        /// </summary>
        /// <param name="numDims">The spatial dimensionality.</param>
        public SpatialConditioner(int numDims)
        {
            this.NumDims = numDims;
            this.Origin = new float[numDims];
            this.Spread = new float[numDims];
        }
        
        /// <summary>
        /// Create a deep copy of the instance.
        /// </summary>
        /// <returns>A deep copy.</returns>
        public SpatialConditioner Clone()
        {
            SpatialConditioner output = new SpatialConditioner(this.NumDims);
            Array.Copy(this.Origin, output.Origin, this.NumDims);
            Array.Copy(this.Spread, output.Spread, this.NumDims);
            return output;
        }
        
        /// <summary>
        /// Conditions the data.
        /// </summary>
        /// <param name="data">The data to be conditioned.</param>
        public void Condition([NotNull] CategorizedData data)
        {
            for (int iCat = 0; iCat < data.NumCats; iCat++)
            {
                for (int iRow = 0; iRow < data.NumEach[iCat]; iRow++)
                {
                    this.ConditionInPlace(data.X[iCat][iRow]);
                }
            }
        }
        
        /// <summary>
        /// Conditions a coordinate by subtracting the <see cref="Origin"/> and then dividing by the <see cref="Spread"/> for each dimension.
        /// In other words, this outputs something like a z-score for each dimension.
        /// </summary>
        /// <param name="x">The input coordinate.</param>
        /// <returns></returns>
        public float[] Condition([NotNull] IReadOnlyList<float> x)
        {
            float[] output = new float[x.Count];

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = (x[i] - this.Origin[i]) / this.Spread[i];
            }

            return output;
        }
        
        /// <summary>
        /// Conditions a coordinate by subtracting the <see cref="Origin"/> and then dividing by the <see cref="Spread"/> for each dimension.
        /// In other words, this outputs something like a z-score for each dimension.
        /// </summary>
        /// <param name="x">On input, the unconditioned coordinate.  On output, the conditioned coordinate.</param>
        public void ConditionInPlace([NotNull] float[] x)
        {
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = (x[i] - this.Origin[i]) / this.Spread[i];
            }
        }
        
        /// <summary>
        /// Deconditions the data.
        /// </summary>
        /// <param name="data">The data to be deconditioned.</param>
        public void Decondition([NotNull] CategorizedData data)
        {
            for (int iCat = 0; iCat < data.NumCats; iCat++)
            {
                for (int iRow = 0; iRow < data.NumEach[iCat]; iRow++)
                {
                    this.DeconditionInPlace(data.X[iCat][iRow]);
                }
            }
        }
        
        /// <summary>
        /// Deconditions a coordinate by multiplying by the <see cref="Spread"/> and then adds the <see cref="Origin"/> for each dimension.
        /// In other words, this takes something like a z-score (for each dimension) and outputs the original coordinate.
        /// </summary>
        /// <param name="x">The conditioned input coordinate.</param>
        /// <returns>The deconditioned (original) coordinate.</returns>
        public float[] Decondition([NotNull] IReadOnlyList<float> x)
        {
            float[] output = new float[x.Count];

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = x[i] * this.Spread[i] + this.Origin[i];
            }

            return output;
        }
        
        /// <summary>
        /// Deconditions a coordinate by multiplying by the <see cref="Spread"/> and then adds the <see cref="Origin"/> for each dimension.
        /// In other words, this takes something like a z-score (for each dimension) and outputs the original coordinate.
        /// </summary>
        /// <param name="x">On input, the conditioned coordinate.  On output, the deconditioned (original) coordinate.</param>
        public void DeconditionInPlace([NotNull] float[] x)
        {
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = x[i] * this.Spread[i] + this.Origin[i];
            }
        }

        /// <summary>
        /// Creates a spatial conditioner for a subset of the dimensions included in this one.
        /// </summary>
        /// <param name="dims">The spatial dimensions to be selected.  Zero-based indices.</param>
        /// <returns>The spatial conditioner for the given subspace.</returns>
        [return: NotNull]
        public SpatialConditioner ForSubspace(
            [NotNull] int[] dims)
        {
            SpatialConditioner output = new SpatialConditioner(dims.Length);

            for (int i = 0; i < dims.Length; i++)
            {
                int d = dims[i];
                Chk.LessOrEqual(0, d, "A negative spatial dimension, {0}, was given.", d);
                Chk.Less(d, this.NumDims, "A given dimension, {0}, was greater or equal to the limit of {1}.",
                    d,
                    this.NumDims);

                output.Origin[i] = this.Origin[d];
                output.Spread[i] = this.Spread[d];
            }

            return output;
        }
    }
}
