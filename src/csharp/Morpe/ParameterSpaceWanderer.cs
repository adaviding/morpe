using System.Diagnostics.CodeAnalysis;
using Morpe.Validation;

using D = Morpe.Numerics.D;
using D1 = Morpe.Numerics.D1;
using F = Morpe.Numerics.F;

namespace Morpe
{
    /// <summary>
    /// This is used during optimization to take "random" steps in the parameter space.
    /// </summary>
    public class ParameterSpaceWanderer
    {
        /// <summary>
        /// If <see cref="NumParams"/> is less or equal to this value, then random rotation matrices are used as an
        /// orthonormal basis for the parameter search; otherwise, we just generate random Gaussian deviates.
        /// </summary>
        public static readonly int UpperLimitForActualOrthonormalBasis = 50;
        
        /// <summary>
        /// The number of free parameters.  This is equal to the number of polynomials times the number of parameters
        /// per polynomial.  It is also the number of floating point values returned by <see cref="NextBasis"/>;
        /// </summary>
        public readonly int NumParams;

        /// <summary>
        /// The number of free parameters (a.k.a. coefficients) per polynomial.
        /// </summary>
        public int NumParamsPerPoly => paramScale.Length;
        
        /// <summary>
        /// The number of polynomials for this classification problem.
        /// </summary>
        public readonly int NumPolys;
        
        /// <summary>
        /// Constructs a new instance.
        /// </summary>
        /// <param name="numPolys">The number of polynomial functions.</param>
        /// <param name="paramScale">The scale of parameters (coefficients) for each polynomial function.  This is
        /// typically just <see cref="PreTrainingAnalysis.ParamScale"/>.</param>
        public ParameterSpaceWanderer(int numPolys, [NotNull] float[] paramScale)
        {
            Chk.Less(0, numPolys, "The number of polynomials must be positive and non-zero.");
            Chk.NotNull(paramScale, nameof(paramScale));

            this.NumPolys = numPolys;
            this.paramScale = paramScale;
            this.NumParams = numPolys * this.NumParamsPerPoly;

            if (this.NumParams < UpperLimitForActualOrthonormalBasis)
            {
                // Initialize the basis.
                this.initBasis();
            }
        }

        /// <summary>
        /// Fetch the next basis.
        /// </summary>
        /// <param name="scale">The returned basis is scaled by this value.  This will be the L2 norm of the returned
        /// values.</param>
        /// <returns>The next basis (scaled).</returns>
        public float[][] NextBasis(float scale)
        {
            if (this.basis == null)
            {
                float[][] output = Util.NewArrays<float>(this.NumPolys, this.NumParamsPerPoly);

                if (this.NumParams == 1)
                {
                    output[0][0] = Util.Rand.NextDouble() > 0.5
                        ? scale
                        : -scale;

                    return output;
                }

                for (int iPoly = 0; iPoly < this.NumPolys; iPoly++)
                    for (int iCoeff = 0; iCoeff < this.NumParamsPerPoly; iCoeff++)
                        output[iPoly][iCoeff] = (float)D1.GaussianDistribution.Rand();

                double norm = F.Util.NormL2(output);
                F.Util.Scale(output, (float)(scale/norm));
                return output;
            }
            else
            {
                if (this.idxNextBasis >= this.basis.Length)
                    this.initBasis();

                float[][] output = this.basis[this.idxNextBasis++];
                F.Util.Scale(output, scale);
                return output;
            }
        }

        /// <summary>
        /// This is the random orthonormal basis which is used when the number of free parameters is less than or equal
        /// to <see cref="UpperLimitForActualOrthonormalBasis"/>.
        ///
        /// The size of this cube is:  [num params][num polys][num params per poly]
        /// </summary>
        private float[][][] basis;
        
        /// <summary>
        /// The index of the next basis vector to be retrieved.
        /// </summary>
        private int idxNextBasis = 0;
        
        /// <summary>
        /// The scale of parameters (coefficients) for each polynomial function.  This is
        /// typically just <see cref="PreTrainingAnalysis.ParamScale"/>.
        /// </summary>
        private readonly float[] paramScale;

        /// <summary>
        /// This initializes the search basis by using a new random rotation matrix.  It is only used when
        /// <see cref="NumParams"/> is less or equal to <see cref="UpperLimitForActualOrthonormalBasis"/>.
        /// </summary>
        private void initBasis()
        {
            this.idxNextBasis = 0;

            this.basis = Util.NewArrays<float>(this.NumParams, this.NumPolys, this.NumParamsPerPoly);
            double[,] ortho = D.Util.RandomRotationMatrix(this.NumParams);
            
            for (int iBasis = 0; iBasis < this.NumParams; iBasis++)
            {
                int jBasis = 0;
                
                for (int iPoly = 0;iPoly < this.NumPolys; iPoly++)
                {
                    for (int iCoeff = 0;iCoeff < this.paramScale.Length; iCoeff++)
                    {
                        this.basis[iBasis][iPoly][iCoeff] = (float)(ortho[iBasis,jBasis] * this.paramScale[iCoeff]);
                        jBasis++;
                    }
                }

                float scalar = (float)(1.0/F.Util.NormL2(this.basis[iBasis]));
                F.Util.Scale(this.basis[iBasis], scalar);
            }
        }
    }
}
