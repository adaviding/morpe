using System;
using System.Diagnostics.CodeAnalysis;
using System.Threading.Tasks;
using Morpe.Validation;
using D1 = Morpe.Numerics.D1;

namespace Morpe
{
    /// <summary>
    ///    Represents an instance of the MoRPE classifier.
    /// </summary>
    public class Classifier
    {
        /// <summary>
        /// Allows a multivariate polynomial to be constructed.
        /// </summary>
        public Polynomial Coeffs;

        /// <summary>
        /// This is used to condition the data prior to a polynomial expansion.
        /// </summary>
        public SpatialConditioner Conditioner;

        /// <summary>
        /// The number of categories.
        /// </summary>
        public int NumCats;

        /// <summary>
        /// The spatial dimensionality.
        /// </summary>
        public int NumDims;

        /// <summary>
        /// The number of polynomials.
        ///
        /// For the 2-category problem, this is equal to 1; otherwise it is <see cref="NumCats"/>.
        /// </summary>
        public int NumPoly;

        /// <summary>
        /// The number of quantiles used to "bin" probability for each polynomial.
        /// </summary>
        public int NumQuantiles;

        /// <summary>
        /// The model parameters.  Each row specifies coefficients for a polynomial, so this array has <see cref="NumPoly"/>
        /// rows and <see cref="Coeffs"/>.Ncoeff columns (See <see cref="Polynomial.NumCoeff"/>).
        /// </summary>
        public float[][] Params;

        /// <summary>
        /// Holds data related to the quantization of the decision variable (for each polynomial) throughout the training sample.
        /// Creates a monotonic mapping from the decision variable to a probability of category membership;
        /// </summary>
        public Quantization[] Quant;

        /// <summary>
        /// Initializes an untrained instance of the <see cref="Classifier"/> class.  The classifier should be trained
        /// (see <see cref="Train"/>) before it can be used for classification.
        /// </summary>
        /// <param name="numCats"><see cref="NumCats"/></param>
        /// <param name="numDims"><see cref="NumDims"/></param>
        /// <param name="rank"><see cref="Rank"/></param>
        /// <param name="numQuantiles"><see cref="NumQuantiles"/></param>
        public Classifier(int numCats, int numDims, int rank, int numQuantiles)
        {
            this.NumCats = numCats;
            this.NumDims = numDims;
            this.Coeffs = new Polynomial (numDims, rank);
            this.NumQuantiles = numQuantiles;
            this.NumPoly = 1;
            if (this.NumCats > 2)
                this.NumPoly = this.NumCats;
            this.Params = Util.NewArrays<float>(this.NumPoly, this.Coeffs.NumCoeffs);
            this.Quant = new Quantization[this.NumPoly];

            for (int i = 0; i < this.NumPoly; i++)
            {
                this.Quant[i] = new Quantization(numQuantiles);
            }
        }

        public Classifier(
            int numCats,
            int numDims,
            int rank,
            int numQuantiles,
            [NotNull] D1.Range probabilityRange,
            [NotNull] float[][] parameters)
        {
            this.NumCats = numCats;
            this.NumDims = numDims;
            this.Coeffs = new Polynomial (numDims, rank);
            this.NumQuantiles = numQuantiles;
            this.NumPoly = 1;
            if (this.NumCats > 2)
                this.NumPoly = this.NumCats;
            this.Params = parameters;
            this.Quant = new Quantization[this.NumPoly];

            Chk.NotNull(parameters, nameof(parameters));
            Chk.Equal(parameters.Length, this.NumPoly, "The number of parameter arrays {0} must be equal to the number of polynomials {1}.",
                parameters.Length,
                this.NumPoly);

            for (int i = 0; i < this.NumPoly; i++)
            {
                this.Quant[i] = new Quantization(numQuantiles, probabilityRange);

                Chk.NotNull(parameters[i], "{0}[{1}]", nameof(parameters), nameof(i));
                Chk.Equal(parameters[i].Length, this.Coeffs.NumCoeffs, "The parameters for polynomial {0} had an incorrect length {1}, expected {2}.",
                    i,
                    this.Params[i].Length,
                    this.Coeffs.NumCoeffs);
            }
        }

        /// <summary>
        /// This is used by <see cref="GetDual"/> and <see cref="GetDuals"/>.
        /// </summary>
        /// <param name="toCopy">A classifier having more than 1 polynomial.</param>
        /// <param name="targetPoly">The target polynomial function.</param>
        private Classifier(Classifier toCopy, int targetPoly)
        {
            this.NumCats = 2;
            this.NumDims = toCopy.NumDims;
            this.Coeffs = new Polynomial(this.NumDims, toCopy.Coeffs.Rank);
            this.NumPoly = 1;

            this.Params = Util.NewArrays<float>(this.NumPoly, this.Coeffs.NumCoeffs);
            Array.Copy(toCopy.Params[targetPoly], this.Params[0], this.Coeffs.NumCoeffs);

            this.Quant = new Quantization[this.NumPoly];
            Quantization q;
            q = toCopy.Quant[targetPoly];
            if (q != null)
                this.Quant[0] = q.Clone();
        }

        /// <summary>
        /// Classifies the multivariate coordinate.
        /// </summary>
        /// <param name="x">The multivariate coordiante.</param>
        /// <returns>A vector of length <see cref="NumCats"/> giving the conditional probability of category membership for each category.
        /// Sums to exactly 1.0 (guaranteed).</returns>
        public double[] Classify([NotNull] float[] x)
        {
            Chk.NotNull(this.Conditioner, nameof(this.Conditioner));
            Chk.NotNull(this.Coeffs, nameof(this.Coeffs));

            float[] conditioned = this.Conditioner.Condition(x);
            float[] expanded = this.Coeffs.Expand(conditioned);

            return this.ClassifyExpanded(expanded);
        }

        /// <summary>
        /// Classifies the expanded multivariate coordinate.
        /// </summary>
        /// <param name="x">The expanded multivariate coordinate.  For more information, see <see cref="Polynomial.Expand"/>.</param>
        /// <returns>A vector of length <see cref="NumCats"/> giving the conditional probability of category membership for each category.
        /// Sums to exactly 1.0 (guaranteed).</returns>
        public double[] ClassifyExpanded(float[] x)
        {
            double p = 0.0;
            if (this.NumPoly == 1)
            {
                double y = this.EvalPolyFromExpanded(0, x);
                Quantization q = this.Quant[0];
                if (y < q.Ymid[0])
                    p = q.P[0];
                else if (y > q.Ymid[q.NumQuantiles - 1])
                    p = q.P[q.NumQuantiles - 1];
                else
                    p = D1.Util.Linterp(q.Ymid, q.P, y);
                return new double[] { p, 1.0 - p };
            }
            else
            {
                double[] output = new double[this.NumCats];
                double pSum = 0.0;
                for (int iCat = 0; iCat < this.NumCats; iCat++)
                {
                    double y = this.EvalPolyFromExpanded(iCat, x);
                    Quantization q = this.Quant[iCat];
                    if (y < q.Ymid[0])
                        p = q.P[0];
                    else if (y > q.Ymid[q.NumQuantiles - 1])
                        p = q.P[q.NumQuantiles - 1];
                    else
                        p = D1.Util.Linterp(q.Ymid, q.P, y);
                    output[iCat] = p;
                    pSum += p;
                }
                for (int iCat = 0; iCat < this.NumCats; iCat++)
                    output[iCat] /= pSum;
                return output;
            }
        }

        /// <summary>
        /// Classifies the polynomial outputs.
        /// </summary>
        /// <param name="y">The output of each polynomial.</param>
        /// <returns>A vector of length <see cref="NumCats"/> giving the conditional probability of category membership for each category.
        /// Sums to exactly 1.0 (guaranteed).</returns>
        public double[] ClassifyPolynomialOutputs(float[] y)
        {
            double p = 0.0;
            if (this.NumPoly == 1)
            {
                Quantization q = this.Quant[0];
                if (y[0] < q.Ymid[0])
                    p = q.P[0];
                else if (y[0] > q.Ymid[q.NumQuantiles - 1])
                    p = q.P[q.NumQuantiles - 1];
                else
                    p = D1.Util.Linterp(q.Ymid, q.P, y[0]);
                return new double[] { p, 1.0 - p };
            }
            else
            {
                double[] output = new double[this.NumCats];
                double pSum = 0.0;
                for (int iCat = 0; iCat < this.NumCats; iCat++)
                {
                    Quantization q = this.Quant[iCat];
                    if (y[iCat] < q.Ymid[0])
                        p = q.P[0];
                    else if (y[iCat] > q.Ymid[q.NumQuantiles - 1])
                        p = q.P[q.NumQuantiles - 1];
                    else
                        p = D1.Util.Linterp(q.Ymid, q.P, y[iCat]);
                    output[iCat] = p;
                    pSum += p;
                }
                for (int iCat = 0; iCat < this.NumCats; iCat++)
                    output[iCat] /= pSum;
                return output;
            }
        }

        /// <summary>
        /// Creates a deep copy.
        /// </summary>
        /// <returns>The deep copy.</returns>
        [return: NotNull]
        public Classifier Clone()
        {
            Classifier output = (Classifier)this.MemberwiseClone();
            output.Coeffs = this.Coeffs?.Clone();
            output.Conditioner = this.Conditioner?.Clone();
            output.Params = Util.Clone(this.Params);
            output.Quant = Util.Clone(this.Quant);
            return output;
        }

        /// <summary>
        /// Evaluates the specified polynomial function for an expanded multivariate coordinate.
        /// </summary>
        /// <param name="iPoly">The polynomial to evaluate.  This is a zero-based index into a row of <see cref="Params"/>.</param>
        /// <param name="x">The expanded multivariate coordinate.  For more information, see <see cref="Polynomial.Expand"/>.</param>
        /// <returns>The value of the polynomial expression.</returns>
        public double EvalPolyFromExpanded(int iPoly, float[] x)
        {
            double output = 0.0;
            float[] poly = this.Params[iPoly];
            for (int i = 0; i < poly.Length; i++)
                output += poly[i] * x[i];
            return output;
        }

        /// <summary>
        /// If this classifier has more than 1 polynomial, this function returns the "dual" classifier consisting of 1 polynomial.
        /// It is named "dual" because it deals only with 2 categories.
        /// </summary>
        /// <param name="targetPoly">The target polynomial, also the target category.  All other polynomials are discarded, and all
        /// other categories are merged into 1.</param>
        /// <returns>The dual classifier.</returns>
        public Classifier GetDual(int targetPoly)
        {
            if (this.NumPoly == 1)
                return null;
            return new Classifier(this, targetPoly);
        }

        /// <summary>
        /// If this classifier has more than 1 polynomial, this function returns all "dual" classifiers.  See <see cref="GetDual"/> for more information.
        /// </summary>
        /// <returns>All dual classifiers.</returns>
        public Classifier[] GetDuals()
        {
            if (this.NumPoly == 1)
                return null;
            Classifier[] output = new Classifier[this.NumPoly];
            for (int iPoly = 0; iPoly < this.NumPoly; iPoly++)
                output[iPoly] = this.GetDual(iPoly);
            return output;
        }
    }
}
