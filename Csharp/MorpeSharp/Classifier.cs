using System;

namespace Morpe
{
	/// <summary>
	///	Represents an instance of the MoRPE classifier.
	/// </summary>
	public class Classifier
	{
		/// <summary>
		/// The number of categories.
		/// </summary>
		public readonly int Ncats;
		/// <summary>
		/// The spatial dimensionality.
		/// </summary>
		public readonly int Ndims;
		/// <summary>
		/// The number of polynomials.
		/// </summary>
		public readonly int Npoly;
		/// <summary>
		/// Allows a multivariate polynomial to be constructed.
		/// </summary>
		public readonly Poly Coeffs;
		/// <summary>
		/// The model parameters.  Each row specifies coefficients for a polynomial, so this array has <see cref="Npoly"/>
		/// rows and <see cref="Coeffs"/>.Ncoeff columns (See <see cref="Poly.Ncoeff"/>).
		/// </summary>
		public readonly float[][] Params;
		/// <summary>
		/// Holds data related to the quantization of the decision variable (for each polynomial) throughout the training sample.
		/// Creates a monotonic mapping from the decision variable to a probability of category membership;
		/// </summary>
		public readonly Quantization[] Quant;
		/// <summary>
		/// Initializes an untrained instance of the <see cref="Morpe.Classifier"/> class.  The classifier should be trained
		/// (see <see cref="Train"/>) before it can be used for classification.
		/// </summary>
		/// <param name="nCats"><see cref="Ncats"/></param>
		/// <param name="nDims"><see cref="Ndims"/></param>
		/// <param name="rank"><see cref="Rank"/></param>
		public Classifier(int nCats, int nDims, int rank)
		{
			this.Ncats = nCats;
			this.Ndims = nDims;
			this.Coeffs = new Poly (nDims, rank);
			this.Npoly = 1;
			if (this.Ncats > 2)
				this.Npoly = this.Ncats;
			this.Params = Static.NewArrays<float>(this.Npoly, this.Coeffs.Ncoeffs);
			this.Quant = new Quantization[this.Npoly];
		}
		protected Classifier(Classifier toCopy)
		{
			this.Ncats = toCopy.Ncats;
			this.Ndims = toCopy.Ndims;
			this.Coeffs = new Poly(this.Ndims, toCopy.Coeffs.Rank);
			this.Npoly = 1;
			if (this.Ncats > 2)
				this.Npoly = this.Ncats;
			this.Params = Static.NewArrays<float>(this.Npoly, this.Coeffs.Ncoeffs);
			Static.Copy(toCopy.Params, this.Params);
		}
		/// <summary>
		/// Classifies the multivariate coordinate.  (It should not be expanded.)
		/// </summary>
		/// <param name="x">The multivariate coordiante.</param>
		/// <returns>A vector of length <see cref="Ncats"/> giving the conditional probability of category membership for each category.
		/// Sums to exactly 1.0 (guaranteed).</returns>
		public double[] Classify(float[] x)
		{
			return this.ClassifyExpanded(this.Coeffs.Expand(x));
		}
		/// <summary>
		/// Classifies the expanded multivariate coordinate.
		/// </summary>
		/// <param name="x">The expanded multivariate coordinate.  For more information, see <see cref="Poly.Expand"/>.</param>
		/// <returns>A vector of length <see cref="Ncats"/> giving the conditional probability of category membership for each category.
		/// Sums to exactly 1.0 (guaranteed).</returns>
		public double[] ClassifyExpanded(float[] x)
		{
			double p = 0.0;
			if (this.Npoly == 1)
			{
				double y = this.EvalPolyFromExpanded(0, x);
				Quantization q = this.Quant[0];
				if (y < q.Ymid[0])
					p = q.P[0];
				else if (y > q.Ymid[q.Nquantiles - 1])
					p = q.P[q.Nquantiles - 1];
				else
					p = Static.Linterp(q.Ymid, q.P, y);
				return new double[] { p, 1.0 - p };
			}
			else
			{
				double[] output = new double[this.Ncats];
				double pSum = 0.0;
				for (int iCat = 0; iCat < this.Ncats; iCat++)
				{
					double y = this.EvalPolyFromExpanded(iCat, x);
					Quantization q = this.Quant[iCat];
					if (y < q.Ymid[0])
						p = q.P[0];
					else if (y > q.Ymid[q.Nquantiles - 1])
						p = q.P[q.Nquantiles - 1];
					else
						p = Static.Linterp(q.Ymid, q.P, y);
					output[iCat] = p;
					pSum += p;
				}
				for (int iCat = 0; iCat < this.Ncats; iCat++)
					output[iCat] /= pSum;
				return output;
			}
		}
		/// <summary>
		/// Creates a copy of the classifier which utilizes totally different memory resources.
		/// </summary>
		/// <returns>A copy of the classifier, with identical parameter values, but utilizing completely different memory resources.</returns>
		public Classifier Copy()
		{
			return new Classifier(this);
		}
		/// <summary>
		/// Evaluates the specified polynomial function for an expanded multivariate coordinate.
		/// </summary>
		/// <param name="iPoly">The polynomial to evaluate.  This is a zero-based index into a row of <see cref="Params"/>.</param>
		/// <param name="x">The expanded multivariate coordinate.  For more information, see <see cref="Poly.Expand"/>.</param>
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
		/// Trains the classifier by optimizing parameters based on the training data using default solver options.
		/// </summary>
		/// <param name="ops">Sets the <see cref="Trainer.Options"/> of the trainer.  If a value is not provided, default options are used.</param>
		public Trainer Train(CategorizedData data, SolverOptions ops)
		{
			if (data == null)
				throw new ArgumentException("Argument cannot be null.");

			Trainer trainer = new Trainer(this);
			trainer.Train(data, ops);
			return trainer;
		}
	}
}