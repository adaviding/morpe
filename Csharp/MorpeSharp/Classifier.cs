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
		/// The coefficients of a multivariate polynomial.
		/// </summary>
		public readonly Poly Coeffs;
		/// <summary>
		/// The solver parameters.  Each row specifies parameters for a polynomial, so this array has <see cref="Npoly"/> rows and Coeffs.<see cref="Poly.Ncoeff"/> columns.
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
		/// Creates a copy of the classifier which utilizes totally different memory resources.
		/// </summary>
		/// <returns>A copy of the classifier, with identical parameter values, but utilizing completely different memory resources.</returns>
		public Classifier Copy()
		{
			return new Classifier(this);
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