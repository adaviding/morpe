using System;
using System.Linq;

using F1 = Morpe.Numerics.F1;

namespace Morpe
{
	/// <summary>
	/// Multivariate data that is segregated by category membership.
	/// </summary>
	public class CategorizedData
	{
		/// <summary>
		/// The number of categories, also equal to <see cref="NumEach.Length"/>.
		/// </summary>
		public readonly int NumCats;
		
		/// <summary>
		/// The spatial dimensionality, and the number of rows of <see cref="X"/>.
		/// </summary>
		public readonly int NumDims;
		
		/// <summary>
		/// The number of data points in each category.  The number of categories is equal to <see cref="NumEach.Length"/>;
		/// </summary>
		public readonly int[] NumEach;
		
		/// <summary>
		/// The total number of data points (across all categories).
		/// </summary>
		public readonly int NumTotal;
		
		/// <summary>
		/// The data, indexed as [c][i][j].  Each page c is a category.  Each row i is a unique data point.  Each column
		/// j is a spatial axis.
		/// </summary>
		public readonly float[][][] X;
		
		/// <summary>
		/// Initializes a new instance of the <see cref="CategorizedData"/> class.  Allocates memory for pages, rows,
		/// and columns of <see cref="X"/>.
		/// </summary>
		/// <param name="numEach"><see cref="NumEach"/></param>
		/// <param name="numDims"><see cref="NumDims"/></param>
		public CategorizedData(int[] numEach, int numDims)
		{
			this.NumEach = numEach;
			this.NumDims = numDims;
			this.NumCats = this.NumEach.Length;
			this.NumTotal = this.NumEach.Sum();
			this.X = new float[this.NumCats][][];
			for (int iCat=0; iCat<this.NumCats; iCat++)
			{
				this.X[iCat] = new float[this.NumEach[iCat]][];
				int n = this.NumEach[iCat];
				for(int i=0; i<n; i++)
					this.X[iCat][i] = new float[this.NumDims];
			}
		}
		
		/// <summary>
		/// Gets the categorized data to train a "dual" classifier.  Dual classifiers are only used when the number
		/// of categories is greater than 2.  Each dual classifier is used to discriminate between the target
		/// category and the sum of all other categories.
		/// </summary>
		/// <param name="dualViewable">The original categorized data which should have more than 2 categories.</param>
		/// <param name="targetCat">The target category.</param>
		protected CategorizedData(CategorizedData dualViewable, int targetCat)
		{
			this.NumEach = new int[] { dualViewable.NumEach[targetCat], dualViewable.NumTotal - dualViewable.NumEach[targetCat] };
			this.NumDims = dualViewable.NumDims;
			this.NumCats = this.NumEach.Length;
			this.NumTotal = this.NumEach.Sum();

			//	Allocate page holders
			this.X = new float[this.NumCats][][];

			//	Allocate row holders
			for (int iCat=0; iCat<this.NumCats; iCat++)
				this.X[iCat] = new float[this.NumEach[iCat]][];

			//	Fill rows of category 0 with targetCat
			int n = this.NumEach[0];
			for (int i = 0; i < n; i++)
				this.X[0][i] = dualViewable.X[targetCat][i];

			//	Fill rows of category 1 with remaining data.
			int iDatum = 0;
			for (int iCat = 0; iCat < dualViewable.NumCats; iCat++)
			{
				if (iCat != targetCat)
				{
					n = dualViewable.NumEach[iCat];
					for (int iSamp = 0; iSamp < n; iSamp++)
						this.X[1][iDatum++] = dualViewable.X[iCat][iSamp];
				}
			}
		}
		
		/// <summary>
		/// This reverses the effects of <see cref="Expand"/>.  Expanded terms are removed if they exist.
		/// </summary>
		public void Contract()
		{
			for (int iCat = 0; iCat < this.NumCats; iCat++)
			{
				int nSamp = this.NumEach[iCat];
				for (int iSamp = 0; iSamp < nSamp; iSamp++)
				{
					float[] x = this.X[iCat][iSamp];
					if (x.Length > this.NumDims)
						this.X[iCat][iSamp] = F1.Util.GetSubarray(x, 0, this.NumDims - 1);
				}
			}
		}
		
		/// <summary>
		/// Computes the polynomial expansion in place.  You can undo the expansion by calling <see cref="Contract"/>
		/// </summary>
		/// <param name="poly">The defiinition of the polynomial expansion.</param>
		public void Expand(Poly poly)
		{
			for (int iCat = 0; iCat < this.NumCats; iCat++)
			{
				int nSamp = this.NumEach[iCat];
				for (int iSamp = 0; iSamp < nSamp; iSamp++)
				{
					float[] x = this.X[iCat][iSamp];
					if (x.Length > this.NumDims)
						this.X[iCat][iSamp] = poly.Expand(F1.Util.GetSubarray(x, 0, this.NumDims - 1));
					else
						this.X[iCat][iSamp] = poly.Expand(x);
				}
			}
		}
		
		/// <summary>
		/// If the data consists of more than 2 categories, this function returns a "dual" view of the classification data.  This
		/// views the data as if it consisted only two categories.  In the dual view, the targetCat is represented as category 0 while
		/// all other data are represented as category 1.
		/// </summary>
		/// <param name="targetCat">The target category to be represented as category 0 in the dual view.</param>
		/// <returns>The dual view of this data set.</returns>
		public CategorizedData GetDual(int targetCat)
		{
			if (this.NumCats <= 2)
				return null;

			return new CategorizedData(this, targetCat);
		}
		
		/// <summary>
		/// If the data consists of more than 2 cateogries, this function returns all "dual" views of the classification data.
		/// See <see cref="GetDual"/> for more information.
		/// </summary>
		/// <returns>All dual views of the data.</returns>
		public CategorizedData[] GetDuals()
		{
			if (this.NumCats <= 2)
				return null;

			CategorizedData[] output = new CategorizedData[this.NumCats];
			for (int iCat = 0; iCat < this.NumCats; iCat++)
				output[iCat] = this.GetDual(iCat);
			return output;
		}
	}
}