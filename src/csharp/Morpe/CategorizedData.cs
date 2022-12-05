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
		/// The data, indexed as [c][i][j].  Each page c is a category.  Each row i is a unique data point.  Each column j is a spatial axis.
		/// </summary>
		public readonly float[][][] X;
		
		/// <summary>
		/// The total number of data points (across all categories).
		/// </summary>
		public readonly int Ntotal;
		
		/// <summary>
		/// The number of data points in each category.  The number of categories is equal to <see cref="Neach.Length"/>;
		/// </summary>
		public readonly int[] Neach;
		
		/// <summary>
		/// The spatial dimensionality, and the number of rows of <see cref="X"/>.
		/// </summary>
		public readonly int Ndims;
		
		/// <summary>
		/// The number of categories, also equal to <see cref="Neach.Length"/>.
		/// </summary>
		public readonly int Ncats;
		
		/// <summary>
		/// Initializes a new instance of the <see cref="CategorizedData"/> class.  Allocates memory for pages, rows, and columns of <see cref="X"/>.
		/// </summary>
		/// <param name="Neach"><see cref="Neach"/></param>
		/// <param name="Ndims"><see cref="Ndims"/></param>
		public CategorizedData(int[] Neach, int Ndims)
		{
			this.Neach = Neach;
			this.Ndims = Ndims;
			this.Ncats = this.Neach.Length;
			this.Ntotal = this.Neach.Sum();
			this.X = new float[this.Ncats][][];
			for (int iCat=0; iCat<this.Ncats; iCat++)
			{
				this.X[iCat] = new float[this.Neach[iCat]][];
				int n = this.Neach[iCat];
				for(int i=0; i<n; i++)
					this.X[iCat][i] = new float[this.Ndims];
			}
		}
		
		protected CategorizedData(CategorizedData dualViewable, int targetCat)
		{
			this.Neach = new int[] { dualViewable.Neach[targetCat], dualViewable.Ntotal - dualViewable.Neach[targetCat] };
			this.Ndims = dualViewable.Ndims;
			this.Ncats = this.Neach.Length;
			this.Ntotal = this.Neach.Sum();

			//	Allocate page holders
			this.X = new float[this.Ncats][][];

			//	Allocate row holders
			for (int iCat=0; iCat<this.Ncats; iCat++)
				this.X[iCat] = new float[this.Neach[iCat]][];

			//	Fill rows of category 0 with targetCat
			int n = this.Neach[0];
			for (int i = 0; i < n; i++)
				this.X[0][i] = dualViewable.X[targetCat][i];

			//	Fill rows of category 1 with remaining data.
			int iDatum = 0;
			for (int iCat = 0; iCat < dualViewable.Ncats; iCat++)
			{
				if (iCat != targetCat)
				{
					n = dualViewable.Neach[iCat];
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
			for (int iCat = 0; iCat < this.Ncats; iCat++)
			{
				int nSamp = this.Neach[iCat];
				for (int iSamp = 0; iSamp < nSamp; iSamp++)
				{
					float[] x = this.X[iCat][iSamp];
					if (x.Length > this.Ndims)
						this.X[iCat][iSamp] = F1.Util.GetSubarray(x, 0, this.Ndims - 1);
				}
			}
		}
		
		/// <summary>
		/// Computes the polynomial expansion in place.  You can undo the expansion by calling <see cref="Contract"/>
		/// </summary>
		/// <param name="poly">The defiinition of the polynomial expansion.</param>
		public void Expand(Poly poly)
		{
			for (int iCat = 0; iCat < this.Ncats; iCat++)
			{
				int nSamp = this.Neach[iCat];
				for (int iSamp = 0; iSamp < nSamp; iSamp++)
				{
					float[] x = this.X[iCat][iSamp];
					if (x.Length > this.Ndims)
						this.X[iCat][iSamp] = poly.Expand(F1.Util.GetSubarray(x, 0, this.Ndims - 1));
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
			if (this.Ncats <= 2)
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
			if (this.Ncats <= 2)
				return null;

			CategorizedData[] output = new CategorizedData[this.Ncats];
			for (int iCat = 0; iCat < this.Ncats; iCat++)
				output[iCat] = this.GetDual(iCat);
			return output;
		}
	}
}