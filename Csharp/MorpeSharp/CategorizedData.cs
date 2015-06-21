using System;

namespace Morpe
{
	/// <summary>
	/// Multivariate data that is segregated by category membership.
	/// </summary>
	public class CategorizedData
	{
		/// <summary>
		/// The data, indexed as [c][i][j].  Each page c is a category.  Each column i is a unique data point.  Each row j is a spatial axis.
		/// </summary>
		public readonly float[][][] X;
		/// <summary>
		/// The total number of data points (across all categories).
		/// </summary>
		public readonly int Ntotal;
		/// <summary>
		/// The number of data points in each category.  The number of categories is equal to <see cref="MorpeSharp.Neach.Length"/>;
		/// </summary>
		public readonly int[] Neach;
		/// <summary>
		/// The spatial dimensionality, and the number of rows of <see cref="MorpeSharp.X"/>.
		/// </summary>
		public readonly int Ndims;
		/// <summary>
		/// The number of categories, also equal to <see cref="MorpeSharp.Neach.Length"/>.
		/// </summary>
		public readonly int Ncats;
		/// <summary>
		/// Initializes a new instance of the <see cref="MorpeSharp.CategorizedData"/> class.  Allocates memory for pages, rows, and columns of <see cref="MorpeSharp.X"/>.
		/// </summary>
		/// <param name="Neach"><see cref="MorpeSharp.Neach"/></param>
		/// <param name="Ndims"><see cref="MorpeSharp.Ndims"/></param>
		public CategorizedData(int[] Neach, int Ndims)
		{
			this.Neach = Neach;
			this.Ndims = Ndims;
			this.Ncats = this.Neach.Length;
			this.Ntotal = Static.Sum(this.Neach);
			this.X = new float[this.Ncats][][];
			for (int c=0; c<this.Ncats; c++)
			{
				this.X[c] = new float[this.Neach[c]][];
				int n = this.Neach[c];
				for(int i=0; i<n; i++)
					this.X[c][i] = new float[this.Ndims];
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
						this.X[iCat][iSamp] = poly.Expand(Static.GetSubarray(x, 0, this.Ndims - 1));
					else
						this.X[iCat][iSamp] = poly.Expand(x);
				}
			}
		}
		/// <summary>
		/// This reverses the effects of <see cref="Expand"/>.  Expanded terms are removed if they exist.
		/// </summary>
		/// <param name="poly">The defiinition of the polynomial expansion.</param>
		public void Contract()
		{
			for (int iCat = 0; iCat < this.Ncats; iCat++)
			{
				int nSamp = this.Neach[iCat];
				for (int iSamp = 0; iSamp < nSamp; iSamp++)
				{
					float[] x = this.X[iCat][iSamp];
					if (x.Length > this.Ndims)
						this.X[iCat][iSamp] = Static.GetSubarray(x, 0, this.Ndims - 1);
				}
			}
		}
	}
}