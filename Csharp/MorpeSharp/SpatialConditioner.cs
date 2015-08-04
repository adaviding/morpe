using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
		public int Ndims;
		/// <summary>
		/// Allocates memory for a new Spatial conditioner with the given dimensionality.
		/// </summary>
		/// <param name="nDims">The spatial dimensionality.</param>
		public SpatialConditioner(int nDims)
		{
			this.Ndims = nDims;
			this.Origin = new float[nDims];
			this.Spread = new float[nDims];
		}
		/// <summary>
		/// Conditions the data.
		/// </summary>
		/// <param name="data"></param>
		public void Condition(CategorizedData data)
		{
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				for (int iRow = 0; iRow < data.Neach[iCat]; iRow++)
				{
					float[] xOld = data.X[iCat][iRow];
					float[] xNew = new float[data.Ndims];
					for (int iDim = 0; iDim < data.Ndims; iDim++)
					{
						xNew[iDim] = (xOld[iDim] - this.Origin[iDim]) / this.Spread[iDim];
					}
				}
			}
		}
	}
}
