using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe
{
	/// <summary>
	/// The SpatialConditionMeasurer is measured at the beginning of optimization, prior to performing the first polynomial expansion on the training data.
	/// It is used to condition the feature space.  In other words, it is used to ensure that the feature space is centered on 0 and has
	/// spread of 1 in all directions.  (The polynomial coefficients are defined for the conditioned feature space.)
	/// </summary>
	public class SpatialConditionMeasurer
	{
		/// <summary>
		/// The average medain of the unconditioned training data.  This variable has this.Ndims columns.  It is the average of this.Medians.
		/// </summary>
		public float[] AvgMedian;
		/// <summary>
		/// The measured spread of the unconditioned training data.  This variable has this.Ndims columns.  It is the RMS deviation from the AvgMedian
		/// for all data.  All categories have equal influence over the Spread, regardless of their base rates.
		/// </summary>
		public float[] Spread;
		/// <summary>
		/// The medians of the data within each category, prior to conditioning.
		/// This variable has this.Ncats rows and this.Ndims columns.
		/// </summary>
		public float[][] Medians;
		/// <summary>
		/// The RMS deviation from the median for each category, prior to conditioning.
		/// This variable has this.Ncats rows and this.Ndims columns.
		/// </summary>
		public float[][] Spreads;
		/// <summary>
		/// The number of categories.
		/// </summary>
		public int Ncats;
		/// <summary>
		/// The number of spatial dimensions.
		/// </summary>
		public int Ndims;
		/// <summary>
		/// Allocates memory for a specified number of categories and dimensions.
		/// </summary>
		/// <param name="nCats">The number of categories.</param>
		/// <param name="nDims">The number of spatial dimensions (unexpanded).</param>
		public SpatialConditionMeasurer(int nCats, int nDims)
		{
			this.Ncats = nCats;
			this.Ndims = nDims;
			this.AvgMedian = new float[nDims];
			this.Spread = new float[nDims];
			this.Medians = Static.NewArrays<float>(nCats, nDims);
			this.Spreads = Static.NewArrays<float>(nCats, nDims);
		}
		/// <summary>
		/// Measures a space conditioner for the training data.
		/// </summary>
		/// <param name="data">The training data.</param>
		/// <returns>The space conditioner which has been measured for the training data.</returns>
		public static SpatialConditionMeasurer Measure(CategorizedData data)
		{
			if (data == null)
				return null;
			SpatialConditionMeasurer output = new SpatialConditionMeasurer(data.Ncats, data.Ndims);

			float temp;

			int[] idxVec = null;
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				if (idxVec == null || idxVec.Length != data.Neach[iCat])
					idxVec = new int[data.Neach[iCat]];
				bool isOdd = data.Neach[iCat] % 2 == 1;
				int iMed = data.Neach[iCat] / 2;
				for (int iCol = 0; iCol < data.Ndims; iCol++)
				{
					Static.FillSeries(idxVec);
					Static.QuickSortIndex(idxVec, data.X[iCat], iCol, 0, data.Neach[iCat] - 1);

					if (isOdd)
						output.Medians[iCat][iCol] = temp = data.X[iCat][iMed][iCol];
					else
						output.Medians[iCat][iCol] = temp = (data.X[iCat][iMed-1][iCol] + data.X[iCat][iMed][iCol]) / 2.0f;
					output.AvgMedian[iCol] += temp;
				}
			}
			for (int iCol = 0; iCol < data.Ndims; iCol++)
			{
				output.AvgMedian[iCol] /= (float)data.Ncats;
				for(int iCat = 0; iCat < data.Ncats; iCat++)
				{
					double ssMedian = 0.0f;
					double ssOrigin = 0.0f;
					double x, dx;
					int nRows = data.Neach[iCat];
					for(int iRow=0; iRow<nRows; iRow++)
					{
						x = data.X[iCat][iRow][iCol];

						dx = x - output.Medians[iCat][iCol];
						ssMedian += dx * dx;

						dx = x - output.AvgMedian[iCol];
						ssOrigin += dx * dx;
					}
					dx = 1.0/(double)(nRows-1);
					output.Spreads[iCat][iCol] = (float)Math.Sqrt(dx * ssMedian);
					output.Spread[iCol] += (float)Math.Sqrt(dx * ssOrigin)/(float)data.Ncats;
				}
			}
			return output;
		}
		public SpatialConditioner Conditioner()
		{
			SpatialConditioner output = new SpatialConditioner(this.Ndims);
			Array.Copy(this.AvgMedian, output.Origin, this.Ndims);
			Array.Copy(this.Spread, output.Spread, this.Ndims);
			return output;
		}
	}
}
