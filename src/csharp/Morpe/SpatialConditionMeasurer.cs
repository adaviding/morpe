using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using F1 = Morpe.Numerics.F1;
using I1 = Morpe.Numerics.I1;

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
        /// Measures a space conditioner for the training data.
        /// </summary>
        /// <param name="data">The training data.</param>
        /// <returns>The space conditioner which has been measured for the training data.</returns>
        [return: NotNull]
        public static SpatialConditionMeasurer Measure(
            [NotNull] CategorizedData data)
        {
            SpatialConditionMeasurer output = new SpatialConditionMeasurer(data.NumCats, data.NumDims);

            float temp;

            int[] idxVec = null;
            for (int iCat = 0; iCat < data.NumCats; iCat++)
            {
                if (idxVec == null || idxVec.Length != data.NumEach[iCat])
                    idxVec = new int[data.NumEach[iCat]];
                bool isOdd = data.NumEach[iCat] % 2 == 1;
                int iMed = data.NumEach[iCat] / 2;
                for (int iCol = 0; iCol < data.NumDims; iCol++)
                {
                    I1.Util.FillSeries(idxVec);
                    F1.Util.QuickSortIndex(idxVec, data.X[iCat], iCol, 0, data.NumEach[iCat] - 1);

                    if (isOdd)
                        output.Medians[iCat][iCol] = temp = data.X[iCat][iMed][iCol];
                    else
                        output.Medians[iCat][iCol] = temp = (data.X[iCat][iMed-1][iCol] + data.X[iCat][iMed][iCol]) / 2.0f;
                    output.AvgMedian[iCol] += temp;
                }
            }
            for (int iCol = 0; iCol < data.NumDims; iCol++)
            {
                output.Spread[iCol] = 0.0f;
                output.AvgMedian[iCol] /= (float)data.NumCats;
                for(int iCat = 0; iCat < data.NumCats; iCat++)
                {
                    double ssMedian = 0.0f;
                    double ssOrigin = 0.0f;
                    double x, dx;
                    int nRows = data.NumEach[iCat];
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
                    output.Spread[iCol] += (float)Math.Sqrt(dx * ssOrigin)/(float)data.NumCats;
                }
            }
            return output;
        }
        
        /// <summary>
        /// The average medain of the unconditioned training data.  This variable has this.Ndims columns.  It is the average of this.Medians.
        /// </summary>
        public float[] AvgMedian;
        
        /// <summary>
        /// The medians of the data within each category, prior to conditioning.
        /// This variable has this.Ncats rows and this.Ndims columns.
        /// </summary>
        public float[][] Medians;
        
        /// <summary>
        /// The measured spread of the unconditioned training data.  This variable has this.Ndims columns.  It is the RMS deviation from the AvgMedian
        /// for all data.  All categories have equal influence over the Spread, regardless of their base rates.
        /// </summary>
        public float[] Spread;
        
        /// <summary>
        /// The RMS deviation from the median for each category, prior to conditioning.
        /// This variable has this.Ncats rows and this.Ndims columns.
        /// </summary>
        public float[][] Spreads;
        
        /// <summary>
        /// The number of categories.
        /// </summary>
        public int NumCats;
        
        /// <summary>
        /// The number of spatial dimensions.
        /// </summary>
        public int NumDims;
        
        /// <summary>
        /// Allocates memory for a specified number of categories and dimensions.
        /// </summary>
        /// <param name="numCats">The number of categories.</param>
        /// <param name="numDims">The number of spatial dimensions (unexpanded).</param>
        public SpatialConditionMeasurer(int numCats, int numDims)
        {
            this.NumCats = numCats;
            this.NumDims = numDims;
            this.AvgMedian = new float[numDims];
            this.Spread = new float[numDims];
            this.Medians = Util.NewArrays<float>(numCats, numDims);
            this.Spreads = Util.NewArrays<float>(numCats, numDims);
        }
        
        /// <summary>
        /// Creates a spatial conditioner based on the measurement.
        /// </summary>
        /// <returns>The new spatial conditioner.</returns>
        public SpatialConditioner Conditioner()
        {
            SpatialConditioner output = new SpatialConditioner(this.NumDims);
            Array.Copy(this.AvgMedian, output.Origin, this.NumDims);
            Array.Copy(this.Spread, output.Spread, this.NumDims);
            return output;
        }
        
        /// <summary>
        /// Deep copy.
        /// </summary>
        /// <returns>A deep copy.</returns>
        public SpatialConditionMeasurer Clone()
        {
            SpatialConditionMeasurer output = new SpatialConditionMeasurer(this.NumCats, this.NumDims);
            Array.Copy(this.AvgMedian, output.AvgMedian, this.NumDims);
            Array.Copy(this.Spread, output.Spread, this.NumDims);
            Util.Copy<float>(this.Medians, output.Medians);
            Util.Copy<float>(this.Spreads, output.Spreads);
            return output;
        }
    }
}
