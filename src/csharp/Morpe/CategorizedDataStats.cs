using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe
{
	public class CategorizedDataStats
	{
		/// <summary>
		/// The grand mean of all data (for each column of data).
		/// </summary>
		public double[] Mean;
		/// <summary>
		/// The grand variance of all data (for each column of data).
		/// </summary>
		public double[] Var;
		/// <summary>
		/// The category means of all data (for each column of data).  Indexed [iCat][iCoeff].
		/// </summary>
		public double[][] Means;
		/// <summary>
		/// The category variances of all data (for each column of data).  Indexed [iCat][iCoeff].
		/// </summary>
		public double[][] Vars;
		/// <summary>
		/// Calculates means and variances for the given data and weighting rule.
		/// </summary>
		/// <param name="data">Categorized data.</param>
		/// <param name="wRule">Weighting rule, provides an option to enforce equal priors.</param>
		public CategorizedDataStats(CategorizedData data, WeightingRule wRule)
		{
			int nCoeffs = data.X[0][0].Length;

			this.Mean = new double[nCoeffs];
			this.Var = new double[nCoeffs];
			this.Means = Util.NewArrays<double>(data.Ncats, nCoeffs);
			this.Vars = Util.NewArrays<double>(data.Ncats, nCoeffs);

			//-----------------------
			//	Compute CatWeights.
			//-----------------------
			double cwTotal = 0.0;
			double totalWeight = 0.0;
			double[] catWeights = new double[data.Ncats];
			if (wRule==WeightingRule.EqualPriors)
			{
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					double w = (double)data.Ntotal / (double)data.Neach[iCat] / (double)data.Ncats;
					catWeights[iCat] = w;
					cwTotal += w;
					totalWeight += w * (double)data.Neach[iCat];
				}
			}
			else if (wRule == WeightingRule.ObservedPriors)
			{
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					totalWeight += (float)data.Neach[iCat];
					catWeights[iCat] = 1.0f;
					cwTotal += 1.0f;
				}
			}
			else
				throw new ApplicationException("Unhandled weighting rule.");

			int iDatum;
			double mean, del, sum;

			//-----------------------
			//	Compute the means and variances.
			//-----------------------
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				for (int iCoeff = 0; iCoeff < nCoeffs; iCoeff++)
				{
					int nRows = data.X[iCat].Length;
					sum = 0.0;
					for (int iRow = 0; iRow < nRows; iRow++)
						sum += data.X[iCat][iRow][iCoeff];
					mean = sum / (double)nRows;
					sum = 0.0;
					for (int iRow = 0; iRow < nRows; iRow++)
					{
						del = data.X[iCat][iRow][iCoeff] - mean;
						sum += del * del;
					}
					this.Means[iCat][iCoeff] = mean;
					this.Vars[iCat][iCoeff] = sum / (double)(nRows - 1);

					this.Mean[iCoeff] += mean * catWeights[iCat];
				}
			}

			//-----------------------
			//	Finish computing the grand mean and variance.  Also get the parameter scale.
			//-----------------------
			float x;
			float[] xVec = new float[data.Ntotal];
			int[] idxVec = new int[data.Ntotal];
			int i025 = (int)(0.5 + 0.025 * (double)(data.Ntotal - 1));
			int i975 = (int)(0.5 + 0.025 * (double)(data.Ntotal - 1));
			sum = 0.0;
			for (int iCoeff = 0; iCoeff < nCoeffs; iCoeff++)
			{
				//	Grand mean
				mean = this.Mean[iCoeff] / cwTotal;
				this.Mean[iCoeff] = mean;

				//	Grand variance
				iDatum = 0;
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					int nRows = data.X[iCat].Length;
					for (int iRow = 0; iRow < nRows; iRow++)
					{
						x = data.X[iCat][iRow][iCoeff];
						del = x - mean;
						sum += del * catWeights[iCat];
						xVec[iDatum++] = x;
					}
				}
				this.Var[iCoeff] = sum / (totalWeight - 1.0f);
			}
		}
	}
}
