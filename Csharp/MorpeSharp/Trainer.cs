using System;
using System.Collections.Generic;
using System.Threading;
using System.Text;

namespace Morpe
{
	/// <summary>
	/// Trains an instance of the MoRPE classifier using the training data provided.
	/// </summary>
	public class Trainer
	{
		/// <summary>
		/// Used to perform monotonic regression of a tabulated function.
		/// </summary>
		public static readonly ThreadLocal<F.MonotonicRegressor> Regressor = new ThreadLocal<F.MonotonicRegressor>(() => new F.MonotonicRegressor());
		/// <summary>
		/// The classifier that is trained.
		/// </summary>
		public readonly Classifier Classifier;
		/// <summary>
		/// The relative importance of each category on the solution.  This is important when calculating accuracy and entropy.
		/// </summary>
		public readonly double[] CatWeights;
		/// <summary>
		/// Solver options used during training.  Prior to training, this member contains the default options.
		/// </summary>
		public SolverOptions Options { get { return this.options; } }
		protected SolverOptions options = null;
		/// <summary>
		/// Construct a trainer for a given classifier.
		/// </summary>
		/// <param name="classifier">The classifier to be trained.</param>
		public Trainer(Classifier classifier)
		{
			this.Classifier = classifier;
			this.CatWeights = new double[classifier.Ncats];

			//	Initialize thread-local storage for training.
			this.localSortIdx = new ThreadLocal<int[]>[classifier.Npoly];
			for (int i = 0; i < classifier.Npoly; i++)
				this.localSortIdx[i] = new ThreadLocal<int[]>();
		}
		/// <summary>
		/// The total amount of weight across all training data.
		/// </summary>
		protected double totalWeight;
		/// <summary>
		/// The number of quantiles used.
		/// </summary>
		public int Nquantiles = 15;
		/// <summary>
		/// The scale of each polynomial coefficient.  The scales are the same across all polynomials.
		/// </summary>
		public float[] ParamScales;
		/// <summary>
		/// The initial values of the polynomial coefficients.  Indexed as [iPoly][iCoeff].
		/// </summary>
		public float[][] ParamInit;
		/// <summary>
		/// The grand mean of all training data (for each column of data).
		/// </summary>
		public double[] Xmean;
		/// <summary>
		/// The grand variance of all training data (for each column of data).
		/// </summary>
		public double[] Xvar;
		/// <summary>
		/// The scale of the training data (for each column of data).
		/// </summary>
		public double[] Xscale;
		/// <summary>
		/// The category means of all training data (for each column of data).  Indexed [iCat][iCoeff].
		/// </summary>
		public double[][] Xmeans;
		/// <summary>
		/// The category variances of all training data (for each column of data).  Indexed [iCat][iCoeff].
		/// </summary>
		public double[][] Xvars;
		/// <summary>
		/// The scales of the training data for each category and each column of data.  Indexed as [iCat][iCoeff].
		/// </summary>
		public double[][] Xscales;
		/// <summary>
		/// Univariate classifiers for each category and each column of data.  Indexed as [iCat][iCoeff].
		/// </summary>
		public UniCrit[][] Xcrits;
		/// <summary>
		/// Trains the classifier by optimizing parameters based on the training data using specified solver options.
		/// </summary>
		/// <param name="data">The training data.</param>
		/// <param name="ops">Sets the member variable <see cref="Classifier.Options"/>.  If a value is not provided, default options are used.</param>
		public void Train(CategorizedData data, SolverOptions ops)
		{
			//-----------------------
			//	Set the solver otions.
			//-----------------------
			if (ops == null)
				ops = new SolverOptions();
			this.options = ops;

			//-----------------------
			//	Compute CatWeights.
			//-----------------------
			double cwTotal = 0.0;
			this.totalWeight = 0.0;  // Total weight across all training data.
			if (this.options.WeightingRule == WeightingRule.EqualPriors)
			{
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					double w = (double)data.Ntotal / (double)data.Neach[iCat] / (double)data.Ncats;
					this.CatWeights[iCat] = w;
					cwTotal += w;
					this.totalWeight += w * (double)data.Neach[iCat];
				}
			}
			else if (this.options.WeightingRule == WeightingRule.ObservedPriors)
			{
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					this.totalWeight += (float)data.Neach[iCat];
					this.CatWeights[iCat] = 1.0f;
					cwTotal += 1.0f;
				}
			}
			else
				throw new ApplicationException("Unhandled weighting rule.");

			//-----------------------
			//	Perform the polynomial expansion of the data if necessary.
			//-----------------------
			if (data.X[0][0].Length < this.Classifier.Coeffs.Ncoeffs)
				data.Expand(this.Classifier.Coeffs);

			this.ParamScales = new float[this.Classifier.Coeffs.Ncoeffs];

			this.Xmean = new double[this.Classifier.Coeffs.Ncoeffs];
			this.Xvar = new double[this.Classifier.Coeffs.Ncoeffs];
			this.Xscale = new double[this.Classifier.Coeffs.Ncoeffs];
			this.Xmeans = Static.NewArrays<double>(this.Classifier.Ncats, this.Classifier.Coeffs.Ncoeffs);
			this.Xvars = Static.NewArrays<double>(this.Classifier.Ncats, this.Classifier.Coeffs.Ncoeffs);
			this.Xscales = Static.NewArrays<double>(this.Classifier.Ncats, this.Classifier.Coeffs.Ncoeffs);
			this.Xcrits = Static.NewArrays<UniCrit>(this.Classifier.Ncats, this.Classifier.Coeffs.Ncoeffs);

			double mean, del, sum;

			//-----------------------
			//	Compute the means and variances.
			//-----------------------
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				for(int iCoeff=0; iCoeff<this.Classifier.Coeffs.Ncoeffs; iCoeff++)
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
						sum += del*del;
					}
					this.Xmeans[iCat][iCoeff] = mean;
					this.Xvars[iCat][iCoeff] = sum / (double)(nRows-1);

					this.Xmean[iCoeff] += mean * this.CatWeights[iCat];
				}
			}

			//-----------------------
			//	Prepare category labels for each datum.
			//-----------------------
			int[] catVec = new int[data.Ntotal];
			int iDatum = 0;
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				int nRows = data.X[iCat].Length;
				for (int iRow = 0; iRow < nRows; iRow++)
					catVec[iDatum++] = iCat;
			}

			//-----------------------
			//	Compute the expected minimum value for the univariate 2-category classification accuracy.
			//-----------------------
			double[] accMin = new double[data.Ncats];
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				accMin[iCat] = (this.totalWeight - this.CatWeights[iCat] * (double)data.Neach[iCat]) / this.totalWeight;
				if (accMin[iCat] < 0.5)
					accMin[iCat] = 1.0 - accMin[iCat];
			}

			//-----------------------
			//	Finish computing the grand mean and variance.  Also get the parameter scale.
			//-----------------------
			double scale;
			float x;
			float[] xVec = new float[data.Ntotal];
			int[] idxVec = new int[data.Ntotal];
			int i025 = (int)(0.5 + 0.025 * (double)(data.Ntotal - 1));
			int i975 = (int)(0.5 + 0.025 * (double)(data.Ntotal - 1));
			sum = 0.0;
			for (int iCoeff = 0; iCoeff < this.Classifier.Coeffs.Ncoeffs; iCoeff++)
			{
				//	Grand mean
				mean = this.Xmean[iCoeff] / cwTotal;
				this.Xmean[iCoeff] = mean;

				//	Grand variance
				iDatum = 0;
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					int nRows = data.X[iCat].Length;
					for (int iRow = 0; iRow < nRows; iRow++)
					{
						x = data.X[iCat][iRow][iCoeff];
						del = x - mean;
						sum += del * this.CatWeights[iCat];
						xVec[iDatum++] = x;
					}
				}
				this.Xvar[iCoeff] = sum / (this.totalWeight - 1.0f);

				//	Quantiles determine the parameter scale.
				Static.FillSeries(idxVec);
				Static.QuickSortIndex(idxVec, xVec, 0, xVec.Length - 1);
				scale = xVec[idxVec[i975]] - xVec[idxVec[i025]];
				this.Xscale[iCoeff] = scale;
				this.ParamScales[iCoeff] = (float)(1.0 / scale);

				//	Get univariate classification criteria to get a first-order clue about the saliency of each feature.
				for (int iCat = 0; iCat < data.Ncats; iCat++)
					this.Xcrits[iCat][iCoeff] = UniCrit.MaximumAccuracy(iCat, catVec, xVec, idxVec, this.CatWeights);
			}

			//-----------------------
			//	Compute initial params.
			//-----------------------
			if (ops.InitializeParams)
			{
				//	The magnitude of each parameter is a function of univariate classification accuracy for the corresponding spatial dimension.

				double invNtotal = 1.0 / data.Ntotal;
				for (int iCoeff = 0; iCoeff < this.Classifier.Coeffs.Ncoeffs; iCoeff++)
				{
					if (this.Classifier.Npoly == 1)
					{
						double tAcc = 0.5 *
							(
								Math.Max(0.0,this.Xcrits[0][iCoeff].Accuracy - accMin[0])
								+
								Math.Max(0.0,this.Xcrits[1][iCoeff].Accuracy - accMin[1])
							);
						tAcc /= (1.0 - 0.5 * (accMin[0] + accMin[1]) + invNtotal);
						if (this.Xcrits[0][iCoeff].TargetUpper)
							this.ParamInit[0][iCoeff] =  (float)tAcc * this.ParamScales[iCoeff];
						else
							this.ParamInit[0][iCoeff] = -(float)tAcc * this.ParamScales[iCoeff];
							
					}
					else
					{
						for (int iPoly = 0; iPoly < this.Classifier.Npoly; iPoly++)
						{
							double tAcc = Math.Max(0.0, this.Xcrits[iPoly][iCoeff].Accuracy - accMin[iPoly]);
							tAcc /= (1.0 - accMin[iPoly] + invNtotal);
							if (this.Xcrits[iPoly][iCoeff].TargetUpper)
								this.ParamInit[iPoly][iCoeff] =  (float)tAcc * this.ParamScales[iCoeff];
							else
								this.ParamInit[iPoly][iCoeff] = -(float)tAcc * this.ParamScales[iCoeff];
						}
					}
				}
			}
			else
			{
				//	Inherit parameters passed in by the classifier.  We assume the classifier was already initialized with parameters.
				Static.Copy(this.Classifier.Params, this.ParamInit);
			}

			//	Prepare optimization memory.
			byte[] cats = new byte[data.Ntotal];
			double[] y = new double[data.Ntotal];

			//	Set the category labels for each datum.
			int jSamp = 0;
			for (int iCat = 0; iCat < data.Ncats; iCat++)
			{
				int n = data.Neach[iCat];
				for (int iSamp = 0; iSamp < n; iSamp++)
					cats[jSamp++] = (byte)iCat;
			}

			//	TODO:  Optimize parameters.
		}

		/// <summary>
		/// Internal memory resource.  Used for sorting y-values during the quantization procedure.  Indexed for every [iPoly].
		/// </summary>
		protected ThreadLocal<int[]>[] localSortIdx;
		/// <summary>
		/// This performs the quantization procedure for the data provided.
		/// </summary>
		/// <param name="yVals">The y-value of each datum in the sample.  These values are the output
		/// of the polynomial function for the target category.</param>
		/// <param name="cats">The category label of each datum.</param>
		/// <param name="catWeight">The weight assigned to each category.</param>
		/// <param name="targetCat">The target category for the y-values provided.</param>
		protected Quantization quantize(float[] yVals, byte[] cats, double[] catWeight, byte targetCat)
		{
			//	Prepare output
			Quantization output = new Quantization(this.Nquantiles);
			double wPerBin = this.totalWeight / (double)this.Nquantiles;
			output.Pmin = 0.5 / wPerBin;
			output.Pmax = 1.0 - output.Pmin;

			//	Get the memory resource for sorting yVals.
			ThreadLocal<int[]> localIdx = this.localSortIdx[targetCat];
			if (localIdx.Value == null || localIdx.Value.Length < yVals.Length)
			{
				localIdx.Value = new int[yVals.Length];
				Static.FillSeries(localIdx.Value);
			}
			//	IMPORTANT PERFORMANCE NOTE:
			//	Each time we enter this function, the sort index is preserved from the prior function call.
			//	This typically leads to faster sort times during optimization because typically the list is sorted
			//	already (or partially sorted).
			int[] idx = localIdx.Value;

			Static.QuickSortIndex(idx, yVals, 0, yVals.Length);

			int iBin = 0;  // The current bin.
			double ySum = 0.0; // Weighted sum of y-values for the current bin.
			double correctWeight = 0.0; // Correct weight for the current bin.
			double errorWeight = 0.0;  // Incorrect weight for the current bin.
			double binWeight;
			
			double wNextBin = wPerBin; // The amount of cucmulative weight that separates the current bin from the next.
			double wThis = 0.0; // The current accumulated weight.
			double dwThis; // The amount of weight added by the current sample.
			double wLast;	//	The accumulated weight prior to the current sample.

			//	Quantize the  samples.
			for (int iSamp = 0; iSamp < yVals.Length; iSamp++)
			{
				int iSort = idx[iSamp];
				byte idCat = cats[iSort];
				wLast = wThis;
				dwThis = catWeight[idCat];
				wThis += dwThis;
				if (wThis > wNextBin || iSamp == yVals.Length-1)
				{
					//---------------------------
					//	It is time for a new bin.
					//---------------------------
					if ( iBin < this.Nquantiles-1 && wNextBin - wLast > wThis - wNextBin)
					{
						//	Rewind to last sample.
						iSort = idx[--iSamp];
						wThis = wLast;
					}
					else
					{
						//	Process this sample  as normal.
						ySum += dwThis * yVals[iSort];
						if (idCat == targetCat)
							correctWeight += dwThis;
						else
							errorWeight += dwThis;
					}

					//---------------------------
					//	Special processing for the last bin.
					//---------------------------
					if (iBin == this.Nquantiles - 1)
					{
						//	There should not be any more samples, but we are doing this just to be sure.
						while (++iSamp < yVals.Length)
						{
							iSort = idx[iSamp];
							idCat = cats[iSort];
							dwThis = catWeight[idCat];
							wThis += dwThis;
							ySum += dwThis * yVals[iSort];
							if (idCat == targetCat)
								correctWeight += dwThis;
							else
								errorWeight += dwThis;
						}
					}

					//---------------------------
					//	Close out the current bin.
					//---------------------------
					binWeight = Math.Min(0.00001,correctWeight + errorWeight);
					output.P[iBin] = correctWeight / binWeight;
					output.Ymid[iBin] = ySum / binWeight;
					if (iBin < this.Nquantiles-1)
						output.Ysep[iBin] = (yVals[iSort] + yVals[idx[iSamp + 1]]) / 2.0;
					
					//---------------------------
					//	Start a next bin.
					//---------------------------
					iBin++;
					correctWeight = errorWeight = 0.0;
					wNextBin += wPerBin;
				}
				else
				{
					//	Process each sample.
					ySum += dwThis * yVals[iSort];
					if (idCat == targetCat)
						correctWeight += dwThis;
					else
						errorWeight += dwThis;
				}
			}
			return output;
		}
	}
}
