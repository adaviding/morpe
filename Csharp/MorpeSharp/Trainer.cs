using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System.Text;

using Morpe.Distributions.D1;

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
		public static readonly ThreadLocal<MonotonicRegressor> Regressor = new ThreadLocal<MonotonicRegressor>(() => new MonotonicRegressor());
		/// <summary>
		/// The classifier that is trained.
		/// </summary>
		public readonly Fitted<Classifier> Classifier;
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
			this.Classifier = new Fitted<Classifier>(classifier, double.NaN);
			this.CatWeights = new double[classifier.Ncats];
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
		/// The initial values of the polynomial coefficients.  Indexed as [iPoly][iCoeff].
		/// </summary>
		public float[][] ParamInit;
		/// <summary>
		/// Returns true if the trainer is currently training a classifier.
		/// </summary>
		public bool IsTraining { get { return this.trainerIsRunning; } }
		private bool trainerIsRunning = false;
        /// <summary>
        /// A container for analyses performed prior to optimization.
        /// </summary>
        public PreOptimizationAnalysis Analysis;
		/// <summary>
		/// Trains the classifier by optimizing parameters based on the training data using specified solver options.
		/// This can only be called once.
		/// </summary>
		/// <param name="data">The training data.
        /// WARNING:  In order to save memory, this data is altered in palce (instead of copying a new object).</param>
		/// <param name="ops">Sets the member variable <see cref="Classifier.Options"/>.  If a value is not provided, default options are used.</param>
		public async Task<Trainer> Train(CategorizedData data, SolverOptions ops)
		{
			lock (this)
			{
				if (this.trainerIsRunning)
					throw new ApplicationException("The Train method can only be called once at a time.  You might consider creating multiple Trainer objects.");
				this.trainerIsRunning = true;
			}
			try
			{
				if (data.Ncats != this.Classifier.Instance.Ncats)
					throw new ArgumentException("The number of categories in the classifier must be equal to the number of categories in the training set.");

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
				//	Condition the data and perform a polynomial expansion.
				//-----------------------
                this.Analysis = new PreOptimizationAnalysis();
                this.Analysis.Conditioner = SpaceConditioner.Measure(data);
                this.Analysis.Conditioner.Condition(data);
                data.Expand(this.Classifier.Instance.Coeffs);

				this.Analysis.ParamScale = new float[this.Classifier.Instance.Coeffs.Ncoeffs];
                this.Analysis.Xscale = new double[this.Classifier.Instance.Coeffs.Ncoeffs];

				if (ops.InitializeParams)
					this.Analysis.Crits = Static.NewArrays<UniCrit>(this.Classifier.Instance.Ncats, this.Classifier.Instance.Coeffs.Ncoeffs);

				//-----------------------
				//	Prepare category labels for each datum.
				//-----------------------
				byte[] catVec = new byte[data.Ntotal];
				int iDatum = 0;
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					int nRows = data.X[iCat].Length;
					for (int iRow = 0; iRow < nRows; iRow++)
						catVec[iDatum++] = (byte)iCat;
				}

				//-----------------------
				//	Get the parameter scale.
				//-----------------------
				double scale;
				float[] xVec = new float[data.Ntotal];
				int[] idxVec = new int[data.Ntotal];
				int i025 = (int)(0.5 + 0.025 * (double)(data.Ntotal - 1));
				int i500 = (int)(0.5 + 0.500 * (double)(data.Ntotal - 1));
				int i975 = (int)(0.5 + 0.975 * (double)(data.Ntotal - 1));
				for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
				{
					//	Quantiles determine the parameter scale.
					Static.FillSeries(idxVec);
					iDatum = 0;
					for (int iCat = 0; iCat < data.Ncats; iCat++)
					{
						int nSamp = data.Neach[iCat];
						for (int iSamp = 0; iSamp < nSamp; iSamp++)
						{
							xVec[iDatum++] = data.X[iCat][iSamp][iCoeff];
						}
					}
					Static.QuickSortIndex(idxVec, xVec, 0, xVec.Length - 1);
					scale = xVec[idxVec[i975]] - xVec[idxVec[i025]];
                    this.Analysis.Xscale[iCoeff] = scale;
                    this.Analysis.ParamScale[iCoeff] = (float)(1.0 / scale);

					if (ops.InitializeParams)
					{
						//	TODO:  This would work better if the median of x was subtracted out before the polynomial expansion was performed.
						//	I need to compute the median and then perform manual expansions.

						//	Get univariate classification criteria to get a first-order clue about the saliency of each feature.
						for (int iCat = 0; iCat < data.Ncats; iCat++)
							this.Analysis.Crits[iCat][iCoeff] = UniCrit.MaximumAccuracy(iCat, catVec, xVec, idxVec, this.CatWeights);
					}
				}
				xVec = null;


				//-----------------------
				//	Compute initial params.
				//-----------------------
				if (ops.InitializeParams)
				{
					//	Compute the expected minimum value for the univariate 2-category classification accuracy.
					double[] accMin = new double[data.Ncats];
					for (int iCat = 0; iCat < data.Ncats; iCat++)
					{
						accMin[iCat] = (this.totalWeight - this.CatWeights[iCat] * (double)data.Neach[iCat]) / this.totalWeight;
						if (accMin[iCat] < 0.5)
							accMin[iCat] = 1.0 - accMin[iCat];
					}

					//	The magnitude of each parameter is a function of univariate classification accuracy for the corresponding spatial dimension.
					double invNtotal = 1.0 / data.Ntotal;
					for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
					{
						if (this.Classifier.Instance.Npoly == 1)
						{
							double tAcc = 0.5 *
								(
									Math.Max(0.0, this.Analysis.Crits[0][iCoeff].Accuracy - accMin[0])
									+
                                    Math.Max(0.0, this.Analysis.Crits[1][iCoeff].Accuracy - accMin[1])
								);
							tAcc /= (1.0 - 0.5 * (accMin[0] + accMin[1]) + invNtotal);
                            if (this.Analysis.Crits[0][iCoeff].TargetUpper)
								this.ParamInit[0][iCoeff] = (float)tAcc * this.Analysis.ParamScale[iCoeff];
							else
                                this.ParamInit[0][iCoeff] = -(float)tAcc * this.Analysis.ParamScale[iCoeff];

						}
						else
						{
							for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
							{
                                double tAcc = Math.Max(0.0, this.Analysis.Crits[iPoly][iCoeff].Accuracy - accMin[iPoly]);
								tAcc /= (1.0 - accMin[iPoly] + invNtotal);
                                if (this.Analysis.Crits[iPoly][iCoeff].TargetUpper)
									this.ParamInit[iPoly][iCoeff] = (float)tAcc * this.Analysis.ParamScale[iCoeff];
								else
                                    this.ParamInit[iPoly][iCoeff] = -(float)tAcc * this.Analysis.ParamScale[iCoeff];
							}
						}
					}

					Static.Copy<float>(this.ParamInit, this.Classifier.Instance.Params);
				}
				else
				{
					//	Inherit parameters passed in by the classifier.  We assume the classifier was already initialized with parameters.
					Static.Copy<float>(this.Classifier.Instance.Params, this.ParamInit);
				}

				//-----------------------
				//	Perform dual optimizations.
				//-----------------------
				if (this.Classifier.Instance.Npoly > 2 && this.options.InitializeParams)
				{
					SolverOptions dualOps = this.options.Copy();
					dualOps.InitializeParams = false;
					dualOps.WeightingRule = WeightingRule.EqualPriors;
					Task<Trainer>[] dualTasks = new Task<Trainer>[this.Classifier.Instance.Npoly];
					for(int iPoly=0; iPoly<this.Classifier.Instance.Npoly; iPoly++)
					{
						Trainer t = new Trainer(this.Classifier.Instance.GetDual(iPoly));
						dualTasks[iPoly] = t.Train(data.GetDual(iPoly), dualOps);
					}
					for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
					{
						Trainer t = await dualTasks[iPoly];
						Array.Copy(t.Classifier.Instance.Params, this.ParamInit[iPoly], this.Classifier.Instance.Coeffs.Ncoeffs);
					}
					Static.Copy<float>(this.ParamInit, this.Classifier.Instance.Params);
				}

				//-----------------------
				//	Construct the classifiers based on parameters already specified.
				//-----------------------
				Static.FillSeries(idxVec);
				throw new ApplicationException("TO DO");

				//-----------------------
				//	Measure the classifier's fit.
				//-----------------------
				throw new ApplicationException("TO DO");

				//-----------------------
				//	Prepare optimization memory.
				//-----------------------
				int iBasis = 0;
				float[][][] basis = null;
				int nParams = this.Classifier.Instance.Npoly * this.Classifier.Instance.Coeffs.Ncoeffs;
				if (nParams <= 100)
					basis = this.randomDeviates();

				//-----------------------
				//	Optimize.
				//-----------------------
				throw new ApplicationException("TO DO");
			}
			finally
			{
				this.trainerIsRunning = false;
			}
		}
		
		/// <summary>
		/// This performs the quantization procedure for the data provided.
		/// </summary>
		/// <param name="yVals">The y-value of each datum in the sample.  These values are the output
		/// of the polynomial function for the target category.</param>
		/// <param name="yIdx">On output, provides the zero-based index into yVals that rank-orders the yVals.  On output,
		/// this order is calculated by calling <see cref="Static.QuickSortIndex"/> to sort yVals[yIdx[:]].
		/// On input, the order from the previous optimization step is passed in.  This means that the list is
		/// mostly sorted (for small parameter changes), and a mostly sorted list will make the optimization go
		/// faster than anothr random ordering such as 0...(yVals.Length-1).</param>
		/// <param name="cats">The category label of each datum.</param>
		/// <param name="catWeight">The weight assigned to each category.</param>
		/// <param name="targetCat">The target category for the y-values provided.</param>
		protected Quantization quantize(float[] yVals, int[] yIdx, byte[] catVec, double[] catWeight, byte targetCat)
		{
			//	Prepare output
			Quantization output = new Quantization(this.Nquantiles);
			double wPerBin = this.totalWeight / (double)this.Nquantiles;
			output.Pmin = 0.5 / wPerBin;
			output.Pmax = 1.0 - output.Pmin;

			//	IMPORTANT PERFORMANCE NOTE:
			//	Each time we enter this function, the sort index is preserved from the prior function call.
			//	This typically leads to faster sort times during optimization because typically the list is sorted
			//	already (or partially sorted).
			Static.QuickSortIndex(yIdx, yVals, 0, yVals.Length);

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
				int iSort = yIdx[iSamp];
				byte idCat = catVec[iSort];
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
						iSort = yIdx[--iSamp];
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
							iSort = yIdx[iSamp];
							idCat = catVec[iSort];
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
						output.Ysep[iBin] = (yVals[iSort] + yVals[yIdx[iSamp + 1]]) / 2.0;
					
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
		/// <summary>
		/// Generates an orthonormal basis of scaled random deviates based on a proper orthonormal matrix and the scale of each parameter.
		/// </summary>
		/// <returns>The orthonormal basis of random deviates, with entries multiplied by appropriate scale factors.</returns>
		protected float[][][] randomDeviates()
		{
			int nParams = this.Classifier.Instance.Npoly * this.Classifier.Instance.Coeffs.Ncoeffs;
			double[,] ortho = Static.RandomRotationMatrix(nParams);
			float[][][] output = Static.NewArrays<float>(nParams, this.Classifier.Instance.Npoly, this.Classifier.Instance.Coeffs.Ncoeffs);
			for (int iBasis = 0; iBasis < nParams; iBasis++)
			{
				for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
				{
					for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
					{
						int j = iPoly * this.Classifier.Instance.Coeffs.Ncoeffs + iCoeff;
						output[iBasis][iPoly][iCoeff] = (float)(ortho[iBasis, j] * this.Analysis.ParamScale[iCoeff]);
					}
				}
			}
			return output;
		}
		/// <summary>
		/// Generates a random deviate in the parameter space.
		/// </summary>
		/// <returns>The random deviate.</returns>
		protected float[][] randomDeviate()
		{
			int nParams = this.Classifier.Instance.Npoly * this.Classifier.Instance.Coeffs.Ncoeffs;
			float[][] output = Static.NewArrays<float>(this.Classifier.Instance.Npoly, this.Classifier.Instance.Coeffs.Ncoeffs);
			double sumsq = 0.0;
			double r;
			for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
			{
				for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
				{
					r = Gaussian.InvCdf(Static.Rand.NextDouble());
					sumsq += r * r;
					output[iPoly][iCoeff] = (float)r;
				}
			}
			float norm = (float)(1.0/sumsq);
			for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
			{
				for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
				{
					output[iPoly][iCoeff] *= this.Analysis.ParamScale[iCoeff] * norm;
				}
			}
			return output;
		}
	}
}
