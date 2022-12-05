using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System.Text;

using D  = Morpe.Numerics.D;
using D1 = Morpe.Numerics.D1;
using F1 = Morpe.Numerics.F1;
using I1 = Morpe.Numerics.I1;

namespace Morpe
{
	/// <summary>
	/// Trains an instance of the MoRPE classifier using the training data provided.
	/// </summary>
	public class Trainer
	{
		/// <summary>
		/// The relative importance of each category on the solution.  This is important when calculating accuracy and entropy.
		/// </summary>
		public double[] CatWeights;
		/// <summary>
		/// The trained classifier.
		/// </summary>
		public Fitted<Classifier> Classifier;
		/// <summary>
		/// Returns true if the trainer is currently training a classifier.
		/// </summary>
		public bool IsTraining { get { return this.trainerIsRunning; } }
		private bool trainerIsRunning = false;
		/// <summary>
		/// The number of quantiles used.
		/// </summary>
		public int Nquantiles = 15;
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
		/// A container for the output of analyses performed prior to optimization.
		/// </summary>
		public PreOptimizationAnalysis Analysis;
		/// <summary>
		/// Trains the classifier by optimizing parameters based on the training data using specified solver options.
		/// This can only be called once.
		/// </summary>
		/// <param name="data">The training data.
		/// WARNING:  In order to save memory, this data is altered in palce (instead of copying a new object).</param>
		/// <param name="ops">Sets the member variable <see cref="Classifier.Options"/>.  If a value is not provided, default options are used.</param>
		/// <param name="ops">The pre-optimization analysis <see cref="Classifier.Analysis"/>.  If a value is not provided, default options are used.</param>
		public async Task<Trainer> Train(CategorizedData data, SolverOptions ops=null, PreOptimizationAnalysis analysis=null)
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

				//	[iDatum] Used for identifying the category label for each datum.
				byte[] catVec = new byte[data.Ntotal];
				//	[iDatum] Used for storing the polynomial calculation for each datum.
				float[][] yVec = Util.NewArrays<float>(this.Classifier.Instance.Npoly, data.Ntotal);
				//	[iPoly][iDatum]  Used for sorting data rows.  Sort order is preserved for each polynomial.
				int[][] idxVec = Util.NewArrays<int>(this.Classifier.Instance.Npoly,data.Ntotal);

				//-----------------------
				//	Prepare category labels for each datum.
				//-----------------------
				int iDatum = 0;
				for (int iCat = 0; iCat < data.Ncats; iCat++)
				{
					int nRows = data.X[iCat].Length;
					for (int iRow = 0; iRow < nRows; iRow++)
						catVec[iDatum++] = (byte)iCat;
				}

				if (analysis != null)
					this.Analysis = analysis;
				else
				{
					//-----------------------
					//	Condition the data and perform a polynomial expansion.
					//-----------------------
					this.Analysis = new PreOptimizationAnalysis();
					this.Analysis.ConditionMeasurer = SpatialConditionMeasurer.Measure(data);
					this.Analysis.Conditioner = this.Analysis.ConditionMeasurer.Conditioner();
					this.Analysis.Conditioner.Condition(data);
					data.Expand(this.Classifier.Instance.Coeffs);

					this.Analysis.ParamScale = new float[this.Classifier.Instance.Coeffs.Ncoeffs];
					this.Analysis.ParamInit = Util.NewArrays<float>(this.Classifier.Instance.Ncats, this.Classifier.Instance.Coeffs.Ncoeffs);

					//if (ops.InitializeParams)
					this.Analysis.Crits = Util.NewArrays<UniCrit>(this.Classifier.Instance.Ncats, this.Classifier.Instance.Coeffs.Ncoeffs);

					//-----------------------
					//	Get the parameter scale.
					//-----------------------
					float scale;
					int i15 = (int)(0.5 + 0.15 * (double)(data.Ntotal - 1));
					int i50 = (int)(0.5 + 0.50 * (double)(data.Ntotal - 1));
					int i85 = (int)(0.5 + 0.85 * (double)(data.Ntotal - 1));
					float[] xVec = yVec[0];
					for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
					{
						//	Quantiles determine the parameter scale.
						I1.Util.FillSeries(idxVec[0]);
						iDatum = 0;
						for (int iCat = 0; iCat < data.Ncats; iCat++)
						{
							int nSamp = data.Neach[iCat];
							for (int iSamp = 0; iSamp < nSamp; iSamp++)
								xVec[iDatum++] = data.X[iCat][iSamp][iCoeff];
						}
						F1.Util.QuickSortIndex(idxVec[0], xVec, 0, xVec.Length - 1);
						scale = xVec[idxVec[0][i85]] - xVec[idxVec[0][i15]];
						this.Analysis.ParamScale[iCoeff] = (float)(1.0 / scale);

						//if (ops.InitializeParams)
						//{
						//	Get univariate classification criteria to get a first-order clue about the saliency of each feature.
						for (int iCat = 0; iCat < data.Ncats; iCat++)
							this.Analysis.Crits[iCat][iCoeff] = UniCrit.MaximumAccuracy(iCat, catVec, xVec, idxVec[0], this.CatWeights);
						//}
					}

					//-----------------------
					//	Compute initial params.
					//-----------------------
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
								this.Analysis.ParamInit[0][iCoeff] = (float)tAcc * this.Analysis.ParamScale[iCoeff];
							else
								this.Analysis.ParamInit[0][iCoeff] = -(float)tAcc * this.Analysis.ParamScale[iCoeff];

						}
						else
						{
							for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
							{
								double tAcc = Math.Max(0.0, this.Analysis.Crits[iPoly][iCoeff].Accuracy - accMin[iPoly]);
								tAcc /= (1.0 - accMin[iPoly] + invNtotal);
								if (this.Analysis.Crits[iPoly][iCoeff].TargetUpper)
									this.Analysis.ParamInit[iPoly][iCoeff] = (float)tAcc * this.Analysis.ParamScale[iCoeff];
								else
									this.Analysis.ParamInit[iPoly][iCoeff] = -(float)tAcc * this.Analysis.ParamScale[iCoeff];
							}
						}
					}
				}

				//-----------------------
				//	Set classifier parameters.
				//-----------------------
				if (ops.InitializeParams)
				{
					Util.Copy<float>(this.Analysis.ParamInit, this.Classifier.Instance.Params);
				}
				else
				{
					//	Inherit parameters passed in by the classifier.  We assume the classifier was already initialized with parameters.
					Util.Copy<float>(this.Classifier.Instance.Params, this.Analysis.ParamInit);
				}

				//-----------------------
				//	Perform dual optimizations.
				//-----------------------
				if (this.Classifier.Instance.Npoly > 2 && this.options.InitializeParams)
				{
					SolverOptions dualOps = this.options.Copy();
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
						Array.Copy(t.Classifier.Instance.Params, this.Analysis.ParamInit[iPoly], this.Classifier.Instance.Coeffs.Ncoeffs);
					}
					Util.Copy<float>(this.Analysis.ParamInit, this.Classifier.Instance.Params);
				}

				//-----------------------
				//	Finish constructing the classifier for the first time.
				//-----------------------
				D1.MonotonicRegressor regressor = new D1.MonotonicRegressor();
				for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
					I1.Util.FillSeries(idxVec[iPoly]);

				for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
				{
					//	Set the quantized probability limits.
					double nPerQuantile = (double)Math.Min(data.Neach[iPoly], data.Ntotal - data.Neach[iPoly]) / (double)this.Classifier.Instance.Quant[iPoly].Nquantiles;
					this.Classifier.Instance.Quant[iPoly].Pmin = 1.0 / nPerQuantile;
					this.Classifier.Instance.Quant[iPoly].Pmax = 1.0-this.Classifier.Instance.Quant[iPoly].Pmin;

					//	Evaluate the polynomial expression for each datum.
					iDatum = 0;
					for (int iCat = 0; iCat < data.Ncats; iCat++)
					{
						int nSamp = data.Neach[iCat];
						for (int iSamp = 0; iSamp < nSamp; iSamp++)
						{
							//	Evaluate the polynomial expression.
							yVec[iPoly][iDatum++] = (float)this.Classifier.Instance.EvalPolyFromExpanded(iPoly, data.X[iCat][iSamp]);
						}
					}

					//	Sort the output.  Indexes are preserved to speed up subsequent sorts.
					F1.Util.QuickSortIndex(idxVec[iPoly], yVec[iPoly], 0, data.Ntotal-1);

					//	Quantize the output.
					this.Classifier.Instance.Quant[iPoly].Measure(idxVec[iPoly], yVec[iPoly], catVec, (byte)iPoly, this.CatWeights, totalWeight, regressor);
				}
								
				//-----------------------
				//	Measure the conditional entropy.
				//-----------------------
				float[] y = new float[this.Classifier.Instance.Npoly];
				double[] p;
				double fit = 0.0;
				byte c;
				for(iDatum=0; iDatum < data.Ntotal; iDatum++)
				{
					for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
						y[iPoly] = yVec[iPoly][iDatum];
					p = this.Classifier.Instance.ClassifyPolynomialOutputs(y);
					c = catVec[iDatum];
					fit += Math.Log(p[c]) * this.CatWeights[c];
				}
				//	Change logarithm base and normalize by total weight.
				this.Classifier.Fit = -fit / totalWeight / Math.Log((double)data.Ncats);	//	<-- The conditional entropy... we want to minimize it.
				
				//-----------------------
				//	Prepare optimization memory.
				//-----------------------
				
				//	Initialize an orthonormal if the parameter space is small enough.
				float[][][] ortho = null;
				int nParams = this.Classifier.Instance.Npoly * this.Classifier.Instance.Coeffs.Ncoeffs;
				if (nParams <= 100)
					ortho = this.randomDeviates();
				int iOrtho = 0;
				int ctOrtho = 0;
				int ctSteps = 0;
				
				//	Every `modOrtho` steps through the orthonormal basis, we try a gradient search.
				int modOrtho = Math.Min(10, nParams);

				//	For a trip through `modOrtho` bases, changes in entropy are partialled across the parameter space.
				float[][] dhOrtho = Util.NewArrays<float>(this.Classifier.Instance.Npoly, this.Classifier.Instance.Coeffs.Ncoeffs);

				//-----------------------
				//	Optimize.
				//-----------------------
				//	The attempted classifier.
				Classifier cTry = this.Classifier.Instance.Copy();
				//	The initial step size.
				float stepSize = this.Options.ParamDiffMax;
				//	The optimization mode.  We start by iterating through the orthogonal bases.
				OptimizationMode mode = OptimizationMode.Ortho;

				throw new NotImplementedException("TO DO");
				
				bool keepOptimizing = true;
				while (keepOptimizing)
				{
					if(mode==OptimizationMode.Ortho)
					{
						if(iOrtho >= modOrtho)
						{
							iOrtho = 0;
						}
					}
					else if(mode==OptimizationMode.Gradient)
					{
					}
					else
						throw new ApplicationException("Unhandled optimization mode.");
				}
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
		/// this order is calculated by calling <see cref="Util.QuickSortIndex"/> to sort yVals[yIdx[:]].
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
			F1.Util.QuickSortIndex(yIdx, yVals, 0, yVals.Length);

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
			double[,] ortho = D.Util.RandomRotationMatrix(nParams);
			float[][][] output = Util.NewArrays<float>(nParams, this.Classifier.Instance.Npoly, this.Classifier.Instance.Coeffs.Ncoeffs);
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
			float[][] output = Util.NewArrays<float>(this.Classifier.Instance.Npoly, this.Classifier.Instance.Coeffs.Ncoeffs);
			double sumsq = 0.0;
			double r;
			for (int iPoly = 0; iPoly < this.Classifier.Instance.Npoly; iPoly++)
			{
				for (int iCoeff = 0; iCoeff < this.Classifier.Instance.Coeffs.Ncoeffs; iCoeff++)
				{
					r = D1.GaussianDistribution.InvCdf(Util.Rand.NextDouble());
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
