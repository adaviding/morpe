using System;
using System.Threading;

namespace Morpe
{
	/// <summary>
	/// This represents the output of a quantization operation where decision values for an entire training sample are quantized
	/// to generate a quantized probability.  Decision values are calculated for each training datum using one of the polynomial functions.
	/// This instance refers to an entire training sample, but it is only relevant for one of the polynomial functions.
	/// </summary>
	public class Quantization
	{
		/// <summary>
		/// The number of quantiles.
		/// </summary>
		public readonly int Nquantiles;
		/// <summary>
		/// The quantized probability.
		/// </summary>
		public readonly double[] P;
		/// <summary>
		/// The minimum probability limit.  This is computed as a function of the average weight per bin.
		/// </summary>
		public double Pmin = 0.0;
		/// <summary>
		/// The maximum probability limit.  This is computed as a function of the average weight per bin.
		/// </summary>
		public double Pmax = 1.0;
		/// <summary>
		/// The averate Y-value for data inside each quantile.
		/// </summary>
		public readonly double[] Ymid;
		/// <summary>
		/// The boundaries that separate the quantile.  This array has 1 fewer element than <see cref="Ymid"/>.
		/// </summary>
		public readonly double[] Ysep;
		/// <summary>
		/// Constructs a new container for quantized data.
		/// </summary>
		/// <param name="Nquantiles">The number of quantiles.</param>
		public Quantization(int Nquantiles)
		{
			this.Nquantiles = Nquantiles;
			this.Ymid = new double[Nquantiles];
			this.P = new double[Nquantiles];
			this.Ysep = new double[Nquantiles - 1];
		}
		/// <summary>
		/// Copies the quantization data.
		/// </summary>
		/// <returns></returns>
		public Quantization Copy()
		{
			Quantization output = new Quantization(this.Nquantiles);
			output.Pmin = this.Pmin;
			output.Pmax = this.Pmax;
			for (int i = 0; i < this.Nquantiles; i++)
			{
				output.P[i] = this.P[i];
				output.Ymid[i] = this.Ymid[i];
				if(i < this.Nquantiles-1)
					output.Ysep[i] = this.Ysep[i];
			}
			return output;
		}
		/// <summary>
		/// Begins measuring the quantiles using an intermediate result calculated from training data.  This method is
		/// called repeatedly during classifier optimization.
		/// </summary>
		/// <param name="yIdx">[iDatum] Contains unique zero-based indices into yValues such that yValues[yIdx[iDatum]]
		/// increases as iDatum increases.</param>
		/// <param name="yValues">[iDatum] The y-value calculated for each datum.</param>
		/// <param name="cat">[iDatum] The category label of each datum.</param>
		/// <param name="targetCat">[iDatum] The target category of each datum.</param>
		/// <param name="catWeight">The weight assigned to each category label.</param>
		/// <param name="totalWeight">The total weight for the entire sample.</param>
		public void Measure(int[] yIdx, float[] yValues, byte[] cat, byte targetCat, double[] catWeight, double totalWeight,
			MonotonicRegressor regressor)
		{
			//	Target weight per bin.
			double wPerBin = totalWeight / (double)(this.Nquantiles + 0.01);
			//	Keep track of the cumulative weight 
			double wNextBin=wPerBin;
			double w=0.0,dwThisBin;
			double wLastDatum=0.0, wLastBin=0.0, dwThis=0.0, wcBin=0.0, yBin=0.0;
			bool doRewind = false;
			//	Keep track of the current bin number.
			int iBin = 0;
			for(int iDatum=0; iDatum<yValues.Length; iDatum++)
			{
				int iiDatum = yIdx[iDatum];
				byte c = cat[iiDatum];
				//	Update the weight.
				dwThis = catWeight[c];
				wLastDatum = w;
				w += dwThis;
				//	Is the bin finished?
				if( w>= wNextBin )
				{
					//	Is the current datum closest to the bin boundary?  Or the last datum?
					doRewind = wNextBin-wLastDatum < w-wNextBin;
					if(doRewind)
					{
						//	The last datum is closer.  Rewind.
						iiDatum = yIdx[--iDatum];
						w = wLastDatum;
					}
					else
					{
						//	The current datum is closer.  Accumulate totals.
						yBin += dwThis * yValues[iiDatum];
						if (c==targetCat) wcBin += dwThis;
					}
					//	Finalize the current bin.
					dwThisBin = w - wLastBin;
					this.P[iBin] = wcBin / dwThisBin;
					this.Ymid[iBin] = yBin / dwThisBin;
					if(iBin<this.Ysep.Length)
					{
						float ysep = yValues[iiDatum];
						if(doRewind)
							ysep = (ysep + yValues[yIdx[iDatum+1]]) / 2.0f;
						else if(iDatum>0)
							ysep = (ysep + yValues[yIdx[iDatum-1]]) / 2.0f;
						this.Ysep[iBin] = ysep;
					}
					//	Advance to the next bin
					iBin++;
					yBin = wcBin = 0.0;
					wLastBin = w;
					wNextBin += wPerBin;
					//	The last bin is special.
					if(iBin==this.Nquantiles-1)
					{
						while(++iDatum<yValues.Length)
						{
							iiDatum = yIdx[iDatum];
							c = cat[iiDatum];
							//	Update the weight.
							dwThis = catWeight[c];
							w += dwThis;
							//	Accumulate totals.
							yBin += dwThis * yValues[iiDatum];
							if (c == targetCat) wcBin += dwThis;
						}
						//	Finalize the last bin.
						dwThisBin = w - wLastBin;
						this.P[iBin] = wcBin / dwThisBin;
						this.Ymid[iBin] = yBin / dwThisBin;
					}
				}
				else
				{
					//	Accumulate totals.
					yBin += dwThis * yValues[iiDatum];
					if (c == targetCat) wcBin += dwThis;
				}
			}
			//	Perform monotonic regression.
			regressor.Run(this.P, (double[])this.P.Clone());
			//	Range limit
			for(iBin=0; iBin<this.P.Length; iBin++)
				this.P[iBin] = Math.Max(this.Pmin, Math.Min(this.Pmax, this.P[iBin]));
		}
	}
}