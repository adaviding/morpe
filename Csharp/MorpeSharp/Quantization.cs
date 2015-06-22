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
	}
}