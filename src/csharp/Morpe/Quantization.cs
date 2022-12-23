using System;
using System.Diagnostics.CodeAnalysis;
using System.Threading;

using D1 = Morpe.Numerics.D1;

namespace Morpe
{
    /// <summary>
    /// This represents the output of a quantization operation where decision values for an entire training sample are quantized
    /// to generate a quantized probability.  Decision values are calculated for each training datum using one of the polynomial functions.
    /// This instance refers to an entire training sample, but it is only relevant for one of the polynomial functions.
    /// </summary>
    public class Quantization : ICloneable
    {
        /// <summary>
        /// The number of quantiles.
        /// </summary>
        public int NumQuantiles { get; private set; }
        
        /// <summary>
        /// The quantized probability.
        /// </summary>
        public double[] P  { get; private set; }

        /// <summary>
        /// The range of P-values that can be delivered by this quantization.
        /// </summary>
        public D1.Range ProbabilityRange { get; private set; } 
        
        /// <summary>
        /// The average Y-value for data inside each quantile.
        /// </summary>
        public double[] Ymid  { get; private set; }
        
        /// <summary>
        /// The boundaries that separate the quantile.  This array has 1 fewer element than <see cref="Ymid"/>.
        /// </summary>
        public double[] Ysep  { get; private set; }
        
        /// <summary>
        /// Constructs a new container for quantized data.
        /// </summary>
        /// <param name="numQuantiles"><see cref="NumQuantiles"/></param>
        public Quantization(int numQuantiles)
        {
            this.NumQuantiles = numQuantiles;
            this.Ymid = new double[numQuantiles];
            this.P = new double[numQuantiles];
            this.ProbabilityRange = new D1.Range(0.0, 1.0);
            this.Ysep = new double[numQuantiles - 1];
        }
        
        /// <summary>
        /// Constructs a new container for quantized data.
        /// </summary>
        /// <param name="numQuantiles"><see cref="NumQuantiles"/></param>
        /// <param name="probabilityRange"><see cref="ProbabilityRange"/></param>
        public Quantization(
            int numQuantiles,
            [NotNull] D1.Range probabilityRange)
        {
            this.NumQuantiles = numQuantiles;
            this.Ymid = new double[numQuantiles];
            this.P = new double[numQuantiles];
            this.ProbabilityRange = probabilityRange.Clone();
            this.Ysep = new double[numQuantiles - 1];
        }
        
        /// <summary>
        /// Creates a deep copy of the quantization data.
        /// </summary>
        /// <returns>The deep copy.</returns>
        public Quantization Clone()
        {
            Quantization output = (Quantization)this.MemberwiseClone();
            output.P = (double[])this.P.Clone();
            output.Ymid = (double[])this.Ymid.Clone();
            output.Ysep = (double[])this.Ysep.Clone();
            
            return output;
        }
        object ICloneable.Clone()
        {
            return this.Clone();
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
        /// <param name="catWeights">The weight assigned to each category label.</param>
        /// <param name="regressor">If specified, monotonic regression is performed with this object; otherwise it is
        /// not performed.</param>
        public void Measure(
            CancellationToken cancellationToken,
            [NotNull] int[] yIdx,
            [NotNull] float[] yValues,
            [NotNull] byte[] cat,
            byte targetCat,
            [NotNull] CategoryWeights catWeights,
            [MaybeNull] D1.MonotonicRegressor regressor)
        {
            //    Target weight per bin.
            double wPerBin = catWeights.TotalWeight / (double)(this.NumQuantiles + 0.01);
            //    Keep track of the cumulative weight 
            double wNextBin=wPerBin;
            double w=0.0,dwThisBin;
            double wLastDatum=0.0, wLastBin=0.0, dwThis=0.0, wcBin=0.0, yBin=0.0;
            bool doRewind = false;
            //    Keep track of the current bin number.
            int iBin = 0;
            for(int iDatum=0; iDatum<yValues.Length; iDatum++)
            {
                // Check for cancellation.
                cancellationToken.ThrowIfCancellationRequested();
                
                int iiDatum = yIdx[iDatum];
                byte c = cat[iiDatum];
                
                //    Update the weight.
                dwThis = catWeights.Weights[c];
                wLastDatum = w;
                w += dwThis;
                
                //    Is the bin finished?
                if( w>= wNextBin )
                {
                    //    Is the current datum closest to the bin boundary?  Or the last datum?
                    doRewind = wNextBin-wLastDatum < w-wNextBin;
                    if(doRewind)
                    {
                        //    The last datum is closer.  Rewind.
                        iiDatum = yIdx[--iDatum];
                        w = wLastDatum;
                    }
                    else
                    {
                        //    The current datum is closer.  Accumulate totals.
                        yBin += dwThis * yValues[iiDatum];
                        if (c==targetCat) wcBin += dwThis;
                    }
                    //    Finalize the current bin.
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
                    //    Advance to the next bin
                    iBin++;
                    yBin = wcBin = 0.0;
                    wLastBin = w;
                    wNextBin += wPerBin;
                    //    The last bin is special.
                    if(iBin==this.NumQuantiles-1)
                    {
                        while(++iDatum<yValues.Length)
                        {
                            iiDatum = yIdx[iDatum];
                            c = cat[iiDatum];
                            //    Update the weight.
                            dwThis = catWeights.Weights[c];
                            w += dwThis;
                            //    Accumulate totals.
                            yBin += dwThis * yValues[iiDatum];
                            if (c == targetCat) wcBin += dwThis;
                        }
                        //    Finalize the last bin.
                        dwThisBin = w - wLastBin;
                        this.P[iBin] = wcBin / dwThisBin;
                        this.Ymid[iBin] = yBin / dwThisBin;
                    }
                }
                else
                {
                    //    Accumulate totals.
                    yBin += dwThis * yValues[iiDatum];
                    if (c == targetCat) wcBin += dwThis;
                }
            }
            
            //    Perform monotonic regression.
            if (regressor != null)
                regressor.Run(cancellationToken, this.P, (double[])this.P.Clone());
            
            //    Range limit
            for(iBin=0; iBin<this.P.Length; iBin++)
                this.P[iBin] = this.ProbabilityRange.Clamp(this.P[iBin]);
        }
    }
}
