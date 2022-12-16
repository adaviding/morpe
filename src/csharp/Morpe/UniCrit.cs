using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;

namespace Morpe
{
    /// <summary>
    /// Represents a single unidimensional criterion for classification.
    /// </summary>
    public class UniCrit : ICloneable
    {
        /// <summary>
        /// Finds the unidimensional criterion that maximizes accuracy for detecting the target category in a sample.
        /// </summary>
        /// <param name="targetCat">The target category.</param>
        /// <param name="cat">The category labels for each datum in the sample.</param>
        /// <param name="x">The univariate value for each datum in the sample.</param>
        /// <param name="idx">The index into x[.] and cat[.] which rank-orders the sample in terms of increasing x.</param>
        /// <param name="catWeight">The weight of each category on the accuracy.</param>
        /// <returns>The unidimensional criterion.</returns>
        [return: NotNull]
        public static UniCrit MaximumAccuracy(
            byte targetCat,
            [NotNull] byte[] cat,
            [NotNull] float[] x,
            [NotNull] int[] idx,
            [NotNull] double[] catWeight)
        {
            UniCrit output = new UniCrit();
            output.MaximizeAccuracy(targetCat, cat, x, idx, catWeight);
            return output;
        }
        
        /// <summary>
        /// The classification accuracy.
        /// </summary>
        public double Accuracy = float.NaN;
        
        /// <summary>
        /// The index where the criterion is placed (with respect to a training sample).
        /// Typically, this index is halfway between two integers.
        /// </summary>
        public float CritIndex = float.NaN;
        
        /// <summary>
        /// The value that determines the placement of a univariate criterion.
        /// </summary>
        public float CritValue = float.NaN;
        
        /// <summary>
        /// If true, only the values above CritValue are assigned the label of the target category; otherwise, if false,
        /// only the values below CritValue are assigned the label of the target category.
        /// </summary>
        public bool TargetUpper = true;
        
        /// <summary>
        /// Creates a deep copy.  See <see cref="ICloneable.Clone"/>.
        /// </summary>
        /// <returns>A deep copy.</returns>
        [return: NotNull]
        public object Clone()
        {
            UniCrit output = (UniCrit)this.MemberwiseClone();
            return output;
        }
        
        /// <summary>
        /// Optimizes this criterion to maximize accuracy for detecting the target category in a sample.
        /// </summary>
        /// <param name="targetCat">The target category.</param>
        /// <param name="cat">The category labels for each datum in the sample.</param>
        /// <param name="x">The univariate value for each datum in the sample.</param>
        /// <param name="idx">The index into x[.] and cat[.] which rank-orders the sample in terms of increasing x.</param>
        /// <param name="catWeight">The weight of each category on the accuracy.</param>
        public void MaximizeAccuracy(
            byte targetCat,
            [NotNull] byte[] cat,
            [NotNull] float[] x,
            [NotNull] int[] idx,
            [NotNull] double[] catWeight)
        {
            //double cwTotal = Static.Sum(catWeight);
            int i, idCat;
            double w;
            int nSamp = cat.Length;

            //--------------------------------------------------------
            //    Compute the total weight for the target-category and the non-target-category samples.
            //--------------------------------------------------------
            double wt0 = 0.0; // non-target weight
            double wt1 = 0.0; //     target weight
            for (i = 0; i < nSamp; i++)
            {
                idCat = cat[idx[i]];
                if (idCat == targetCat)
                    wt1 += catWeight[idCat];
                else
                    wt0 += catWeight[idCat];
            }
            double wt = wt0 + wt1;    // total weight

            //--------------------------------------------------------
            //    Simulate the criteria by stepping through x[idx[.]].
            //    For each criterion value, two decision rules are actually simulated:
            //        1.  where TargetUpper=true
            //        2.  where TargetUpper=false
            //    Maximize accuracy by maximizing the difference between the correctly and incorrectly classified weight.
            //--------------------------------------------------------

            //    Initialize variables to i == -1
            double wPos = wt1;            // The correctly classified weight for TargetUpper=true
            double wNeg = wt0;            // The correctly classified weight for TargetUpper=false

            double wPosMax = wPos;        // The maximum of wPos over all i
            int ilPos = -1;                // The lowest value of i for which wPos==wPosMax.
            int cPos = 0;                // The number of values of i for which wPos==wPosMax.

            double wNegMax = wNeg;        // The maximum of wNeg over all i
            int ilNeg = -1;                // The lowest value of i for which wNeg==wNegMax.
            int cNeg = 0;                // The number of values of i for which wNeg==wNegMax.

            //    Try shifting the criterion at each value of the sample (for samples in ascending order).
            for (i = 0; i < nSamp; i++)
            {
                idCat = cat[idx[i]];
                w = catWeight[idCat];
                //    Increment or decrement the correctly classified weight for each decision rule.
                if (idCat == targetCat)
                {
                    //    Positive rule declines
                    wPos -= w;
                    //    Negative rule improves
                    wNeg += w;
                }
                else
                {
                    //    Positive rule improves
                    wPos += w;
                    //    Negative rule declines
                    wNeg -= w;
                }
                //    Track the best-performing positive and negative rules.
                //    If this Z-value is the same as the last Z-value, then a criterion cannot be placed at that i, so the rule does not get tracked.
                if (i == 0 || x[idx[i]] != x[idx[i - 1]])
                {
                    //    A criterion can be placed here because this Z-value increases.
                    if (wPos > wPosMax)
                    {
                        wPosMax = wPos;
                        ilPos = i;
                        cPos = 1;
                    }
                    else if (wPos == wPosMax)
                    {
                        cPos++;
                    }
                    if (wNeg > wNegMax)
                    {
                        wNegMax = wNeg;
                        ilNeg = i;
                        cNeg = 1;
                    }
                    else if (wNeg == wNegMax)
                    {
                        cNeg++;
                    }
                }
            }

            //--------------------------------------------------------
            //    Recover the output and criterion values.
            //--------------------------------------------------------
            if( wPosMax >= wNegMax )
            {
                this.TargetUpper = true;
                //    Compute the accuracy
                this.Accuracy = (float)(wPosMax/wt);
                if (wPos >= wPosMax)
                {
                    this.CritIndex = (float)nSamp - 0.5f;
                    this.CritValue = x[idx[nSamp - 1]] + 1.0f;
                }
                else if (cPos <= 2 || ilPos == -1)
                {
                    this.CritIndex = (float)ilPos + 0.5f;
                    if (ilPos == -1)
                        this.CritValue = x[idx[0]] - 1.0f;
                    else if(ilPos==nSamp-1)
                        this.CritValue = x[idx[nSamp - 1]] + 1.0f;
                    else
                        this.CritValue = (x[idx[ilPos]] + x[idx[ilPos+1]])/2.0f;
                }
                else
                {
                    //------------------------------------------------------------
                    //    Search for the cPos/2 maximum somewhere after ilPos
                    //------------------------------------------------------------
                    cPos /= 2;
                    wPos = wPosMax;
                    i = ilPos;
                    while (cPos > 0)
                    {
                        idCat = cat[idx[i++]];
                        w = catWeight[idCat];
                        if (idCat == targetCat)
                            wPos -= w;
                        else
                            wPos += w;
                        //    Assert
                        if (i == nSamp)
                        {
                            //    Assertion failed.  This should never happen.
                            i = ilPos;
                            //    Exit loop
                            cPos = 0;
                        }
                        //    If this Z-value is the same as the last Z-value, then a criterion cannot be placed at that i.
                        if (i == 0 || x[idx[i]] != x[idx[i - 1]])
                        {
                            //    A criterion can be placed here because this Z-value increases.
                            if (wPos >= wPosMax)
                                cPos--;
                        }
                    }
                    this.CritIndex = (float)i + 0.5f;
                    this.CritValue = (x[idx[i]] + x[idx[i + 1]]) / 2.0f;
                    //------------------------------------------------------------
                }
            }
            else
            {
                this.TargetUpper = false;
                this.Accuracy = wNegMax/wt;
                if( wNeg >= wNegMax )
                {
                    this.CritIndex = (float)nSamp - 0.5f;
                    this.CritValue = x[idx[nSamp - 1]] + 1.0f;
                }
                else if( ilNeg==-1 )
                {
                    this.CritIndex = -0.5f;
                    this.CritValue = x[idx[0]]-1.0f;
                }
                else if (cNeg <= 2)
                {
                    this.CritIndex = ilNeg + 0.5f;
                    if(ilNeg==nSamp-1)
                        this.CritValue = x[idx[nSamp - 1]] + 1.0f;
                    else
                        this.CritValue = (x[idx[ilNeg]] + x[idx[ilNeg + 1]]) / 2.0f;
                }
                else
                {
                    //------------------------------------------------------------
                    //    Search for the cNeg/2 minimum somewhere after ilNeg
                    //------------------------------------------------------------
                    cNeg /= 2;
                    wNeg = wNegMax;
                    i = ilNeg;
                    while (cNeg > 0)
                    {
                        idCat = cat[idx[i++]];
                        w = catWeight[idCat];
                        if (idCat == targetCat)
                            wNeg += w;
                        else
                            wNeg -= w;
                        //    Assert
                        if (i == nSamp)
                        {
                            //    Assertion failed.  This should never happen.
                            i = ilNeg;
                            //    Exit loop
                            cNeg = 0;
                        }
                        //    If this Z-value is the same as the last Z-value, then a criterion cannot be placed at that i.
                        if (i == 0 || x[idx[i]] != x[idx[i - 1]])
                        {
                            //    A criterion can be placed here because this Z-value increases.
                            if (wNeg >= wNegMax)
                                cNeg--;
                        }
                    }
                    this.CritIndex = (float)i + 0.5f;
                    this.CritValue = (x[idx[i]] + x[idx[i + 1]]) / 2.0f;
                    //------------------------------------------------------------
                }
            }
        }
    }
}
