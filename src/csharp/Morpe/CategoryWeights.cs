using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Morpe.Validation;

namespace Morpe
{
    /// <summary>
    /// This allows the <see cref="Trainer"/> to assign proper weight to each training datum.
    /// </summary>
    public class CategoryWeights
    {
        /// <summary>
        /// The weight assignable to a single training datum from each category. 
        /// </summary>
        [NotNull]
        public IReadOnlyList<double> Weights => weights;
        private double[] weights;
        
        /// <summary>
        /// The total weight contained in the training data.
        /// </summary>
        public double TotalWeight { get; private set; }

        internal CategoryWeights(int numCats)
        {
            Chk.Less(1, numCats, "The number of categories must be 2 or greater.");
            
            this.weights = new double[numCats];
        }

        /// <summary>
        /// Constructs an instance which is appropriate for the given information.
        ///
        /// This can be used for regular classifiers and "dual" classifiers.  "Dual" classifiers are constructed for
        /// classification problems having M categories, where M > 2.  Each dual classifier is identified by its target
        /// category.  It is a 2-category classifier that segregates the target category (as category 0) and the sum of
        /// all other categories (as category 1).
        /// </summary>
        /// <param name="numEach">The number of training data in each category.</param>
        /// <param name="rule">The category weighting rule.</param>
        /// <param name="targetCategory">If this measurement is for a "dual" classifier, then this is the zero-based
        /// index of its target category; otherwise it is null.</param>
        /// <returns>The new instance.</returns>
        public static CategoryWeights Measure(
            [NotNull] int[] numEach,
            CategoryWeightingRule rule,
            int? targetCategory)
        {
            Chk.NotNull(numEach, nameof(numEach));

            int numCats = numEach.Length;
            int numTotal = numEach.Sum();
            
            CategoryWeights output = new CategoryWeights(numCats);

            if (rule == CategoryWeightingRule.EqualPriors)
            {
                for (int iCat = 0; iCat < numCats; iCat++)
                {
                    double weight = (double)numTotal / numEach[iCat] / numCats;

                    if (targetCategory != null)
                    {
                        if (iCat == targetCategory)
                            weight *= numCats / 2.0;
                        else
                            weight *= numCats / (2.0 * (numCats - 1));
                    }

                    output.weights[iCat] = weight;
                    output.TotalWeight += weight * numEach[iCat];
                }
            }
            else if (rule == CategoryWeightingRule.ObservedPriors)
            {
                for (int iCat = 0; iCat < numCats; iCat++)
                {
                    output.weights[iCat] = 1.0;
                    output.TotalWeight += numEach[iCat];
                }
            }
            else
            {
                throw new NotImplementedException($"{nameof(CategoryWeightingRule)}.{rule}");
            }

            return output;
        }

        /// <summary>
        /// Measures the weights for all "dual" classifiers.
        /// </summary>
        /// <param name="numEach">The number of training data in each category.</param>
        /// <param name="rule">The category weighting rule.</param>
        /// <returns>The new instance.</returns>
        public static CategoryWeights[] MeasureAllDuals(
            [NotNull] int[] numEach,
            CategoryWeightingRule rule)
        {
            CategoryWeights[] output = new CategoryWeights[numEach.Length];

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = CategoryWeights.Measure(numEach, rule, i);
            }

            return output;
        }
    }
}
