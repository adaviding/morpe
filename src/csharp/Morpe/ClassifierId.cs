using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Morpe
{
    /// <summary>
    /// Identifies a classifier based on polynomial rank and dimensionality, so that its performance may be compared
    /// to other classifiers having slightly different attributes.
    /// </summary>
    public class ClassifierId
    {
        /// <summary>
        /// The spatial dimensions included by the classifier, listed in ascending order without duplicates.
        /// </summary>
        public int[] Dims;

        /// <summary>
        /// The rank of the polynomial used by the classifier.
        /// </summary>
        public int Rank;

        /// <summary>
        /// If this is a "dual" classifier, then this is the zero-based index of the target category; otherwise this is
        /// null.
        ///
        /// Dual classifiers are only used when the number of categories is greater than 2.  
        /// </summary>
        public int? TargetCat;

        /// <summary>
        /// Constructs an instance having the specified properties.
        /// </summary>
        /// <param name="dims"><see cref="Dims"/></param>
        /// <param name="rank"><see cref="Rank"/></param>
        /// <param name="targetCat"><see cref="TargetCat"/></param>
        public ClassifierId(
            [NotNull] int[] dims,
            int rank,
            int? targetCat)
        {
            this.Dims = dims.ToArray();
            this.Rank = rank;
            this.TargetCat = targetCat;
        }

        [return: NotNull]
        public ClassifierId Clone()
        {
            return new ClassifierId(
                dims: this.Dims.ToArray(),
                rank: this.Rank,
                targetCat: this.TargetCat);
        }
        
        [return: MaybeNull]
        public ClassifierId NextLowerSubrank()
        {
            ClassifierId output = null;
            
            if (this.Rank > 1)
            {
                output = this.Clone();
                output.Rank--;
            }

            return output;
        }
    }
}