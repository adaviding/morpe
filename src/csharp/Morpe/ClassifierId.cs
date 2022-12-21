using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Runtime.InteropServices;
using Morpe.Validation;

namespace Morpe
{
    /// <summary>
    /// Identifies a classifier based on polynomial rank and dimensionality, so that its performance may be compared
    /// to other classifiers having slightly different attributes.
    ///
    /// This can be used as a dictionary key.
    /// </summary>
    public class ClassifierId : IEquatable<ClassifierId>
    {
        /// <summary>
        /// The spatial dimensions included by the classifier, listed in ascending order without duplicates.
        /// </summary>
        public IReadOnlyList<int> Dims => this.dims;

        /// <summary>
        /// The rank of the polynomial used by the classifier.
        /// </summary>
        public int Rank { get; private set; }

        /// <summary>
        /// If this is a "dual" classifier, then this is the zero-based index of the target category; otherwise this is
        /// null.
        ///
        /// Dual classifiers are only used when the number of categories is greater than 2.  
        /// </summary>
        public int? TargetCat { get; private set; }

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
            Chk.NotNull(dims, nameof(dims));

            this.dims = (int[])dims.Clone();
            this.Rank = rank;
            this.TargetCat = targetCat;
        }

        [return: NotNull]
        public ClassifierId Clone()
        {
            ClassifierId output = (ClassifierId)this.MemberwiseClone();
            output.dims = this.dims.ToArray();
            
            return output;
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

        /// <summary>
        /// Compares two instances by value.
        /// </summary>
        /// <param name="that">The other instance.</param>
        /// <returns>Returns true if this instance is equal to the other one, false otherwise.</returns>
        public bool Equals(ClassifierId that)
        {
            if (that == null)
                return false;
            
            if (that.Rank != this.Rank)
                return false;
            
            if (that.TargetCat != this.TargetCat)
                return false;

            if (!Test.EqualListing(that.Dims, this.Dims))
                return false;

            return true;
        }

        /// <summary>
        /// Compares two instances by value.
        /// </summary>
        /// <param name="that">The other instance.</param>
        /// <returns>Returns true if this instance is equal to the other one, false otherwise.</returns>
        public override bool Equals(object that)
        {
            if (that != null && that is ClassifierId cid)
                return Equals(cid);

            return false;
        }

        /// <summary>
        /// Gets the hash code.  The value is cached for better performance.
        /// </summary>
        /// <returns>The hash code.</returns>
        public override int GetHashCode()
        {
            if (!this.hashCode.HasValue)
            {
                HashCode hasher = new HashCode();
                hasher.Add(this.dims.Length);
                hasher.Add(this.Rank);
                hasher.Add(this.TargetCat);
                foreach (int dim in this.dims)
                {
                    hasher.Add(dim);
                }

                this.hashCode = hasher.ToHashCode();
            }

            return this.hashCode.Value;
        }
        
        /// <summary>
        /// The spatial dimensions included by the classifier, listed in ascending order without duplicates.
        /// </summary>
        private int[] dims;

        /// <summary>
        /// The cached hash code.
        /// </summary>
        private int? hashCode;
    }
}
