using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Morpe.Validation;

namespace Morpe
{
    /// <summary>
    /// A trained classifier.
    /// </summary>
    public class TrainedClassifier
    {
        /// <summary>
        /// Selects the best classifier (having the lowest entropy) and updates its fields so that the counts
        /// of steps are the sum of all given classifiers.
        /// </summary>
        /// <param name="args">The trained classifiers.</param>
        /// <returns>The best classifier (having the lowest entropy), with updated counts.</returns>
        [return: NotNull]
        public static TrainedClassifier SelectLowestEntropyAndMerge(params TrainedClassifier[] args)
        {
            return SelectLowestEntropyAndMerge(args);
        }

        /// <summary>
        /// Selects the best classifier (having the lowest entropy) and updates its fields so that the counts
        /// of steps are the sum of all given classifiers.
        /// </summary>
        /// <param name="classifiers">The trained classifiers.</param>
        /// <returns>The best classifier (having the lowest entropy), with updated counts.</returns>
        [return: NotNull]
        public static TrainedClassifier SelectLowestEntropyAndMerge([NotNull] IReadOnlyList<TrainedClassifier> classifiers)
        {
            Chk.NotNull(classifiers, nameof(classifiers));
            Chk.True(classifiers.All(a => a?.Entropy != null), "All classifiers must have an entropy value.");
            Chk.Less(0, classifiers.Count, "At least 1 classifiers must be provided.");

            if (classifiers.Count == 1)
                return classifiers[0];

            TrainedClassifier tc = classifiers[0];
            double maxEntropy = tc.Entropy.Value;
            int iMaxEntropy = 0;

            int numGoodStepsTaken = tc.NumGoodStepsTaken;
            int numStarts = tc.NumAproaches;
            int numStepsTaken = tc.NumStepsTaken;

            for (int i = 1; i < classifiers.Count; i++)
            {
                tc = classifiers[i];

                if (tc.Entropy.Value > maxEntropy)
                {
                    maxEntropy = tc.Entropy.Value;
                    iMaxEntropy = i;
                }
                
                numGoodStepsTaken += tc.NumGoodStepsTaken;
                numStarts += tc.NumAproaches;
                numStepsTaken += tc.NumStepsTaken;
            }

            TrainedClassifier output = classifiers[iMaxEntropy];
            output.NumGoodStepsTaken = numGoodStepsTaken;
            output.NumAproaches = numStarts;
            output.NumStepsTaken = numStepsTaken;

            return output;
        }
        
        /// <summary>
        /// If specified, this gives the classification accuracy.  This is typically measured with respect to the
        /// training data and a <see cref="CategoryWeightingRule"/>.
        /// </summary>
        public double? Accuracy;
        
        /// <summary>
        /// The instance that was trained.
        /// </summary>
        public Classifier Classifier;
        
        /// <summary>
        /// If specified, this gives the conditional entropy.  This is typically measured with respect to the
        /// training data and a <see cref="CategoryWeightingRule"/>.
        /// </summary>
        public double? Entropy;
        
        /// <summary>
        /// An identifier of the trained classifier.
        /// </summary>
        [NotNull]
        public ClassifierId Id { get; private set; }
        
        /// <summary>
        /// This is the number of good steps that were taken through the parameter space for the optimization routine
        /// that yielded the best fit.  (A good step is one where the fit is improved.)
        /// </summary>
        public int NumGoodStepsTaken;
        
        /// <summary>
        /// This is the number of times that the optimization routine was started.  With each optimization routine,
        /// a different path is travelled through parameter space to arrive at the solution.
        /// </summary>
        public int NumAproaches;
        
        /// <summary>
        /// This is the number of steps taken (good and bad) through the parameter space for the optimization routine
        /// that yielded the best fit.
        /// </summary>
        public int NumStepsTaken;

        /// <summary>
        /// The trained classifier.
        /// </summary>
        /// <param name="id"><see cref="Id"/></param>
        public TrainedClassifier([NotNull] ClassifierId id)
            : this(id, classifier: null)
        {
        }
        
        /// <summary>
        /// The trained classifier.
        /// </summary>
        /// <param name="id"><see cref="Id"/></param>
        /// <param name="classifier"><see cref="Classifier"/></param>
        public TrainedClassifier(
            [NotNull]   ClassifierId id,
            [MaybeNull] Classifier classifier)
        {
            Chk.NotNull(id, nameof(id));
            
            this.Id = id;
            this.Classifier = classifier;
        }

        /// <summary>
        /// Increment counters after a bad step is taken in the parameter space.
        /// </summary>
        public void AddBadStep()
        {
            this.NumStepsTaken++;
        }
        
        /// <summary>
        /// Increment counters after a good step is taken in the parameter space.
        /// </summary>
        public void AddGoodStep()
        {
            this.NumGoodStepsTaken++;
            this.NumStepsTaken++;
        }

        /// <summary>
        /// Creates a deep copy.
        /// </summary>
        /// <returns>The deep copy.</returns>
        [return: NotNull]
        public TrainedClassifier Clone()
        {
            TrainedClassifier output = (TrainedClassifier)this.MemberwiseClone();
            output.Classifier = this.Classifier.Clone();
            output.Id = this.Id.Clone();

            return output;
        }
    }
}
