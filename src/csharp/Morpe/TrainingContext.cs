using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Morpe.Validation;

namespace Morpe
{
    /// <summary>
    /// This encapsulates knowledge which is relevant to training.
    /// </summary>
    public class TrainingContext
    {
        /// <summary>
        /// The analysis that was conducted before training for the highest polynomial rank being considered.
        /// </summary>
        public PreOptimizationAnalysis Analysis;
        
        /// <summary>
        /// Gaussian stats for the training data.
        /// </summary>
        public GaussianStats GaussianStats;

        /// <summary>
        /// The classifiers that have been trained up to this point.
        /// </summary>
        public IReadOnlyDictionary<ClassifierId, TrainedClassifier> TrainedClassifiers => this.trainedClassifiers;

        /// <summary>
        /// This will add the trained classifier into <see cref="TrainedClassifiers"/>.  If a classifier with the same
        /// ID already exists, then the one having the best entropy will be saved, and the value of
        /// <see cref="TrainedClassifier.NumAproaches"/> will be incremented appropriately.
        /// </summary>
        /// <param name="trainedClassifier">The classifier being submitted.</param>
        public void Add([NotNull] TrainedClassifier trainedClassifier)
        {
            Chk.NotNull(trainedClassifier, nameof(trainedClassifier));
            
            if (this.trainedClassifiers.TryGetValue(trainedClassifier.Id, out TrainedClassifier tcPrior))
            {
                if (tcPrior.Entropy <= trainedClassifier.Entropy)
                {
                    tcPrior.NumAproaches += trainedClassifier.NumAproaches;
                }
                else
                {
                    trainedClassifier.NumAproaches += tcPrior.NumAproaches;
                    this.trainedClassifiers[trainedClassifier.Id] = trainedClassifier;
                }
            }
            else
            {
                this.trainedClassifiers.Add(trainedClassifier.Id, trainedClassifier);
            }
        }

        /// <summary>
        /// Trains a classifier.  This involves training many classifiers which are less complicated than the one
        /// requested.  The context retains all of the trained classifiers so they may be re-used by the caller.
        ///
        /// The caller may obtain the trained classifier from <see cref="TrainedClassifiers"/> by using the given
        /// classifier ID.
        /// </summary>
        /// <param name="data">The training data.</param>
        /// <param name="id">The id of the most complicated classifier to be trained.</param>
        /// <param name="options">The training options.</param>
        public void Train(
            [NotNull] CategorizedData data,
            [NotNull] ClassifierId id,
            [NotNull] TrainingOptions options)
        {
        }

        
        /// <summary>
        /// The classifiers that have been trained up to this point.
        /// </summary>
        private Dictionary<ClassifierId, TrainedClassifier> trainedClassifiers =
            new Dictionary<ClassifierId, TrainedClassifier>();
    }
}
