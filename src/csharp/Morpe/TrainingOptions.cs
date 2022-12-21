using System;

namespace Morpe
{
    /// <summary>
    /// Training options.
    /// </summary>
    public class TrainingOptions
    {
        /// <summary>
        /// This is called to construct a new instance of the <see cref="Morpe.Options"/> class.
        /// </summary>
        public TrainingOptions() { }

        /// <summary>
        /// If specified, limits the concurrency during training.
        /// </summary>
        public int? ConcurrencyLimit;
        
        /// <summary>
        /// Optimization halts when the change in conditional entropy is less than this amount.  This number must be
        /// positive.  Default = 0.005.
        /// </summary>
        public float EntropyTol = 0.005f;

        /// <summary>
        /// If true, the classifier data is conditioned and then expanded "in place" to save memory.  This alters the
        /// input training data, but it also saves memory.
        ///
        /// If false, new memory is allocated for storing conditioned and expanded data.
        /// </summary>
        public bool ExpandDataInPlace = true;
        
        /// <summary>
        /// If true, the solver will pick initial values for the classifier's parameters; otherwise the classifier's
        /// existing parameters will be used.
        /// </summary>
        public bool InitializeParams = true;
        
        /// <summary>
        /// The optimization routine will run this many times.  Each time, it will approach the solution along a
        /// different path.
        /// </summary>
        public int NumberOfApproaches = 1;
        
        /// <summary>
        ///    The rate at which parameters shrink during optimization.  Must be greater than 1.0.
        /// </summary>
        public float ParamShrinkFactor = 2.0f;
        
        /// <summary>
        /// Must be greater than 0.  Optimization cannot change parameters less than the following magnitude for any [p,i]-th parameter.
        ///        [MAGNITUDE] = ParamDiffMax * Morpe.Trainer.ParamScale[p][i]
        /// </summary>
        public float ParamDiffMax = 2.0f;
        
        /// <summary>
        /// Must be greater than 0.  Optimization cannot change parameters less than the following magnitude for any [p,i]-th parameter.
        ///        [MAGNITUDE] = ParamDiffMin * Morpe.Trainer.ParamScale[p][i]
        /// </summary>
        public float ParamDiffMin = 0.005f;
        
        /// <summary>
        /// Determines how the optimization routine weights data from each category.
        /// </summary>
        public CategoryWeightingRule CategoryWeightingRule = CategoryWeightingRule.EqualPriors;
        
        /// <summary>
        /// Creates a deep copy of the instance.
        /// </summary>
        /// <returns>The deep copy.</returns>
        public TrainingOptions Clone()
        {
            TrainingOptions output = (TrainingOptions)this.MemberwiseClone();
            return output;
        }
    }
}
