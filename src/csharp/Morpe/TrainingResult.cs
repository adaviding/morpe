namespace Morpe
{
    /// <summary>
    /// The result of classifier training.
    /// </summary>
    public class TrainingResult
    {
        /// <summary>
        /// The instance that did the training.
        /// </summary>
        public Trainer Trainer;
        
        /// <summary>
        /// The type of training result.
        /// </summary>
        public TrainingResultType Type;

        /// <summary>
        /// The number of steps which resulted in an improvement in fit (i.e. lower entropy).
        /// </summary>
        public int GoodStepsTaken = 0;

        /// <summary>
        /// The number of steps taken during gradient descent which resulted in an improvement in fit (i.e. lower entropy).
        /// </summary>
        public int GoodStepsTakenFromGradientDescent = 0;
        
        /// <summary>
        /// The number of steps taken during gradient search which resulted in an improvement in fit (i.e. lower entropy).
        /// Gradient search is when we try different orthonormal bases in search of a gradient.
        /// </summary>
        public int GoodStepsTakenFromGradientSearch = 0;
        
        /// <summary>
        /// The total number of steps taken, good and bad.
        /// </summary>
        public int StepsTaken;
    }
}
