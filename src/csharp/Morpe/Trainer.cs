using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Text;
using Morpe.Validation;

using D  = Morpe.Numerics.D;
using D1 = Morpe.Numerics.D1;
using F  = Morpe.Numerics.F;
using F1 = Morpe.Numerics.F1;
using I  = Morpe.Numerics.I;
using I1 = Morpe.Numerics.I1;

namespace Morpe
{
    /// <summary>
    /// Trains an instance of the MoRPE classifier using the training data provided.
    /// </summary>
    public class Trainer
    {
        /// <summary>
        /// The default number of quantiles for the trained MoRPE classifier.
        /// </summary>
        public static readonly int DefaultNumQuantiles = 32;
        
        /// <summary>
        /// Trains a single classifier without the <see cref="TrainingContext"/>.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="data">The classification data.</param>
        /// <param name="id">The ID of the classifier to be trained.</param>
        /// <param name="options">The training options.</param>
        /// <param name="analysis">The analysis of classifier data.</param>
        /// <param name="parameterStarts">The starting points of the given classifiers.</param>
        /// <param name="taskScheduler">The task scheduler which limits the concurrency of classifier training.</param>
        /// <returns>The trained classifier.</returns>
        public static TrainedClassifier Train(
            CancellationToken cancellationToken,
            [NotNull] CategorizedData data,
            [NotNull] ClassifierId id,
            [NotNull] TrainingOptions options,
            [NotNull] PreOptimizationAnalysis analysis,
            [NotNull] List<float[][]> parameterStarts,
            [NotNull] TaskScheduler taskScheduler)
        {
            TaskFactory<TrainedClassifier> taskFactory = new TaskFactory<TrainedClassifier>(taskScheduler);
            
            // Start a task for each parameter start.
            Task<TrainedClassifier>[] tasks = new Task<TrainedClassifier>[parameterStarts.Count];

            for (int i = 0; i < parameterStarts.Count; i++)
            {
                tasks[i] = taskFactory.StartNew(
                    () => Train(cancellationToken, data, id, options, analysis, parameterStarts[i]),
                    cancellationToken);
            }

            Task.WaitAll(tasks);

            // If any task did not complete successfully, we will throw an exception.
            if (tasks.Any(a => !a.IsCompletedSuccessfully))
            {
                cancellationToken.ThrowIfCancellationRequested();
                
                // fixme throw appropriate exception

                Task<TrainedClassifier>[] canceled = tasks.Where(a => a.IsCanceled).ToArray();
                Task<TrainedClassifier>[] failed = tasks.Where(a => a.IsFaulted).ToArray();
                Task<TrainedClassifier>[] succeeded = tasks.Where(a => a.IsCompletedSuccessfully).ToArray();

                string message = $"Of {tasks.Length} tasks, {canceled.Length} were canceled, {succeeded.Length} succeeded, and {failed.Length} failed.";

                if (failed.Length == 0)
                {
                    throw new OperationCanceledException(message);
                }

                throw new AggregateException(message, failed.Select(a => a.Exception));
            }

            // At this point we know that all tasks completed successfully.
            TrainedClassifier[] trained = tasks
                .Select(a => a.Result)
                .ToArray();
            
            TrainedClassifier output = TrainedClassifier.SelectLowestEntropyAndMerge(trained);
            return output;
        }

        /// <summary>
        /// Trains a single classifier without the <see cref="TrainingContext"/>.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="data">The classification data.</param>
        /// <param name="id">The ID of the classifier to be trained.</param>
        /// <param name="options">The training options.</param>
        /// <param name="analysis">The analysis of classifier data.</param>
        /// <param name="parameterStart">The starting point for the classifier parameters.</param>
        /// <returns>The trained classifier</returns>
        public static TrainedClassifier Train(
            CancellationToken cancellationToken,
            [NotNull] CategorizedData data,
            [NotNull] ClassifierId id,
            [NotNull] TrainingOptions options,
            [NotNull] PreOptimizationAnalysis analysis,
            [NotNull] float[][] parameterStart)
        {
            Chk.Equal(id.Dims.Count, data.NumDims, "The number of dimensions is specified inconsistently.");
            Chk.NotNull(options, nameof(options));
            Chk.Less(0, options.NumberOfApproaches, "The number of approaches must be a positive integer.");

            // This is the container for output.
            TrainedClassifier output = new TrainedClassifier(id);

            // This will regress the quantization noise onto an arbitrary monotonic function.
            D1.MonotonicRegressor regressor = new D1.MonotonicRegressor();
            
            // Here we use simple logic to pick the number of quantiles.  In the future we may want to use logic
            // which is based on the sample size.
            int numQuantiles = DefaultNumQuantiles;
            
            // The approximate number of samples per quantile.
            double numPerQuantile = (double)data.NumTotal / numQuantiles;
            
            // This is the range of probabilities that the classifier can produce.
            D1.Range pRange = new D1.Range(
                min: 1.0 / numPerQuantile,
                max: 1.0);
            pRange.Max -= pRange.Min;

            // Figure out how many categories we are training.
            int numCats = data.NumCats;
            if (id.TargetCat != null)
            {
                // This means we are training a "dual" classifier where the data from multiple (more than 2) categories
                // is separated into 2 categories:
                //
                // 1.  The data from the target category
                // 2.  The combined data from all other categories.
                numCats = 2;
            }
            
            // Figure out how many polynomial coefficients we are training.
            int numPoly = numCats;
            if (numCats == 2)
            {
                // For the 2-category problem, there is only 1 set of polynomial coefficients.
                numPoly = 1;
            }

            // The polynomial coefficients are defined like this (for each polynomial).
            Poly poly = new Poly(
                numDims: data.NumDims,
                rank: id.Rank);
            
            // This is how many free parameters we are training.
            int numParams = numPoly * poly.NumCoeffs;
            
            // This is the scale of each parameter (for each polynomial).  It is possible that we have scales
            // for higher rank polynomial coefficients, so truncate.
            float[] paramScale = analysis.ParamScale;
            if (paramScale.Length > poly.NumCoeffs)
            {
                // This creates a new array which is a truncated version of the original.
                Array.Resize(ref paramScale, poly.NumCoeffs);
            }
            
            // Figure out how to weight each training datum.
            CategoryWeights weights = CategoryWeights.Measure(
                numEach: data.NumEach,
                rule: options.CategoryWeightingRule,
                targetCategory: id.TargetCat);
            
            CategoryWeights[] dualWeights = null;
            if (numCats > 2)
            {
                dualWeights = CategoryWeights.MeasureAllDuals(
                    numEach: data.NumEach,
                    rule: options.CategoryWeightingRule);
            }
            
            // [data.NumTotal] The category label of each training datum.
            byte[] catVec = data.GetCategoryVector();
            
            // [numPoly][data.NumTotal] The output of the polynomial function for each polynomial, and each training datum.
            float[][] yVec = F.Util.JaggedSquare(numPoly, data.NumTotal);
            
            // [numPoly][data.NumTotal] The index which sorts the value of 'y' for each polynomial.
            int[][] idxVec = I.Util.JaggedSquare(numPoly, data.NumTotal);
            for (int i = 0; i < numPoly; i++)
                I.Util.FillSeries(idxVec[i]);

            if (numParams == 1)
            {
                // The solution to the 1-parameter problem is already known.  We can just set the parameters without
                // searching and then build the classifier.
                
                // The polynomial function that points to the target cat.
                int iTargetPoly = id.TargetCat ?? 0;
                
                // Allocate a classifier and fill it in manually.
                output.Classifier = new Classifier(
                    numCats: numCats,
                    numDims: poly.NumDims,
                    rank: id.Rank,
                    numQuantiles: numQuantiles,
                    probabilityRange: pRange.Clone(),
                    parameters: new float[][] { new float[]
                        {
                            analysis.Crits[iTargetPoly, 0].TargetUpper  // It just needs to point towards the target category.
                                ? +1f
                                : -1f
                        }});

                
                // The same solution would be found every time, so we don't need any more starts.
                output.NumApproaches = Int32.MaxValue;
                output.AddGoodStep();

                // Build the classifier.
                BuildClassifier(
                    cancellationToken: cancellationToken,
                    data: data,
                    regressor: regressor,
                    catWeights: weights,
                    dualCatWeights: dualWeights,
                    catVec: catVec,
                    yVec: yVec,
                    idxVec: idxVec,
                    training: output);

                // Measure the fit.
                (output.Accuracy, output.Entropy) = MeasureFit(
                    cancellationToken: cancellationToken,
                    data: data,
                    catWeights: weights,
                    catVec: catVec,
                    yVec: yVec,
                    training: output);
                Chk.NotNull(output.Entropy, "{0}.{1}", nameof(output), nameof(output.Entropy));

                return output;
            }

            // Get ready to wander through the parameter space.
            ParameterSpaceWanderer wanderer = new ParameterSpaceWanderer(numPoly, paramScale);

            for (int iApproach = 0; iApproach < options.NumberOfApproaches; iApproach++)
            {
                // Check for cancellation.
                cancellationToken.ThrowIfCancellationRequested();
                
                // Initialize the classifier.
                TrainedClassifier approach = new TrainedClassifier(
                    id: id,
                    classifier: new Classifier(
                        numCats,
                        data.NumDims,
                        id.Rank,
                        numQuantiles,
                        probabilityRange: pRange.Clone(),
                        parameters: Util.Clone(parameterStart)));
                
                // Build the classifier.
                BuildClassifier(
                    cancellationToken: cancellationToken,
                    data: data,
                    regressor: regressor,
                    catWeights: weights,
                    dualCatWeights: dualWeights,
                    catVec: catVec,
                    yVec: yVec,
                    idxVec: idxVec,
                    training: approach);
                
                // Measure the fit.
                (approach.Accuracy, approach.Entropy) = MeasureFit(
                    cancellationToken: cancellationToken,
                    data: data,
                    catWeights: weights,
                    catVec: catVec,
                    yVec: yVec,
                    training: approach);
                Chk.NotNull(approach.Entropy, "{0}.{1}", nameof(approach), nameof(approach.Entropy));
                
                approach.AddGoodStep();
                
                // ------------------------------
                // Inner optimization method
                // ------------------------------
                TrainedClassifier attempt = approach.Clone();
                
                int ctStreakOfWanderingSteps = 0;
                int maxStreakOfWanderingSteps = Math.Min(12, numParams);
                double entropyWhenStartedWandering = approach.Entropy.Value;
                
                int ctStreakOfInsufficientLineSearches = 0;
                int maxStreakOfInsufficientLineSearches = 2;
                
                float[][] del;
                float[][] gradient = Util.NewArrays<float>(numPoly, paramScale.Length);

                float stepSize = options.ParamDiffMax;
                while (stepSize > options.ParamDiffMin)
                {
                    // Wander
                    del = wanderer.NextBasis(stepSize);
                    
                    // This is set to true when we need to shrink the step size on the next round.
                    bool shrinkStep = false;
                    
                    // Start with 'approach' which represents our best fit so far (for this approach)
                    Util.Copy(approach.Classifier.Params, attempt.Classifier.Params);
                    
                    // Wander by adding 'del'
                    F.Util.Add(del, attempt.Classifier.Params);
                    
                    // Norm the parameters. 
                    F.Util.Scale(attempt.Classifier.Params, (float)(1.0/F.Util.NormL2(attempt.Classifier.Params)));
                    
                    // Build the classifier.
                    BuildClassifier(
                        cancellationToken: cancellationToken,
                        data: data,
                        regressor: regressor,
                        catWeights: weights,
                        dualCatWeights: dualWeights,
                        catVec: catVec,
                        yVec: yVec,
                        idxVec: idxVec,
                        training: attempt);
                    
                    // Measure the fit.
                    (attempt.Accuracy, attempt.Entropy) = MeasureFit(
                        cancellationToken: cancellationToken,
                        data: data,
                        catWeights: weights,
                        catVec: catVec,
                        yVec: yVec,
                        training: attempt);
                    Chk.NotNull(attempt.Entropy, "{0}.{1}", nameof(attempt), nameof(attempt.Entropy));

                    // Update the gradient.
                    double dEntropy = attempt.Entropy.Value - approach.Entropy.Value;
                    F.Util.AddScaled(-dEntropy, del, gradient);

                    if (dEntropy < 0.0)
                    {
                        // This is a good step.  Update our value of 'approach'.
                        attempt.AddGoodStep();
                        approach = attempt;
                    }
                    else
                    {
                        // This is a bad step.
                        approach.AddBadStep();
                    }
                    
                    // Do we have enough information to do a gradient line search?
                    if (++ctStreakOfWanderingSteps >= maxStreakOfWanderingSteps)
                    {
                        // See if wandering is reducing entropy enough.
                        dEntropy = approach.Entropy.Value - entropyWhenStartedWandering;
                        Chk.LessOrEqual(dEntropy, 0.0, "The change in entropy should be non-positive after wandering.");
                        shrinkStep &= -dEntropy < options.EntropyTol;
                        
                        // Line search the gradient.
                        dEntropy = LineSearch(
                            cancellationToken,
                            data,
                            regressor,
                            weights,
                            dualWeights,
                            catVec,
                            yVec,
                            idxVec,
                            gradient,
                            minStepSize: stepSize/4f,
                            training: ref approach);
                        
                        Chk.LessOrEqual(dEntropy, 0.0, "The change in entropy should be non-positive after a line search.");
                        
                        // If line searching produced an insufficient entropy change ...
                        if (-dEntropy < options.EntropyTol)
                        {
                            ctStreakOfInsufficientLineSearches++;

                            if (ctStreakOfInsufficientLineSearches >= maxStreakOfInsufficientLineSearches)
                            {
                                // Shrink the step.
                                shrinkStep = true;
                                ctStreakOfInsufficientLineSearches = 0;
                            }
                        }

                        // Reset the gradient.
                        ctStreakOfWanderingSteps = 0;
                        Util.SetValues(0f, gradient);
                        entropyWhenStartedWandering = approach.Entropy.Value;
                    }

                    if (shrinkStep)
                    {
                        stepSize /= options.ParamShrinkFactor;
                        shrinkStep = false;
                    }
                }

                output = TrainedClassifier.SelectLowestEntropyAndMerge(output, approach);
            }

            return output;
        }

        /// <summary>
        /// Builds the classifier using the given information.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="data">The training data.</param>
        /// <param name="regressor">The object that performs monotonic regression.</param>
        /// <param name="catWeights">The category weights.</param>
        /// <param name="dualCatWeights">The category weights for the "dual" problems (if they exist).</param>
        /// <param name="catVec">The category label for each training data.  See <see cref="CategorizedData.GetCategoryVector"/>.</param>
        /// <param name="yVec">On input, the working memory must be initialized with the size:  [numPoly][catVec.Length].
        /// On output, this contains the y-values (the polynomial output) for each datum.</param>
        /// <param name="idxVec">On input, the working memory must be initialized with the size:  [numPoly][catVec.Length],
        /// and each polynomial must be initialized with I.Util.FillSeries, or any randomly shuffled ordering of these
        /// values.  On output, this contains the sort order (ascending) for the y-values for each polynomial.</param>
        /// <param name="training">On input, the value of <see cref="TrainedClassifier.Id"/> must be set.  On output,
        /// the classifier is built out and the fit (accuracy, entropy) have been measured.</param>
        private static void BuildClassifier(
            CancellationToken cancellationToken,
            [NotNull] CategorizedData data,
            [NotNull] D1.MonotonicRegressor regressor,
            [NotNull] CategoryWeights catWeights,
            [MaybeNull] CategoryWeights[] dualCatWeights,
            [NotNull] byte[] catVec,
            [NotNull] float[][] yVec,
            [NotNull] int[][] idxVec,
            [NotNull] TrainedClassifier training)
        {
            int? targetCat = training.Id?.TargetCat;

            Chk.False(training.Classifier.NumPoly > 1 && targetCat != null,
                "For multiple polynomials, the target category should be null.");
            
            // First we need to calculate the output of each polynomial expression (y-value) for each datum.
            for (int iPoly = 0; iPoly < training.Classifier.NumPoly; iPoly++)
            {
                int iDatum = 0;

                for (int iCat = 0; iCat < data.NumCats; iCat++)
                {
                    int numInCat = data.NumEach[iCat];
                    for (int jDatum = 0; jDatum < numInCat; jDatum++)
                    {
                        // Evaluate polynomial to get y-values.
                        yVec[iPoly][iDatum++] =
                            (float)training.Classifier.EvalPolyFromExpanded(iPoly, data.X[iCat][jDatum]);
                        
                        // Check for cancellation.
                        cancellationToken.ThrowIfCancellationRequested();
                    }
                }
                
                // Sort by y-value 
                F1.Util.QuickSortIndex(idxVec[iPoly], yVec[iPoly], left: 0, right: data.NumTotal - 1);
                
                // Check for cancellation.
                cancellationToken.ThrowIfCancellationRequested();

                // Quantize and perform monotonic regression.
                if (targetCat != null)
                {
                    // The 2-category "dual" problem.
                    training.Classifier.Quant[iPoly].Measure(
                        cancellationToken: cancellationToken,
                        yIdx: idxVec[iPoly],
                        yValues: yVec[iPoly],
                        cat: catVec,
                        targetCat: (byte)targetCat.Value,
                        catWeights: dualCatWeights[targetCat.Value],
                        regressor: regressor);
                }
                else if (training.Classifier.NumPoly == 1)
                {
                    // A 2-category problem.
                    
                    // The 2-category "dual" problem.
                    training.Classifier.Quant[iPoly].Measure(
                        cancellationToken: cancellationToken,
                        yIdx: idxVec[iPoly],
                        yValues: yVec[iPoly],
                        cat: catVec,
                        targetCat: (byte)iPoly,
                        catWeights: catWeights,
                        regressor: regressor);
                }
                else
                {
                    // A M-category problem, where M > 2
                    training.Classifier.Quant[iPoly].Measure(
                        cancellationToken: cancellationToken,
                        yIdx: idxVec[iPoly],
                        yValues: yVec[iPoly],
                        cat: catVec,
                        targetCat: (byte)iPoly,
                        catWeights: dualCatWeights[iPoly],
                        regressor: regressor);
                }
            }
        }
        
        /// <summary>
        /// Performs a line search along the gradient in an attempt to improve classifier performance.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="data">The training data.</param>
        /// <param name="regressor">The object that performs monotonic regression.</param>
        /// <param name="catWeights">The category weights.</param>
        /// <param name="dualCatWeights">The category weights for the "dual" problems (if they exist).</param>
        /// <param name="catVec">The category label for each training data.  See <see cref="CategorizedData.GetCategoryVector"/>.</param>
        /// <param name="yVec">On input, the working memory must be initialized with the size:  [numPoly][catVec.Length].
        /// On output, this contains the y-values (the polynomial output) for each datum.</param>
        /// <param name="idxVec">On input, the working memory must be initialized with the size:  [numPoly][catVec.Length],
        /// and each polynomial must be initialized with I.Util.FillSeries, or any randomly shuffled ordering of these
        /// values.  On output, this contains the sort order (ascending) for the y-values for each polynomial.</param>
        /// <param name="gradient">The gradient which should decrease entropy (we think).</param>
        /// <param name="minStepSize">The minimum step size.  See <see cref="TrainingOptions.ParamDiffMin"/> to understand the units
        /// that this is expressed in.</param>
        /// <param name="training">The classifier being trained.  If a better fit is found, this will be replaced by a new instance.</param>
        /// <returns>The change in entropy.  This should always be zero (for no change) or negative (for an improvement).</returns>
        private static double LineSearch(
            CancellationToken cancellationToken,
            [NotNull] CategorizedData data,
            [NotNull] D1.MonotonicRegressor regressor,
            [NotNull] CategoryWeights weights,
            [MaybeNull] CategoryWeights[] dualWeights,
            [NotNull] byte[] catVec,
            [NotNull] float[][] yVec,
            [NotNull] int[][] idxVec,
            [NotNull] float[][] gradient,
            [NotNull] float minStepSize,
            [NotNull] ref TrainedClassifier training)
        {
            float step = 1.0f;
            double entropyInput = training.Entropy.Value;
            
            TrainedClassifier attempt = training.Clone();

            while (true)
            {
                // Check for cancellation.
                cancellationToken.ThrowIfCancellationRequested();
                
                // Modify classifier parameters
                F.Util.AddScaled(step, gradient, attempt.Classifier.Params);
                F.Util.Scale(attempt.Classifier.Params, (float)(1.0 / F.Util.NormL2(attempt.Classifier.Params)));
                
                // Build the classifier.
                BuildClassifier(
                    cancellationToken: cancellationToken,
                    data: data,
                    regressor: regressor,
                    catWeights: weights,
                    dualCatWeights: dualWeights,
                    catVec: catVec,
                    yVec: yVec,
                    idxVec: idxVec,
                    training: attempt);
                    
                // Measure the fit.
                (attempt.Accuracy, attempt.Entropy) = MeasureFit(
                    cancellationToken: cancellationToken,
                    data: data,
                    catWeights: weights,
                    catVec: catVec,
                    yVec: yVec,
                    training: attempt);
                Chk.NotNull(attempt.Entropy, "{0}.{1}", nameof(attempt), nameof(attempt.Entropy));

                if (attempt.Entropy.Value < training.Entropy.Value)
                {
                    training = attempt;
                    training.AddGoodStep();
                    step *= 1.4142f;  // sqrt(2), Next try a slightly larger step in the same direction
                }
                else
                {
                    training.AddBadStep();
                    step *= -0.5f;   // Try a smaller step in the opposite direction
                }

                if (Math.Abs(step) < minStepSize)
                {
                    // We break out when the step size gets too small.
                    break;
                }
            }

            // The change in entropy should be non-positive.
            double output = training.Entropy.Value - entropyInput;

            return output;
        }
        
        /// <summary>
        /// Measures the accuracy and entropy of the classifier in training.
        /// </summary>
        /// <param name="cancellationToken">If triggered, this thread will throw a <see cref="OperationCanceledException"/>.</param>
        /// <param name="data">The training data.</param>
        /// <param name="catWeights">The category weights.</param>
        /// <param name="catVec">The category label for each training data.  See <see cref="CategorizedData.GetCategoryVector"/>.</param>
        /// <param name="yVec">On input, the working memory must be initialized with the size:  [numPoly][catVec.Length].
        /// On output, this contains the y-values (the polynomial output) for each datum.</param>
        /// <param name="training">The classifier being trained.  If a better fit is found, this will be replaced by a new instance.</param>
        /// <returns>Accuracy and entropy.</returns>
        private static (double accuracy, double entropy) MeasureFit(
            CancellationToken cancellationToken,
            [NotNull] CategorizedData data,
            [NotNull] CategoryWeights catWeights,
            [NotNull] byte[] catVec,
            [NotNull] float[][] yVec,
            [NotNull] TrainedClassifier training)
        {
            (double accuracy, double entropy) output = (accuracy: 0.0, entropy: 0.0);
            
            float[] y = new float[training.Classifier.NumPoly];
            double[] p;
            
            int? targetCat = training.Id?.TargetCat;
            byte targetCatByte = (byte)(targetCat ?? 0);
            
            for (int iDatum = 0; iDatum < data.NumTotal; iDatum++)
            {
                // Check for cancellation.
                cancellationToken.ThrowIfCancellationRequested();
                
                for (int iPoly = 0; iPoly < y.Length; iPoly++)
                    y[iPoly] = yVec[iPoly][iDatum];
                    
                p = training.Classifier.ClassifyPolynomialOutputs(y);
                byte cat = catVec[iDatum];

                double weight = catWeights.Weights[cat];

                if (targetCat.HasValue)
                {
                    // 2-category dual problem.
                    if (cat == targetCatByte)
                    {
                        output.entropy += Math.Log(p[0]) * weight;

                        if (p[0] > p[1])
                            output.accuracy += weight;
                    }
                    else
                    {
                        output.entropy += Math.Log(p[1]) * weight;
                        
                        if (p[1] > p[2])
                            output.accuracy += weight;
                    }
                }
                else
                {
                    // Normal classification problem (not dual).
                    output.entropy += Math.Log(p[cat]) * weight;
                    int catSelected = D.Util.ArgMax(p);
                    if (catSelected == cat)
                    {
                        output.accuracy += weight;
                    }
                }
            }
            
            // Divide by total weight to get the proportion of weight that was correctly classified.
            output.accuracy /= catWeights.TotalWeight;

            // The value of log(numCats) is used to normalize entropy so that a value of 1 represents "chance"
            // performance.
            output.entropy /= -Math.Log(data.NumCats) * catWeights.TotalWeight;

            return output;
        }
    }
}
