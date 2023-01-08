using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading;
using Morpe;
using NUnit.Framework;

using I = Morpe.Numerics.I;

namespace Test.Morpe
{
    public class BasicTests
    {

        /// <summary>
        /// This runs a basic test using the given parameters.  For the 2-category problem, the value of d' ("d prime") is sqrt(numDims).
        /// </summary>
        /// <param name="samplesPerCat">The number of data samples per category.</param>
        /// <param name="numCats">The number of categories.</param>
        /// <param name="numDims"></param>
        /// <param name="rank"></param>
        /// <param name="numQuantiles"></param>
        [TestCase(100,   2, 2, 1, 15)]  // The simplest and fastest
        [TestCase(10000, 2, 2, 2, 25)]
        [TestCase(2000,  3, 4, 3, 21)]
        public static void SphericalScenarios(
            int samplesPerCat,
            int numCats,
            int numDims,
            int rank,
            int numQuantiles)
        {
            int[] numEach = new int[numCats];
            Array.Fill(numEach, samplesPerCat);

            RunTest(
                cancellationToken: default,
                mockData: MockGaussianDataSet.CreateSpherical(numEach, numDims, dPrime: 1.0),
                classifierOptions: new ClassifierOptions(
                    rank: rank,
                    numQuantiles: numQuantiles),
                assertFitness: samplesPerCat >= 1000);
        }

        private static void RunTest(
            CancellationToken cancellationToken,
            MockGaussianDataSet mockData,
            ClassifierOptions classifierOptions,
            bool assertFitness)
        {
            CategorizedData data = mockData.Data;
            TrainingOptions trainingOptions = new TrainingOptions();

            // Measure the optimal fit for the given data set.
            CategoryWeights weights = CategoryWeights.Measure(
                numEach: data.NumEach,
                rule: trainingOptions.CategoryWeightingRule,
                targetCategory: null);
            (double accuracy, double entropy) optimalFit = mockData.MeasureOptimalFit(weights);

            Stopwatch stopwatch = Stopwatch.StartNew();

            // Condition and expand the data.
            data.State = new CategorizedDataState(
                conditioner: SpatialConditionMeasurer.Measure(data).Conditioner(),
                polynomial: new Polynomial(mockData.NumDims, classifierOptions.Rank));
            data.State.Condition(data);
            data.State.Expand(data);

            // List the dimensions that will be included in the classifier.
            int[] dims = new int[mockData.NumDims];
            I.Util.FillSeries(dims);

            // PreOptimization analysis.
            PreTrainingAnalysis analysis = PreTrainingAnalysis.BuildFromExpandedData(
                cancellationToken: cancellationToken,
                data: data);

            // Invoke training
            TrainedClassifier trained = Trainer.Train(
                cancellationToken: default,
                data: data,
                id: new ClassifierId(
                    dims: dims,
                    rank: classifierOptions.Rank,
                    targetCat: null),
                options: trainingOptions,
                analysis: analysis,
                parameterStart: analysis.ParamInit[classifierOptions.Rank - 1]);

            TimeSpan elapsed = stopwatch.Elapsed;

            // ------------------------
            // Print stuff before making assertions
            // -------------------------

            double pGoodSteps = (double)trained.NumGoodStepsTaken / trained.NumStepsTaken;
            TestContext.WriteLine($"{pGoodSteps:f4} = {trained.NumGoodStepsTaken} / {trained.NumStepsTaken} = Number of Steps (Good / Total)");
            TestContext.WriteLine($"{trained.NumApproaches} = Number of approaches");
            TestContext.WriteLine($"{elapsed.TotalSeconds} seconds = {elapsed} = Time to train");

            double tPerStep = elapsed.TotalSeconds / trained.NumStepsTaken;
            TestContext.WriteLine($"{tPerStep:f4} = Seconds to train per step");

            TestContext.WriteLine("");
            double ratioAccuracy = (trained.Accuracy ?? double.NaN) / optimalFit.accuracy;
            TestContext.WriteLine($"{ratioAccuracy:f4} = {trained.Accuracy:f4} / {optimalFit.accuracy:f4} = Training Accuracy (Observed / Optimal)");

            double lambdaTrain = Math.Pow(data.NumCats, trained.Entropy ?? double.NaN);
            double lambdaOpt = Math.Pow(data.NumCats, optimalFit.entropy);
            double ratioLambda = lambdaTrain / lambdaOpt;
            TestContext.WriteLine($"{ratioLambda:f4} = {lambdaTrain:f4} / {lambdaOpt:f4} = Training Lambda (Observed / Optimal)");

            TestContext.WriteLine("");
            double diffEntropy = (trained.Entropy ?? double.NaN) - optimalFit.entropy;
            TestContext.WriteLine($"{diffEntropy:f4} = {trained.Entropy:f4} - {optimalFit.entropy:f4} = Training Entropy (Observed - Optimal), lower values are better");

            // ------------------------
            // Assertions
            // ------------------------

            Assert.AreEqual(trained.NumApproaches, trainingOptions.NumberOfApproaches, "The number of approaches was wrong.");

            if (assertFitness)
            {
                // We assert based on fitness when the data set is large enough (i.e. when sampling noise is low enough).

                // We still need to pad our criteria well enough to ensure that the test is not flaky over thousands of iterations.
                Assert.Greater(ratioAccuracy, 0.98, "Training accuracy was lower than expected.");
                Assert.Greater(ratioLambda, 0.98, "Training lambda was lower than expected.");
                Assert.Less(diffEntropy, 0.02, "The entropy was higher than expected.");
            }
        }
    }
}
