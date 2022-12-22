using System.Collections.Generic;
using System.Threading;
using Morpe;
using NUnit.Framework;

using I = Morpe.Numerics.I;

namespace Test.Morpe
{
    public class BasicTests
    {
        [Test]
        public void Small2cat2dim()
        {
            // Generate a very small and simple data set.
            int[] numEach = new[] { 100, 100 };
            int numDims = 2;
            int rank = 2;
            MockGaussianDataSet dataSet = new MockGaussianDataSet(numCats: numEach.Length, numDims: numDims);
            CategorizedData data = dataSet.CreateRandomSample(numEach: numEach);
            
            // Condition and expand the data.
            data.State = new CategorizedDataState(
                conditioner: SpatialConditionMeasurer.Measure(data).Conditioner(),
                poly: new Poly(numDims, rank));
            data.State.Condition(data);
            data.State.Expand(data);
            
            // List the dimensions that will be included in the classifier.
            int[] dims = new int[numDims];
            I.Util.FillSeries(dims);
            
            // Fixme:  PreOptimization analysis.
            PreTrainingAnalysis analysis = PreTrainingAnalysis.BuildFromExpandedData(
                cancellationToken: default,
                data: data);

            // Invoke training
            TrainedClassifier trained = Trainer.Train(
                cancellationToken: default,
                data: data,
                id: new ClassifierId(
                    dims: dims,
                    rank: 2,
                    targetCat: null),
                options: new TrainingOptions(),
                analysis: analysis,
                parameterStart: analysis.ParamInit[rank - 1]);
        }
    }
}
