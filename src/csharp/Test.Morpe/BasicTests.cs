using Morpe;
using NUnit.Framework;

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
            MockGaussianDataSet dataSet = new MockGaussianDataSet(numCats: numEach.Length, numDims: numDims);
            CategorizedData data = dataSet.CreateRandomSample(numEach: numEach);

            ClassifierOptions classifierOps = new ClassifierOptions(
                rank: 1,
                numQuantiles: 14);

            TrainingOptions trainingOps = new TrainingOptions();

            Classifier classifier = new Classifier(
                numCats: numEach.Length,
                numDims: numDims,
                rank: classifierOps.Rank,
                numQuantiles: classifierOps.NumQuantiles);
        }
    }
}
