using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe
{
    /// <summary>
    /// Options for a 2-category MoRPE classifier.
    /// </summary>
    public class ClassifierOptions
    {
        /// <summary>
        /// The rank of the polynomial expansion.
        /// </summary>
        public int Rank;

        /// <summary>
        /// The number of quantiles used for mapping the polynomial output to a conditional probability.
        /// </summary>
        public int NumQuantiles;

        public ClassifierOptions(int rank, int numQuantiles)
        {
            this.Rank = rank;
            this.NumQuantiles = numQuantiles;
        }
    }
}
