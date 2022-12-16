using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe
{
    /// <summary>
    /// Holds a fitted instance of type T.  This is useful in the context of optimization.
    /// </summary>
    /// <typeparam name="T">The type of instance being fitted.</typeparam>
    public class Fitted<T>
    {
        /// <summary>
        /// The instance that was fitted.
        /// </summary>
        public T Instance;
        /// <summary>
        /// The fit.  This is the conditional entropy of the training sample (as defined in the paper).
        /// A value of NaN indicates that a fit has not been performed yet.
        /// </summary>
        public double Fit;
        /// <summary>
        /// Constructs a new wrapper for the instance.
        /// </summary>
        /// <param name="instance">The instance that was fitted.</param>
        public Fitted(T instance, double fit)
        {
            this.Instance = instance;
            this.Fit = double.NaN;
        }
    }
}
