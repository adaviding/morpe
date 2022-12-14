using System;

namespace Morpe.Numerics.D1
{
    public class ExponentialDistribution
    {
        public static double Cdf(double x)
        {
            return 1.0 - Math.Exp(-x);
        }

        public static double Cdf(double x, double lambda)
        {
            return 1.0 - Math.Exp(-x * lambda);
        }
        
        public static double InvCdf(double u)
        {
            return -Math.Log(1.0-u);
        }
        
        public static double InvCdf(double u, double lambda)
        {
            return -Math.Log(1.0-u)/lambda;
        }
        
        public static double InvPdf(double d)
        {
            return -Math.Log(d);
        }
        
        public static double InvPdf(double d, double lambda)
        {
            return -Math.Log(d/lambda)/lambda;
        }
        
        public static double Pdf(double x)
        {
            return Math.Exp(-x);
        }
        
        public static double Pdf(double x, double lambda)
        {
            return lambda * Math.Exp(-x * lambda);
        }
        
        /// <summary>
        /// Solves for lambda given 'x' and 'y'.
        /// </summary>
        /// <param name="x">The independent variable of the cumulative distribution function.  Non-negative.</param>
        /// <param name="u">The output of the cumulative distribution function.  Must be in the range [0, 1].</param>
        /// <returns>The value of lambda that corresponds to a given value of 'x' and 'u'.</returns>
        public static double SolveForLambda(double x, double u)
        {
            if (u >= 1.0)
                return Double.PositiveInfinity;
            return -Math.Log(1.0-u)/x;
        }
    }
}
