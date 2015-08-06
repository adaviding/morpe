using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe.Distributions.D1
{
	/// <summary>
	/// A class for computing things related to the Gaussian (standard normal) probability function.
	/// </summary>
	public class Gaussian
	{
		public static readonly double InvPdfScalar;
		public static readonly double PdfScalar;
		public static readonly double Root2;
		public static readonly double InvRoot2;
		public static readonly double InvRootPi;
		public static readonly double MidCut;
		static Gaussian()
		{
			InvPdfScalar = Math.Sqrt(2.0 * Math.PI);
			PdfScalar = 1.0 / InvPdfScalar;
			Root2 = Math.Sqrt(2);
			InvRoot2 = 1.0 / Root2;
			InvRootPi = 1.0 / Math.Sqrt(Math.PI);
			MidCut = 0.46875 * Root2;
		}
		private static double y, x, qq, u, t;
		#region Normal CDF constants
		static readonly double[] a = new double[] {
			1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
			1.887426188426510e+002,3.209377589138469e+003
		};
		static readonly double[] b = new double[] {
			1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
			1.813893686502485e+003,8.044716608901563e+003
		};
		static readonly double[] c = new double[] {
			2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
			6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
			1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
		};
		static readonly double[] d = new double[] {
			1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
			5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
			4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
		};
		static readonly double[] p = new double[] {
			1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
			1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
		};
		static readonly double[] q = new double[] {
			1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
			5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
		};
		#endregion
		#region Inverse CDF constants
		static readonly double[] ia = new double[] {
			-3.969683028665376e+01,  2.209460984245205e+02,
			-2.759285104469687e+02,  1.383577518672690e+02,
			-3.066479806614716e+01,  2.506628277459239e+00
		};
		static readonly double[] ib = new double[] {
			-5.447609879822406e+01,  1.615858368580409e+02,
			-1.556989798598866e+02,  6.680131188771972e+01,
			-1.328068155288572e+01
		};
		static readonly double[] ic = new double[] {
			-7.784894002430293e-03, -3.223964580411365e-01,
			-2.400758277161838e+00, -2.549732539343734e+00,
			4.374664141464968e+00,  2.938163982698783e+00
		};
		static readonly double[] id = new double[] {
			7.784695709041462e-03,  3.224671290700398e-01,
			2.445134137142996e+00,  3.754408661907416e+00
		};
		#endregion
		/// <summary>
		/// The probability density function of z.
		/// </summary>
		/// <param name="z">A variable in the range [-inf, +inf]</param>
		/// <returns>The probability density.</returns>
		public static double Pdf(double z)
		{
			return Math.Exp(-0.5 * z * z) * PdfScalar;
		}
		/// <summary>
		/// The cumulative probability density function of z.
		/// </summary>
		/// <param name="z">A variable in the range [-inf, +inf]</param>
		/// <returns>The cumulative probability density in the range [0, 1].</returns>
		public static double Cdf(double z)
		{
			y = Math.Abs(z);
			if (y <= MidCut)
			{
				/* evaluate erf() for |z| <= sqrt(2)*0.46875 */
				x = y * y;
				y = z * ((((a[0] * x + a[1]) * x + a[2]) * x + a[3]) * x + a[4]) / ((((b[0] * x + b[1]) * x + b[2]) * x + b[3]) * x + b[4]);
				return 0.5 + y;
			}
			x = 0.5 * Math.Exp(-0.5 * y * y);
			if (y <= 4.0)
			{
				/* evaluate erfc() for sqrt(2)*0.46875 <= |z| <= sqrt(2)*4.0 */
				y *= InvRoot2;
				y = ((((((((c[0] * y + c[1]) * y + c[2]) * y + c[3]) * y + c[4]) * y + c[5]) * y + c[6]) * y + c[7]) * y + c[8])
						/
					((((((((d[0] * y + d[1]) * y + d[2]) * y + d[3]) * y + d[4]) * y + d[5]) * y + d[6]) * y + d[7]) * y + d[8]);
				y *= x;
			}
			else
			{
				/* evaluate erfc() for |z| > sqrt(2)*4.0 */
				x *= Root2 / y;
				y = 2.0 / (y * y);
				y *= (((((p[0] * y + p[1]) * y + p[2]) * y + p[3]) * y + p[4]) * y + p[5])
						/
					(((((q[0] * y + q[1]) * y + q[2]) * y + q[3]) * y + q[4]) * y + q[5]);
				y = x * (InvRootPi - y);
			}
			return (z < 0.0 ? y : 1.0 - y);
		}
		/// <summary>
		/// The inverse cumulative probability density of p.
		/// </summary>
		/// <param name="u">The cumulative probability density in the range [0, 1]</param>
		/// <returns>A variable in the range [-inf, +inf]</returns>
		public static double InvCdf(double p)
		{
			if (p == 0.0)
				return Double.NegativeInfinity;
			if (p == 1.0)
				return Double.PositiveInfinity;
			qq = Math.Min(p, 1.0 - p);
			if (qq > 0.02425)
			{
				/* Rational approximation for central region. */
				u = qq - 0.5;
				t = u * u;
				u = u * (((((ia[0] * t + ia[1]) * t + ia[2]) * t + ia[3]) * t + ia[4]) * t + ia[5])
						/
					(((((ib[0] * t + ib[1]) * t + ib[2]) * t + ib[3]) * t + ib[4]) * t + 1);
			}
			else
			{
				/* Rational approximation for tail region. */
				t = Math.Sqrt(-2.0 * Math.Log(qq));
				u = (((((ic[0] * t + ic[1]) * t + ic[2]) * t + ic[3]) * t + ic[4]) * t + ic[5])
						/
					((((id[0] * t + id[1]) * t + id[2]) * t + id[3]) * t + 1);
			}
			/* The relative error of the approximation has absolute value less
			than 1.15e-9.  One iteration of Halley's rational method (third
			order) gives full machine precision... */
			t = Cdf(u) - qq;    /* error */
			t = t * InvPdfScalar * Math.Exp(0.5 * u * u);   /* f(u)/df(u) */
			u = u - t / (1.0 + 0.5 * u * t);     /* Halley's method */

			return (p > 0.5 ? -u : u);
		}
	}
}
