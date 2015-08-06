using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe
{
	public enum OptimizationMode
	{
		/// <summary>
		/// The optimization algorithm is searching along orthogonal bases.
		/// </summary>
		Ortho,
		/// <summary>
		/// The optimization algorithm is moving along a gradient.
		/// </summary>
		Gradient
	}
}
