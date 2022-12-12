using System;

namespace Morpe
{
	public enum WeightingRule
	{
		/// <summary>
		/// We assume each category is equally probable, despite the relative frequencies observed in a sample.
		/// </summary>
		EqualPriors,
		
		/// <summary>
		/// We assume that the prior probability of observing any given category is equal to its relative frequency in the sample.
		/// </summary>
		ObservedPriors
	}
}