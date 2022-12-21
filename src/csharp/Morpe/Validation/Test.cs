using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;

namespace Morpe.Validation
{
    public class Test
    {
        public static bool EqualListing<T>(
            [MaybeNull] IReadOnlyList<T> a,
            [MaybeNull] IReadOnlyList<T> b) where T : IEquatable<T>
        {
            if (a == null && b == null)
            {
                return true;
            }
            if (a == null || b == null || a.Count != b.Count)
            {
                return false;
            }

            for (int i = 0; i < a.Count; i++)
            {
                T aa = a[i];
                T bb = b[i];

                if (aa == null && bb == null)
                {
                    continue;
                }
                
                if (aa == null || bb == null || !a.Equals(b))
                {
                    return false;
                }
            }

            return true;
        }
    }
}
