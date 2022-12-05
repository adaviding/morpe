using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Morpe.Validation;

namespace Morpe.Numerics.F2
{
    /// <summary>
    /// A rectangle defined on a normalized screen coordinate system (where Height = Bottom - Top).
    /// </summary>
    public struct Rect
    {
        /// <summary>
        /// The default empty rect which consists entire of NaN values.
        /// </summary>
        public static Rect Empty { get { return new Rect(float.NaN, float.NaN, float.NaN, float.NaN); } }
        
        public float Bottom;
        public float Left;
        public float Right;
        public float Top;
        
        /// <summary>
        /// Right - Left
        /// </summary>
        public float Width { get { return this.Right - this.Left; } }
        
        /// <summary>
        /// Bottom - Top
        /// </summary>
        public float Height { get { return this.Bottom - this.Top; } }
        
        /// <summary>
        /// Returns true if any NaN values are present, false otherwise.
        /// </summary>
        public bool IsEmpty { get { return float.IsNaN(this.Bottom) || float.IsNaN(this.Top) || float.IsNaN(this.Left) || float.IsNaN(this.Right); } }
        
        /// <summary>
        /// Constructs a new Rect.
        /// </summary>
        /// <param name="left"></param>
        /// <param name="top"></param>
        /// <param name="right"></param>
        /// <param name="bottom"></param>
        public Rect(float left, float top, float right, float bottom)
        {
            this.Left = left;
            this.Top = top;
            this.Right = right;
            this.Bottom = bottom;
        }
        
        /// <summary>
        /// Returns true if the values are equal, false otherwise.
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals([NotNull] Rect other)
        {
            Chk.NotNull(other, nameof(other));
            
            return this.Bottom == other.Bottom
                && this.Left == other.Left
                && this.Right == other.Right
                && this.Top == other.Top;
        }
    }
}
