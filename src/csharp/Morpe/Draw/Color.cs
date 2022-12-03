using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe.Draw
{
    public struct Color
    {
        public static Color Black { get { return new Color(0.0f, 0.0f, 0.0f, 1.0f); } }
        public static Color Blue { get { return new Color(0.0f, 0.0f, 0.0f, 1.0f); } }
        public static Color Green { get { return new Color(0.0f, 0.0f, 0.8f, 0.0f); } }
        public static Color Red { get { return new Color(0.0f, 1.0f, 0.0f, 0.0f); } }
        /// <summary>
        /// Alpha channel value in the range [0, 1].  Min = 0 = Opaque, Max = 1 = Transparent.
        /// </summary>
        public float A;
        /// <summary>
        /// Red channel value in the range [0, 1].
        /// </summary>
        public float R;
        /// <summary>
        /// Green channel value in the range [0, 1].
        /// </summary>
        public float G;
        /// <summary>
        /// Blue channel value in the range [0, 1].
        /// </summary>
        public float B;
        /// <summary>
        /// Constructs a color having the given arguments.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="R"></param>
        /// <param name="G"></param>
        /// <param name="B"></param>
        public Color(float A, float R, float G, float B)
        {
            this.A = A;
            this.R = R;
            this.G = G;
            this.B = B;
        }
        /// <summary>
        /// Returns a gamma-adjusted color.
        /// </summary>
        /// <param name="gamma">The gamma value.</param>
        /// <returns>The gamma-adjusted color.</returns>
        public Color Pow(double gamma)
        {
            return new Color(this.A, (float)Math.Pow(this.R, gamma), (float)Math.Pow(this.G, gamma), (float)Math.Pow(this.B, gamma));
        }
    }
}
