using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe.Draw
{
    public struct Point
    {
        public static Point Empty { get { return new Point(float.NaN, float.NaN); } }
        public float X;
        public float Y;
        public Point(float X, float Y)
        {
            this.X = X;
            this.Y = Y;
        }
        public bool IsEmpty { get { return float.IsNaN(this.X) || float.IsNaN(this.Y); } }
    }
}
