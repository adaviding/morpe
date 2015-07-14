using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe.Draw
{
    /// <summary>
    /// The vertex of a contour.
    /// </summary>
    public class ContourVertex : IComparable
    {
        /// <summary>
        /// The parent of this object.
        /// </summary>
        public IsoContours Parent;
        /// <summary>
        /// Horizontal zero-based floating point number that floats between the zero-based integer indexes of Z.
        /// </summary>
        public double x;
        /// <summary>
        /// Vertical zero-based floating point number that floats between the zero-based integer indexes of Z.
        /// </summary>
        public double y;
        /// <summary>
        /// The direction towards which z increases (in radians).
        /// </summary>
        public double GradientRadians;
        /// <summary>
        /// The gradient magnitude (this is positive only because the gradient direction is pointed towards the positive gradient).
        /// </summary>
        public double dZ_dGradient;
        /// <summary>
        /// The local change in Z as x increases.
        /// </summary>
        public double dZ_dx;
        /// <summary>
        /// The local change in Z as y increases.
        /// </summary>
        public double dZ_dy;
        /// <summary>
        /// The ID number of the contour to which this vertex belongs.
        /// </summary>
        public int IdContour = int.MinValue;
        /// <summary>
        /// The previous contour vertex.  As the linked list is traversed from Prev to Next,
        /// Z increases in the clockwise direciton and decreases in the counter-clockwise direction.
        /// </summary>
        public ContourVertex Prev;
        /// <summary>
        /// The next contour vertex.
        /// </summary>
        public ContourVertex Next;
        /// <summary>
        /// The higher of two affinities.
        /// </summary>
        public float HighestAffinity = float.MinValue;
        /// <summary>
        /// The lower of two affinities.
        /// </summary>
        public float LowestAffinity = float.MinValue;
        public ContourVertex(IsoContours Parent)
        {
            this.Parent = Parent;
        }
        /// <summary>
        /// A string representation of this thing.
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return "[" +
                "(" + ((float)x).ToString() + "," + ((float)y).ToString() + ")"
                + " , " +
                "(" + ((float)dZ_dx).ToString() + "," + ((float)dZ_dy).ToString() + ")"
                + "]";
        }
        private float CompareTo_LowestAffinity;
        /// <summary>
        /// Sort by the lowest affinity.
        /// </summary>
        /// <param name="o"></param>
        /// <returns></returns>
        int IComparable.CompareTo(object o)
        {
            CompareTo_LowestAffinity = ((ContourVertex)o).LowestAffinity;
            if (this.LowestAffinity > CompareTo_LowestAffinity)
                return 1;
            if (this.LowestAffinity == CompareTo_LowestAffinity)
                return 0;
            return -1;
        }
        /*
        public PointF[] Paintable(PaintingInfo pi)
        {
            PointF[] output = new PointF[2];
            float dwi = pi.Rect.Width / Parent.MeshRect.Width;
            float dhi = pi.Rect.Height / Parent.MeshRect.Height;
            output[0].X = ((float)this.x - Parent.MeshRect.X) * dwi + pi.Rect.X;
            output[0].Y = pi.Rect.Height - ((float)this.y - Parent.MeshRect.Y) * dhi + pi.Rect.Y; // Y is inverted
            output[1].X = ((float)this.x - Parent.MeshRect.X + (float)Math.Cos(GradientRadians) * 5.0f) * dwi + pi.Rect.X;
            output[1].Y = pi.Rect.Height - ((float)this.y - Parent.MeshRect.Y + (float)Math.Sin(GradientRadians) * 5.0f) * dhi + pi.Rect.Y; // Y is inverted
            return output;
        }
         */
    }
}
