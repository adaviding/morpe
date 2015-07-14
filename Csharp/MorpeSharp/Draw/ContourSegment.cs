using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe.Draw
{
    /// <summary>
    /// Represents the affinity with which two points (A and B) can join to form a contour segment.
    /// </summary>
    public class ContourSegment : IComparable
    {
        /// <summary>
        /// The parent of this object.
        /// </summary>
        public IsoContours Parent;
        /// <summary>
        /// The starting vertex.
        /// </summary>
        public ContourVertex A;
        /// <summary>
        /// The ending vertex.
        /// </summary>
        public ContourVertex B;
        /// <summary>
        /// A positive number in [0, inf).  The higher this number, the more likely that two adjacent
        /// vertices should be grouped.
        /// The affinity is a function of the agreement in angles and the distance between vertices.
        /// It's just a formula that I made up ad-hoc.  It should perform pretty well, but there may
        /// be room for improvement.
        /// </summary>
        public double Affinity = Double.NegativeInfinity;
        /// <summary>
        /// The change in x moving from A to B.
        /// </summary>
        public double dx;
        /// <summary>
        /// The change in y moving from A to B.
        /// </summary>
        public double dy;
        /// <summary>
        /// The distance between vertices A and B.
        /// </summary>
        public double ds;
        /// <summary>
        /// The change in x per unit distance moving from A to B.
        /// </summary>
        public double dx_ds;
        /// <summary>
        /// The change in y per unit distance moving from A to B.
        /// </summary>
        public double dy_ds;
        /// <summary>
        /// The change in Z per unit of the vector A-->B rotated by +90 degrees;  This is calculated
        /// by rotating this.dZ_dx and this.dZ_dy relative to A-->B such that A-->B is at 0 degrees.
        /// </summary>
        public double dZ_dCosAB, dZ_dSinAB;
        /// <summary>
        /// This is calculated by averaging dZ_dx and dZ_dy from A and B.
        /// </summary>
        public double dZ_dx, dZ_dy;
        /// <summary>
        /// The contour segment joining vertices A and B moving from A to B.  The constructed segment will
        /// map the arguemnt A to either this.A or this.B and map the argument B to the other member
        /// this.A or this.B such that dZ_dSinAB less or equal 0.0.  Thus, the segment A-->B is always 
        /// compatible with a clockwise polygon (where Z is larger inside the polygon) and where A.Next is B
        /// and B.Prev is A.
        /// </summary>
        /// <param name="A">The starting vertex.</param>
        /// <param name="B">The ending vertex.</param>
        public ContourSegment(IsoContours Parent, ContourVertex A, ContourVertex B)
        {
            this.Parent = Parent;
            dx = B.x - A.x;
            dy = B.y - A.y;
            ds = Math.Sqrt(dx * dx + dy * dy);
            dx_ds = dx / ds;
            dy_ds = dy / ds;
            dZ_dx = 0.5 * (A.dZ_dx + B.dZ_dx);
            dZ_dy = 0.5 * (A.dZ_dy + B.dZ_dy);
            dZ_dCosAB = dZ_dx * dx_ds + dZ_dy * dy_ds;
            dZ_dSinAB = -dZ_dx * dy_ds + dZ_dy * dx_ds;
            if (dZ_dSinAB <= 0.0)
            {
                this.A = A;
                this.B = B;
            }
            else
            {
                this.A = B;
                this.B = A;
                dx = -dx;
                dy = -dy;
                dx_ds = -dx_ds;
                dy_ds = -dy_ds;
                dZ_dx = -dZ_dx;
                dZ_dy = -dZ_dy;
                dZ_dCosAB = -dZ_dCosAB;
                dZ_dSinAB = -dZ_dSinAB;
            }
            this.Affinity = (1.0 + (A.dZ_dx * B.dZ_dx + A.dZ_dy * B.dZ_dy) / (A.dZ_dGradient * B.dZ_dGradient)) / ds;
            //this.Affinity = 1.0/ds;
        }
        /// <summary>
        /// Join vertices by linking them in a double-linked-list and recording the linking affinity.
        /// </summary>
        public void JoinVertices()
        {
            this.A.Next = this.B;
            if (this.Affinity > (double)this.A.HighestAffinity)
            {
                this.A.LowestAffinity = this.A.HighestAffinity;
                this.A.HighestAffinity = (float)this.Affinity;
            }
            else
            {
                if (this.Affinity > (double)this.A.LowestAffinity)
                    this.A.LowestAffinity = (float)this.Affinity;
            }

            this.B.Prev = this.A;
            if (this.Affinity > (double)this.B.HighestAffinity)
            {
                this.B.LowestAffinity = this.B.HighestAffinity;
                this.B.HighestAffinity = (float)this.Affinity;
            }
            else
            {
                if (this.Affinity > (double)this.B.LowestAffinity)
                    this.B.LowestAffinity = (float)this.Affinity;
            }
        }
        /// <summary>
        /// Join vertices if they are currently unjoined.
        /// </summary>
        public void JoinUnjoinedVertices()
        {
            if (this.A.Next == null && this.B.Prev == null)
                this.JoinVertices();
        }
        private double comparisonAffinity;
        int IComparable.CompareTo(object obj)
        {
            comparisonAffinity = ((ContourSegment)obj).Affinity;
            if (Affinity > comparisonAffinity)
                return 1;
            if (Affinity == comparisonAffinity)
                return 0;
            return -1;
        }
        /*
        public PointF[] Paintable(PaintingInfo pi)
        {
            PointF[] output = new PointF[2];
            float dwi = pi.Rect.Width / Parent.MeshRect.Width;
            float dhi = pi.Rect.Height / Parent.MeshRect.Height;
            output[0].X = ((float)A.x - Parent.MeshRect.X) * dwi + pi.Rect.X;
            output[0].Y = pi.Rect.Height - ((float)A.y - Parent.MeshRect.Y) * dhi + pi.Rect.Y; // Y is inverted
            output[1].X = ((float)A.x - Parent.MeshRect.X) * dwi + pi.Rect.X;
            output[1].Y = pi.Rect.Height - ((float)A.y - Parent.MeshRect.Y) * dhi + pi.Rect.Y; // Y is inverted
            return output;
        }
         */
    }
}
