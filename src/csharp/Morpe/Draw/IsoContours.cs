using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using D = Morpe.Numerics.D;
using F2 = Morpe.Numerics.F2;

namespace Morpe.Draw
{
    public class IsoContours : IComparable<IsoContours>
    {
        /// <summary>
        /// The level of 'Z' for which these contours are defined.
        /// </summary>
        public double IsoZ;
        /// <summary>
        /// The data for which the contour set is defined.  Z[i,j] is a tabulation of the function Z(y,x).
        /// </summary>
        public double[,] Z;
        /// <summary>
        /// The complete list contours that trace an IsoZ value sorted from longest to shortest.
        /// </summary>
        public List<Contour> Contours = new List<Contour>();
        /// <summary>
        /// The vertex pairs computed as the contours were being built.
        /// </summary>
        public List<ContourSegment> Segments = new List<ContourSegment>();
        /// <summary>
        /// The vertices computed as the conoturs were being built.
        /// </summary>
        public List<ContourVertex> Vertices = new List<ContourVertex>();
        /// <summary>
        /// The rectangle that delimits the extent of the histogram mesh.
        /// </summary>
        public F2.Rect MeshRect;
        /// <summary>
        /// Constructs the iso-value contours for the function tabulated in Z.
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="IsoZ">The value of Z isolated by the contours.</param>
        public IsoContours(double[,] Z, double IsoZ)
        {
            this.Z = Z;
            this.IsoZ = IsoZ;
            int iLen = Z.GetLength(0);
            int jLen = Z.GetLength(1);
            this.MeshRect = new F2.Rect(0.0f, 0.0f, Z.GetLength(1) - 1, Z.GetLength(0) - 1);
            ContourVertex[,] zCrossings = new ContourVertex[iLen, jLen];
            double u0, v0, u1, v1, z00, z01, z10, z11, dz0, dz1;
            bool[] crosses = new bool[4];
            ContourVertex cv;
            int i, j, ii, jj, ctCrosses, ia, ja;
            for (ii = 1; ii < iLen; ii++)
            {
                i = ii - 1;
                for (jj = 1; jj < jLen; jj++)
                {
                    //	In the comments below, let
                    //		"Left" refers to lower col indexes of 'Z'
                    //		"Right" refers to higher col indexes of 'Z'
                    //		"Top" refers to higher row indexes of 'Z'
                    //		"Bottom" refers to lower row indexes of 'Z'
                    j = jj - 1;
                    ctCrosses = 0;
                    //zmid = 0.25* (Z[ii,jj] + Z[i,jj] + Z[ii,j] + Z[i,j]);
                    z00 = Z[i, j];
                    z10 = Z[ii, j];
                    z01 = Z[i, jj];
                    z11 = Z[ii, jj];
                    //	Left
                    crosses[0] = (z00 < IsoZ) != (z10 < IsoZ);
                    //	Right
                    crosses[1] = (z01 < IsoZ) != (z11 < IsoZ);
                    //	Top
                    crosses[2] = (z10 < IsoZ) != (z11 < IsoZ);
                    //	Bottom
                    crosses[3] = (z00 < IsoZ) != (z01 < IsoZ);
                    //	Count crossings
                    if (crosses[0]) ctCrosses++;
                    if (crosses[1]) ctCrosses++;
                    if (crosses[2]) ctCrosses++;
                    if (crosses[3]) ctCrosses++;
                    //  ctCrosses==0 --> No intersection (do nothing)
                    //	ctCrosses==1 --> Not possible (do nothing)
                    //	ctCrosses==3 --> Not possible (do nothing)
                    //	ctCrosses==4 --> Ambiguous  (do nothing).
                    if (ctCrosses == 2)	//	The contour intersects the current cell.
                    {
                        cv = new ContourVertex(this);
                        if (crosses[0] && crosses[1]) //	Horizontal signature
                        {
                            //	Left
                            dz0 = z10 - z00;
                            u0 = Math.Abs(z00 - IsoZ);
                            v0 = Math.Abs(z10 - IsoZ);
                            u0 = u0 / (u0 + v0);
                            //	Right
                            dz1 = z11 - z01;
                            u1 = Math.Abs(z01 - IsoZ);
                            v1 = Math.Abs(z11 - IsoZ);
                            u1 = u1 / (u1 + v1);
                            //	Position
                            cv.x = (double)j + 0.5;
                            cv.y = (double)i + 0.5 * (u0 + u1);
                            //	Angle
                            if (dz0 > 0.0)
                                cv.GradientRadians = Math.Atan2(u1 - u0, 1.0) + D.Util.Pi_1_2;
                            else
                                cv.GradientRadians = Math.Atan2(u1 - u0, 1.0) - D.Util.Pi_1_2;
                            //	Gradient
                            cv.dZ_dy = 0.5 * (dz0 + dz1);
                            cv.dZ_dGradient = cv.dZ_dy / Math.Sin(cv.GradientRadians);
                            cv.dZ_dx = cv.dZ_dGradient * Math.Cos(cv.GradientRadians);
                            //	Storage
                            zCrossings[i, j] = cv;
                            Vertices.Add(cv);
                        }
                        else if (crosses[2] && crosses[3])	//	Vertical signature
                        {
                            //	Bottom
                            dz0 = z01 - z00;
                            u0 = Math.Abs(z00 - IsoZ);
                            v0 = Math.Abs(z01 - IsoZ);
                            u0 = u0 / (u0 + v0);
                            //	Top
                            dz1 = z11 - z10;
                            u1 = Math.Abs(z10 - IsoZ);
                            v1 = Math.Abs(z11 - IsoZ);
                            u1 = u1 / (u1 + v1);
                            //	Position
                            cv.x = (double)j + 0.5 * (u0 + u1);
                            cv.y = (double)i + 0.5;
                            //	Angle
                            if (dz0 > 0.0)
                                cv.GradientRadians = Math.Atan2(u0 - u1, 1.0);
                            else
                                cv.GradientRadians = Math.Atan2(u0 - u1, 1.0) + Math.PI;
                            //	Gradient
                            cv.dZ_dx = 0.5 * (dz0 + dz1);
                            cv.dZ_dGradient = cv.dZ_dx / Math.Cos(cv.GradientRadians);
                            cv.dZ_dy = cv.dZ_dGradient * Math.Sin(cv.GradientRadians);
                            //	Storage
                            zCrossings[i, j] = cv;
                            Vertices.Add(cv);
                        }
                        else //	Mixed signature
                        {
                            if (crosses[0]) // Left
                            {
                                if (crosses[2]) // Top
                                {
                                    //	Top-left corner is cut off z10
                                    //	Top
                                    dz0 = z11 - z10; // dz_dx along top
                                    u0 = u1 = Math.Abs(z10 - IsoZ);
                                    v0 = Math.Abs(z11 - IsoZ);
                                    u0 = u0 / (u0 + v0);
                                    //	Left
                                    dz1 = z10 - z00; // dz_dy along left
                                    //u1 = Math.Abs(z10-IsoZ);
                                    v1 = Math.Abs(z00 - IsoZ);
                                    u1 = u1 / (u1 + v1);
                                    //	Position  (midpoint of hypotenuse)
                                    cv.x = (double)j + 0.5 * u0;  // distance from left
                                    cv.y = (double)ii - 0.5 * u1;  // distance from top
                                    //	Angle
                                    if (dz1 > 0.0)
                                        cv.GradientRadians = Math.Atan2(u1, u0) + D.Util.Pi_1_2;
                                    else
                                        cv.GradientRadians = Math.Atan2(u1, u0) - D.Util.Pi_1_2;
                                    //	Gradient
                                    cv.dZ_dx = dz0;
                                    cv.dZ_dy = dz1;
                                    if (Math.Abs(dz0) > Math.Abs(dz1))
                                        cv.dZ_dGradient = cv.dZ_dx / Math.Cos(cv.GradientRadians);
                                    else
                                        cv.dZ_dGradient = cv.dZ_dy / Math.Sin(cv.GradientRadians);
                                    //	Storage
                                    zCrossings[i, j] = cv;
                                    Vertices.Add(cv);
                                }
                                else if (crosses[3]) // Bottom
                                {
                                    //	Bottom-left corner is cut off  z00
                                    //	Bottom
                                    dz0 = z01 - z00; // dz_dx along bottom
                                    u0 = u1 = Math.Abs(z00 - IsoZ);
                                    v0 = Math.Abs(z01 - IsoZ);
                                    u0 = u0 / (u0 + v0);
                                    //	Left
                                    dz1 = z10 - z00; // dz_dy along left
                                    //u1 = Math.Abs(z00-IsoZ);
                                    v1 = Math.Abs(z10 - IsoZ);
                                    u1 = u1 / (u1 + v1);
                                    //	Position  (midpoint of hypotenuse)
                                    cv.x = (double)j + 0.5 * u0; // distance from left
                                    cv.y = (double)i + 0.5 * u1; // distance from bottom
                                    //	Angle
                                    if (dz1 > 0.0)
                                        cv.GradientRadians = D.Util.Pi_1_2 - Math.Atan2(u1, u0); //Math.PI-Math.Atan2(u1,u0)-RADIANS_1_4;
                                    else
                                        cv.GradientRadians = D.Util.Pi_3_2 - Math.Atan2(u1, u0); //Math.PI-Math.Atan2(u1,u0)+RADIANS_1_4;
                                    //	Gradient
                                    cv.dZ_dx = dz0;
                                    cv.dZ_dy = dz1;
                                    if (Math.Abs(dz0) > Math.Abs(dz1))
                                        cv.dZ_dGradient = cv.dZ_dx / Math.Cos(cv.GradientRadians);
                                    else
                                        cv.dZ_dGradient = cv.dZ_dy / Math.Sin(cv.GradientRadians);
                                    //	Storage
                                    zCrossings[i, j] = cv;
                                    Vertices.Add(cv);
                                }
                            }
                            if (crosses[1]) // Right
                            {
                                if (crosses[2]) // Top
                                {
                                    //	Top-right corner is cut off z11
                                    //	Top
                                    dz0 = z11 - z10; // dz_dx along top
                                    u0 = u1 = Math.Abs(z11 - IsoZ);
                                    v0 = Math.Abs(z10 - IsoZ);
                                    u0 = u0 / (u0 + v0);
                                    //	Right
                                    dz1 = z11 - z01; // dz_dy along right
                                    //u1 = Math.Abs(z11-IsoZ);
                                    v1 = Math.Abs(z01 - IsoZ);
                                    u1 = u1 / (u1 + v1);
                                    //	Position  (midpoint of hypotenuse)
                                    cv.x = (double)jj - 0.5 * u0; // distance from right
                                    cv.y = (double)ii - 0.5 * u1; // distance from top
                                    //	Angle
                                    if (dz1 > 0.0)
                                        cv.GradientRadians = D.Util.Pi_3_2 - Math.Atan2(u1, u0); //Math.PI-Math.Atan2(u1,u0)-RADIANS_1_4;
                                    else
                                        cv.GradientRadians = D.Util.Pi_3_2 - Math.Atan2(u1, u0); //Math.PI-Math.Atan2(u1,u0)+RADIANS_1_4;
                                    //	Gradient
                                    cv.dZ_dx = dz0;
                                    cv.dZ_dy = dz1;
                                    if (Math.Abs(dz0) > Math.Abs(dz1))
                                        cv.dZ_dGradient = cv.dZ_dx / Math.Cos(cv.GradientRadians);
                                    else
                                        cv.dZ_dGradient = cv.dZ_dy / Math.Sin(cv.GradientRadians);
                                    //	Storage
                                    zCrossings[i, j] = cv;
                                    Vertices.Add(cv);
                                }
                                else if (crosses[3]) // Bottom
                                {
                                    //	Bottom-right corner is cut off  z01
                                    //	Bottom
                                    dz0 = z01 - z00; // dz_dx
                                    u0 = u1 = Math.Abs(z01 - IsoZ);
                                    v0 = Math.Abs(z00 - IsoZ);
                                    u0 = u0 / (u0 + v0);
                                    //	Right
                                    dz1 = z11 - z01; // dz_dy
                                    //u1 = Math.Abs(z01-IsoZ);
                                    v1 = Math.Abs(z11 - IsoZ);
                                    u1 = u1 / (u1 + v1);
                                    //	Position  (midpoint of hypotenuse)
                                    cv.x = (double)jj - 0.5 * u0;  // distance from right
                                    cv.y = (double)i + 0.5 * u1;  // distance from bottom
                                    //	Angle
                                    if (dz1 > 0.0)
                                        cv.GradientRadians = Math.Atan2(u1, u0) + D.Util.Pi_1_2;
                                    else
                                        cv.GradientRadians = Math.Atan2(u1, u0) - D.Util.Pi_1_2;
                                    //	Gradient
                                    cv.dZ_dx = dz0;
                                    cv.dZ_dy = dz1;
                                    if (Math.Abs(dz0) > Math.Abs(dz1))
                                        cv.dZ_dGradient = cv.dZ_dx / Math.Cos(cv.GradientRadians);
                                    else
                                        cv.dZ_dGradient = cv.dZ_dy / Math.Sin(cv.GradientRadians);
                                    //	Storage
                                    zCrossings[i, j] = cv;
                                    Vertices.Add(cv);
                                }
                            }
                        }
                        //	ContourSegment relations
                        ia = i - 1;
                        ja = j - 1;
                        if (ia >= 0 && zCrossings[ia, j] != null)
                            Segments.Add(new ContourSegment(this, cv, zCrossings[ia, j]));
                        if (ia >= 0 && ja >= 0 && zCrossings[ia, ja] != null)
                            Segments.Add(new ContourSegment(this, cv, zCrossings[ia, ja]));
                        if (ja >= 0 && zCrossings[i, ja] != null)
                            Segments.Add(new ContourSegment(this, cv, zCrossings[i, ja]));
                    }
                }
            }
            //	Sort the potential contour segments by ascending affinity
            Segments.Sort();
            //	Segment processing logic.
            i = Segments.Count - 1;
            //	Trim any segments that have overlapping vertices
            while (i >= 0 && Segments[i].ds <= 0.0)
                Segments.RemoveAt(i--);
            //	Join the unjoined vertices.
            while (i >= 0)
                Segments[i--].JoinUnjoinedVertices();
            //	Contour formation logic, begin by sorting vertices by descending LowestAffinity.
            Vertices.Sort();
            Vertices.Reverse();
            int nMax = Vertices.Count;
            j = 0;
            for (i = 0; i < Vertices.Count && nMax > 3; i++)
            {
                if (Vertices[i].IdContour == int.MinValue)
                {
                    Contours.Add(new Contour(j, Vertices[i], nMax, MeshRect));
                    nMax -= Contours[j].Vertices.Length;
                    j++;
                }
            }
            //	Sort by contour length descending (longest contour first, shortest last)
            Contours.Sort();
            Contours.Reverse();
            //	Adjust the vertex coordinates to the DataRect
            for (i = 0; i < Contours.Count; i++)
                Contours[i].IdContour = i;
        }
        private double comparisonZ;
        int IComparable<IsoContours>.CompareTo(IsoContours other)
        {
            if (other == null)
                return 1;
            if (IsoZ > other.IsoZ)
                return 1;
            if (IsoZ == other.IsoZ)
                return 0;
            return -1;
        }
        /// <summary>
        /// The pen used for drawing the axis, width is automatically scaled by DrawScale.
        /// </summary>
        public Pen Pen
        {
            get { return pen; }
            set
            {
                pen = value;
                scaledPen = new Pen(pen.Color, pen.Width * drawScale);
            }
        }
        private Pen pen = new Pen(Color.Black, 1.0f);
        private Pen scaledPen = new Pen(Color.Black, 1.0f);
        /// <summary>
        /// The scale with which to draw the axis.  Automatically scales the font and the pen.
        /// </summary>
        public float DrawScale
        {
            get { return drawScale; }
            set
            {
                if (value != drawScale)
                {
                    drawScale = value;
                    scaledPen = new Pen(pen.Color, pen.Width * drawScale);
                }
            }
        }
        private float drawScale = 1.0f;
    }
}
