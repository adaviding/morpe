using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using F2 = Morpe.Numerics.F2;

namespace Morpe.Draw
{
    /// <summary>
    /// A contour.
    /// </summary>
    public class Contour : IComparable
    {
        /// <summary>
        /// The ID tag associated with this contour.  When you set this property, the IdContour property
        /// is automatically set for all elements in this.Vertices.
        /// </summary>
        public int IdContour
        {
            get
            {
                return this.idContour;
            }
            set
            {
                this.idContour = value;
                for (int i = 0; i < Vertices.Length; i++)
                    Vertices[i].IdContour = value;
            }
        }
        private int idContour = int.MinValue;
        public bool IsClosed = false;
        public F2.Rect DataRect;
        public ContourVertex[] Vertices;
        /// <summary>
        /// Construct a new contour with the specified ID 
        /// </summary>
        /// <param name="IdContour">An integer ID number associated with this contour.</param>
        /// <param name="cv">A ContourVertex in the contour.</param>
        /// <param name="nMax">The maximum number of vertices that the contour could contain.</param>
        public Contour(int IdContour, ContourVertex cv, int nMax, F2.Rect DataRect)
        {
            this.DataRect = DataRect;
            this.idContour = IdContour;
            if (cv.IdContour >= 0 && cv.IdContour != idContour)
                Vertices = new ContourVertex[0];
            int ct = 1;
            ContourVertex cvFirst = cv;
            ContourVertex cvLast = cv;
            //	Make cvLast keep going until end.
            while (cvLast.Next != null)
            {
                if (cvLast.Next == cvFirst)
                {
                    this.IsClosed = true;
                    break;
                }
                else
                {
                    ct++;
                    cvLast = cvLast.Next;
                }
                if (ct > nMax)
                    throw new Exception("Loop continues too far.");
            }
            if (!IsClosed)
            {
                while (cvFirst.Prev != null)
                {
                    ct++;
                    cvFirst = cvFirst.Prev;
                    if (ct > nMax)
                        throw new Exception("Loop continues too far.");
                }
            }
            Vertices = new ContourVertex[ct];
            cv = cvFirst;
            for (int i = 0; i < ct; i++)
            {
                Vertices[i] = cv;
                cv.IdContour = IdContour;
                cv = cv.Next;
            }
        }
        private F2.Point[] paintCache;
        private F2.Rect lastPaintRect = F2.Rect.Empty;
        /*
        public PointF[] Paintable(PaintingInfo pi)
        {
            if (paintCache == null || paintCache.Length != Vertices.Length)
            {
                lastPaintRect = RectangleF.Empty;
                paintCache = new PointF[Vertices.Length];
            }
            if (!lastPaintRect.Equals(pi.Rect))
            {
                lastPaintRect = pi.Rect;
                PointF pf = PointF.Empty;
                ContourVertex cv;
                float dwi = pi.Rect.Width / DataRect.Width;
                float dhi = pi.Rect.Height / DataRect.Height;
                for (int i = 0; i < Vertices.Length; i++)
                {
                    cv = Vertices[i];
                    pf.X = ((float)cv.x - DataRect.X) * dwi + pi.Rect.X;
                    pf.Y = pi.Rect.Height - ((float)cv.y - DataRect.Y) * dhi + pi.Rect.Y; // Y is inverted
                    paintCache[i] = pf;
                }
            }
            return paintCache;
        }
        */
        private int comparisonLength;
        int IComparable.CompareTo(object obj)
        {
            comparisonLength = ((Contour)obj).Vertices.Length;
            if (Vertices.Length > comparisonLength)
                return 1;
            if (Vertices.Length == comparisonLength)
                return 0;
            return -1;
        }
    }
}
