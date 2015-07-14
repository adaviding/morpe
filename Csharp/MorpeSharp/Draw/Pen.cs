using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Morpe.Draw
{
    public struct Pen
    {
        public Color Color;
        public float Width;
        public Pen(Color Color, float Width)
        {
            this.Color = Color;
            this.Width = Width;
        }
    }
}
