namespace Morpe.Numerics.D1
{
    public class Range
    {
        public double Min;
        public double Max;

        public Range(double min, double max)
        {
            this.Min = min;
            this.Max = max;
        }
        
        public double Clamp(double x)
        {
            if (x < this.Min)
                return this.Min;
            
            if (x > this.Max)
                return this.Max;

            return x;
        }

        public Range Clone()
        {
            return (Range)this.MemberwiseClone();
        }
    }
}
