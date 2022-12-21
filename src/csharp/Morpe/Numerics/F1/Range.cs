namespace Morpe.Numerics.F1
{
    public class Range
    {
        public float Min;
        public float Max;

        public Range(float min, float max)
        {
            this.Min = min;
            this.Max = max;
        }

        public float Clamp(float x)
        {
            if (x < this.Min)
                return this.Min;
            
            if (x > this.Max)
                return this.Max;

            return x;
        }
    }
}
