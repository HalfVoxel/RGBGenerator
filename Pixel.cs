using System;

namespace FejesJoco.Tools.RGBGenerator
{
	/// <summary>
	/// Represents a pixel in the big image.
	/// </summary>
	public class Pixel
	{
		/// <summary>
		/// Is this pixel empty? If so, it has no color.
		/// </summary>
		public bool Empty;

		//public uint QueueScore;

		public bool inQueue;

		public byte nonEmptyNeigh = 0;
		public int block;
		//public int indexInBlock;


		public RGB avg;
		/// <summary>
		/// Color of this pixel.
		/// </summary>
		public RGB Color;

		/// <summary>
		/// Index of this pixel in the queue (<see cref="PixelList"/>), or -1 if it's not queued.
		/// </summary>
		public int QueueIndex {
			get { return 0; }
			set {}
		}

		/// <summary>
		/// Precalculated array of neighbor pixels.
		/// </summary>
		public Pixel[] Neighbors;

		/// <summary>
		/// A unique weight of thix pixel. Used in comparisons when calculated values are equal. Needs to be randomly distributed, otherwise the
		/// picture cannot grow equally in every direction.
		/// </summary>
		public int Weight;

		#if DEBUG
		// we don't need to know these, just for debugging
		public int X;
		public int Y;
		public override string ToString()
		{
			return X + ";" + Y;
		}
		#endif
	}

	/// <summary>
	/// Represents a pixel with a value. Used to sort and merge parallel run results.
	/// </summary>
	public struct PixelWithValue : IComparable<PixelWithValue>
	{
		public Pixel Pixel;
		public int Value;

		public int CompareTo(PixelWithValue that)
		{
			// a parallel run may have no result at all, in that case we prefer the one that has a value
			if (this.Pixel == null && that.Pixel == null)
				return 0;
			if (this.Pixel == null)
				return 1;
			if (that.Pixel == null)
				return -1;

			// compare the values, or use the weight if they're equal
			var c = this.Value.CompareTo(that.Value);
			if (c == 0)
				c = this.Pixel.Weight.CompareTo(that.Pixel.Weight);
			return c;
		}
	}
}

