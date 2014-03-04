using System;
using System.Drawing;
using System.Drawing.Imaging;

namespace FejesJoco.Tools.RGBGenerator
{
	/// <summary>
	/// Represents an 8-bit RGB color. Faster and smaller than the framework version, <see cref="Color"/>.
	/// </summary>
	public struct RGB
	{
		public byte R;
		public byte G;
		public byte B;

		public RGB ( byte R, byte G, byte B ) {
			this.R = R;
			this.G = G;
			this.B = B;
		}

		public RGB ( int R, int G, int B ) {
			this.R = (byte)R;
			this.G = (byte)G;
			this.B = (byte)B;
		}

		public Color ToColor()
		{
			return Color.FromArgb(R, G, B);
		}

		public override int GetHashCode()
		{
			return R * 256 * 256 + G * 256 + B;
		}

		public override bool Equals(object obj)
		{
			var that = (RGB)obj;
			return this.GetHashCode() == that.GetHashCode();
		}

		public override string ToString ()
		{
			return "(" + R + ", " + G + ", " + B + ")";
		}
	}
}

