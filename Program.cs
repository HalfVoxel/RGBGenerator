#region ©
/*
 * Copyright © Jozsef Fejes, http://joco.name/
 * 
 * This code is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#endregion

/*
 * General optimization notes:
 * - We try hard not to give any job for the GC. We allocate all objects at startup and only use structs later.
 * - To be very fast, we only use arrays as collections, and not even a foreach.
 * - Parallel processing is tricky. Smaller chunks mean more balance but more overhead. Our parallel operations
 *   are very small and we keep the queue balanced, so we aim for the biggest possible chunk size.
 */

using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;

namespace FejesJoco.Tools.RGBGenerator
{
    class Program
    {
        #region settings
        /// <summary>
        /// Number of colors per channel.
        /// </summary>
        static int NumColors;

        /// <summary>
        /// Width of the image.
        /// </summary>
        static int Width;

        /// <summary>
        /// Height of the image.
        /// </summary>
        static int Height;

        /// <summary>
        /// First pixel X coordinate.
        /// </summary>
        static int StartX;

        /// <summary>
        /// First pixel Y coordinate.
        /// </summary>
        static int StartY;

        /// <summary>
        /// Number of image frames to save.
        /// </summary>
        static int NumFrames;

        /// <summary>
        /// Random generator, only used during precalculations in a deterministic way. The same seed awlays results in the same image.
        /// </summary>
        static Random RndGen;

        /// <summary>
        /// Available neighbor X coordinate differences (-1,0,+1).
        /// </summary>
        static int[] NeighX;

        /// <summary>
        /// Available neighbor Y coordinate differences (-1,0,+1).
        /// </summary>
        static int[] NeighY;

        /// <summary>
        /// The chosen color sorting implementation.
        /// </summary>
        static IComparer<RGB> Sorter;

        /// <summary>
        /// The chosen algorithm implementation.
        /// </summary>
        static AlgorithmBase Algorithm;

        /// <summary>
        /// Prints the command line help.
        /// </summary>
        static void printArgsHelp()
        {
            Console.WriteLine();
            Console.WriteLine("You have to run the program like this:");
            Console.WriteLine("{0} [colors] [width] [height] [startx] [starty] [frames] [seed] [neighbors] [sorting] [algo]",
                AppDomain.CurrentDomain.FriendlyName);
            Console.WriteLine();
            Console.WriteLine("For a quick start, try this:");
            Console.WriteLine("{0} 64 512 512 256 256 5 0 11111111 rnd one",
                AppDomain.CurrentDomain.FriendlyName);
            Console.WriteLine();
            Console.WriteLine("Explanation:");
            Console.WriteLine("[colors]: number of colors per channel, at most 256, must be a power of 2");
            Console.WriteLine("[width], [height]: dimensions of the image, each must be a power of 2");
            Console.WriteLine("[startx], [starty]: coordinates of the first pixel");
            Console.WriteLine("[frames]: number of frames to generate, must be positive");
            Console.WriteLine("[seed]: random seed, try 0 for truly random");
            Console.WriteLine("[neighbors]: specify eight 1's or 0's to allow movements in these directions");
            Console.WriteLine("[sorting]: color ordering, can be 'rnd' or 'hue-N' (N=0..360)");
            Console.WriteLine("[algo]: can be 'one' or 'avg' or 'avgsq'");
            Console.WriteLine();
        }

        /// <summary>
        /// Parses the command line arguments.
        /// </summary>
        static bool parseArgs(string[] args)
        {
            if (args.Length != 10)
            {
                Console.WriteLine("There must be exactly 10 arguments given!");
                return false;
            }

            // generate 2^0..2^24 for easier checking
            var twopows = new List<string>();
            var p = 1;
            for (var i = 0; i <= 24; i++)
            {
                twopows.Add(p.ToString());
                p *= 2;
            }

            // [colors]
            if (!twopows.Contains(args[0]) || !int.TryParse(args[0], out NumColors) || NumColors > 256)
            {
                Console.WriteLine("[colors] is an invalid number");
                return false;
            }

            // [width] and [height]
            if (!twopows.Contains(args[1]) || !int.TryParse(args[1], out Width))
            {
                Console.WriteLine("[width] is an invalid number");
                return false;
            }
            if (!twopows.Contains(args[2]) || !int.TryParse(args[2], out Height))
            {
                Console.WriteLine("[height] is an invalid number");
                return false;
            }
            if ((long)Width * Height != NumColors * NumColors * NumColors)
            {
                Console.WriteLine("[width]*[height] must be equal to [colors]*[colors]*[colors]");
                return false;
            }

            // [startx], [starty]
            if (!int.TryParse(args[3], out StartX))
            {
                Console.WriteLine("[startx] is an invalid number");
                return false;
            }
            if (StartX < 0 || StartX >= Width)
            {
                Console.WriteLine("[startx] is out of bounds");
                return false;
            }
            if (!int.TryParse(args[4], out StartY))
            {
                Console.WriteLine("[starty] is an invalid number");
                return false;
            }
            if (StartY < 0 || StartY >= Height)
            {
                Console.WriteLine("[starty] is out of bounds");
                return false;
            }

            // [frames]
            if (!int.TryParse(args[5], out NumFrames))
            {
                Console.WriteLine("[frames] is an invalid number");
                return false;
            }
            if (NumFrames < 1)
            {
                Console.WriteLine("[frames] must be positive");
                return false;
            }

            // [seed]
            int seed;
            if (!int.TryParse(args[6], out seed))
            {
                Console.WriteLine("[seed] is an invalid number");
                return false;
            }
            RndGen = seed == 0 ? new Random() : new Random(seed);

            // [neighbors]
            if (!Regex.IsMatch(args[7], "^[01]{8}$"))
            {
                Console.WriteLine("[neighbors] is not given according to the rules");
                return false;
            }
            var nx = new[] { -1, 0, 1, -1, 1, -1, 0, 1 };
            var ny = new[] { -1, -1, -1, 0, 0, 1, 1, 1 };
            NeighX = Enumerable.Range(0, 8).Where(i => args[7][i] == '1').Select(i => nx[i]).ToArray();
            NeighY = Enumerable.Range(0, 8).Where(i => args[7][i] == '1').Select(i => ny[i]).ToArray();

            // [sorting]
            if (args[8] == "rnd")
                Sorter = new RandomComparer();
            else if (args[8].StartsWith("hue-"))
            {
                int hueshift;
                if (!int.TryParse(args[8].Substring(4), out hueshift) || hueshift < 0 || hueshift > 360)
                {
                    Console.WriteLine("[sorting] has an invalid hue parameter");
                    return false;
                }
                Sorter = new HueComparer(hueshift);
            }
            else
            {
                Console.WriteLine("[sorting] is not one of the allowed values");
                return false;
            }

            // [algo]
            if (args[9] == "one")
                Algorithm = new OneNeighborSqAlgorithm();
            else if (args[9] == "avg")
                Algorithm = new AverageNeighborAlgorithm();
            else if (args[9] == "avgsq")
                Algorithm = new AverageNeighborSqAlgorithm();
            else
            {
                Console.WriteLine("[algo] is not one of the allowed values");
                return false;
            }

            Console.WriteLine("Command line arguments are accepted");
            return true;
        }
        #endregion

        #region color sorting algorithms
        /// <summary>
        /// Totally random (but detereministic) color sorting.
        /// </summary>
        class RandomComparer : IComparer<RGB>
        {
            public int Compare(RGB x, RGB y)
            {
                return RndGen.Next(11) - 5;
            }
        }

        /// <summary>
        /// Compares by hue first, then by brightness, and finally random.
        /// </summary>
        class HueComparer : IComparer<RGB>
        {
            private int shift;

            public HueComparer(int shift)
            {
                this.shift = shift;
            }

            public int Compare(RGB x, RGB y)
            {
                var xc = x.ToColor();
                var yc = y.ToColor();
                var c = ((xc.GetHue() + shift) % 360).CompareTo((yc.GetHue() + shift) % 360);
                if (c == 0)
                    c = xc.GetBrightness().CompareTo(yc.GetBrightness());
                if (c == 0)
                    c = RndGen.Next(11) - 5;
                return c;
            }
        }
        #endregion

        #region pixel management
        /// <summary>
        /// Represents an 8-bit RGB color. Faster and smaller than the framework version, <see cref="Color"/>.
        /// </summary>
        struct RGB
        {
            public byte R;
            public byte G;
            public byte B;

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
        }

        /// <summary>
        /// Represents a pixel in the big image.
        /// </summary>
        class Pixel
        {
            /// <summary>
            /// Is this pixel empty? If so, it has no color.
            /// </summary>
            public bool Empty;

            /// <summary>
            /// Color of this pixel.
            /// </summary>
            public RGB Color;

            /// <summary>
            /// Index of this pixel in the queue (<see cref="PixelList"/>), or -1 if it's not queued.
            /// </summary>
            public int QueueIndex;

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
        struct PixelWithValue : IComparable<PixelWithValue>
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

        /// <summary>
        /// Represents a pixel queue. It's a blend of <see cref="List{T}"/> and <see cref="Dictionary{Tk,Tv}"/> functionality. It allows very quick,
        /// indexed traversal (it just exposes a simple array). It supports O(1) lookups (every pixel contains it's own index in this array). Adding
        /// and removal are also O(1) because we don't usually reallocate the array.
        /// </summary>
        class PixelQueue
        {
            /// <summary>
            /// Array of pixels. Some entries may be null.
            /// </summary>
            public Pixel[] Pixels;

            /// <summary>
            /// Number of entries (non-null elements).
            /// </summary>
            public int Count;

            /// <summary>
            /// The first index where new entries may be added. We don't go past this when reading all elements.
            /// </summary>
            public int UsedUntil;

            public PixelQueue()
            {
                Pixels = new Pixel[1024];
            }

            /// <summary>
            /// Add a pixel to the queue.
            /// </summary>
            public virtual void Add(Pixel p)
            {
                Debug.Assert(p.QueueIndex == -1);
                if (UsedUntil == Pixels.Length)
                    Array.Resize(ref Pixels, Pixels.Length * 2);
                Pixels[UsedUntil] = p;
                p.QueueIndex = UsedUntil;
                UsedUntil++;
                Count++;
            }

            /// <summary>
            /// Remove a pixel from the queue.
            /// </summary>
            public virtual void Remove(Pixel p)
            {
                Debug.Assert(p.QueueIndex > -1);
                Pixels[p.QueueIndex] = null;
                p.QueueIndex = -1;
                Count--;
            }

            /// <summary>
            /// Readds a pixel, while maintaining its associated data in <see cref="PixelQueue{T}"/>.
            /// </summary>
            public virtual void Readd(Pixel p)
            {
                Remove(p);
                Add(p);
            }

            /// <summary>
            /// Compression. Fills up gaps by moving elements to the front. (Gaps add a little overhead to iteration.)
            /// </summary>
            public void Compress()
            {
                // we allow at most 5% to be wasted
                if ((double)UsedUntil / Count < 1.05)
                    return;
                UsedUntil = 0;
                for (var i = 0; UsedUntil < Count; i++)
                {
                    if (Pixels[i] != null)
                    {
                        Readd(Pixels[i]);
                    }
                }
            }
        }

        /// <summary>
        /// The only thing it adds to its ancestor <see cref="PixelQueue"/> is that it also maintains an arbitrary data object for every pixel. It
        /// uses the same array indexing, uses structs and doesn't do cleanups.
        /// </summary>
        class PixelQueue<T> : PixelQueue
            where T : struct
        {
            public T[] Data;

            public PixelQueue()
            {
                Data = new T[Pixels.Length];
            }

            public override void Add(Pixel p)
            {
                base.Add(p);
                // we need to maintain the same array size
                if (Pixels.Length != Data.Length)
                    Array.Resize(ref Data, Pixels.Length);
            }

            public override void Readd(Pixel p)
            {
                // maintain data when moving a pixel
                var data = Data[p.QueueIndex];
                base.Readd(p);
                Data[p.QueueIndex] = data;
            }
        }
        #endregion

        #region pixel placing algorithms
        /// <summary>
        /// Base class of algorithms.
        /// </summary>
        abstract class AlgorithmBase
        {
            /// <summary>
            /// Every implementation has a pixel queue.
            /// </summary>
            public abstract PixelQueue Queue { get; }

            /// <summary>
            /// Places the given color on the image.
            /// </summary>
            public void Place(RGB c)
            {
                // find the next coordinates
                Pixel p;
                if (Queue.Count == 0)
                    p = Image[StartY * Width + StartX];
                else
                    p = placeImpl(c);

                // put the pixel where it belongs
                Debug.Assert(p.Empty);
                p.Empty = false;
                p.Color = c;

                // adjust the queue
                changeQueue(p);
            }

            /// <summary>
            /// Places the given color on the image. Can assume that the queue is not empty.
            /// </summary>
            protected abstract Pixel placeImpl(RGB c);

            /// <summary>
            /// Adjusts the queue after placing the given pixel.
            /// </summary>
            protected abstract void changeQueue(Pixel p);
        }

        /// <summary>
        /// The queue contains filled pixels, which have at least one empty neighbor. In each step, we find the closest match to the new color in the
        /// queue. Then we place the new color into a random empty neighbor. We use a squared difference metric.
        /// </summary>
        class OneNeighborSqAlgorithm : AlgorithmBase
        {
            PixelQueue queue = new PixelQueue();

            public override PixelQueue Queue
            {
                get { return queue; }
            }

            protected override Pixel placeImpl(RGB c)
            {
                // find the best pixel with parallel processing
                var q = queue.Pixels;
                var best = Partitioner.Create(0, queue.UsedUntil, Math.Max(256, queue.UsedUntil / Threads)).AsParallel().Min(range =>
                {
                    var bestdiff = int.MaxValue;
                    Pixel bestpixel = null;
                    for (var i = range.Item1; i < range.Item2; i++)
                    {
                        var qp = q[i];
                        if (qp != null)
                        {
                            var rd = (int)qp.Color.R - c.R;
                            var gd = (int)qp.Color.G - c.G;
                            var bd = (int)qp.Color.B - c.B;
                            var diff = rd * rd + gd * gd + bd * bd;
                            // we have to use the same comparison as PixelWithValue!
                            if (diff < bestdiff || (diff == bestdiff && qp.Weight < bestpixel.Weight))
                            {
                                bestdiff = diff;
                                bestpixel = qp;
                            }
                        }
                    }
                    return new PixelWithValue
                    {
                        Pixel = bestpixel,
                        Value = bestdiff
                    };
                }).Pixel;

                // select a deterministically random empty neighbor
                var shift = best.Weight % best.Neighbors.Length;
                for (var i = 0; i < best.Neighbors.Length; i++)
                {
                    var bestn = best.Neighbors[(i + shift) % best.Neighbors.Length];
                    if (bestn.Empty)
                        return bestn;
                }

                throw new Exception("not reached");
            }

            protected override void changeQueue(Pixel p)
            {
                for (var i = 0; i < p.Neighbors.Length; i++)
                {
                    var np = p.Neighbors[i];
                    if (np.Empty)
                    {
                        // if p has an empty neighbor, it belongs in the queue
                        if (p.QueueIndex == -1)
                            queue.Add(p);
                    }
                    else
                    {
                        // p has a filled neighbor, and p was just filled, so that neighbor may not belong to the queue anymore
                        var stillok = false;
                        for (var j = 0; j < np.Neighbors.Length; j++)
                            if (np.Neighbors[j].Empty)
                            {
                                stillok = true;
                                break;
                            }
                        if (!stillok)
                            queue.Remove(np);
                    }
                }
            }
        }

        /// <summary>
        /// The queue contains empty pixels which have at least one filled neighbor. For every pixel in the queue, we calculate the average color of
        /// its neighbors. In each step, we find the average that matches the new color the most. Uses a linear difference metric. It gives a blurred
        /// effect.
        /// </summary>
        class AverageNeighborAlgorithm : AlgorithmBase
        {
            PixelQueue<RGB> queue = new PixelQueue<RGB>();

            public override PixelQueue Queue
            {
                get { return queue; }
            }

            protected override Pixel placeImpl(RGB c)
            {
                // find the best pixel with parallel processing
                var q = queue.Pixels;
                var best = Partitioner.Create(0, queue.UsedUntil, Math.Max(256, queue.UsedUntil / Threads)).AsParallel().Min(range =>
                {
                    var bestdiff = int.MaxValue;
                    Pixel bestpixel = null;
                    for (var i = range.Item1; i < range.Item2; i++)
                    {
                        var qp = q[i];
                        if (qp != null)
                        {
                            var avg = queue.Data[qp.QueueIndex];
                            var rd = (int)avg.R - c.R;
                            var gd = (int)avg.G - c.G;
                            var bd = (int)avg.B - c.B;
                            var diff = rd * rd + gd * gd + bd * bd;
                            // we have to use the same comparison as PixelWithValue!
                            if (diff < bestdiff || (diff == bestdiff && qp.Weight < bestpixel.Weight))
                            {
                                bestdiff = diff;
                                bestpixel = qp;
                            }
                        }
                    }
                    return new PixelWithValue
                    {
                        Pixel = bestpixel,
                        Value = bestdiff
                    };
                }).Pixel;

                // found the pixel, return it
                queue.Remove(best);
                return best;
            }

            protected override void changeQueue(Pixel p)
            {
                // recalculate the neighbors
                for (var i = 0; i < p.Neighbors.Length; i++)
                {
                    var np = p.Neighbors[i];
                    if (np.Empty)
                    {
                        int r = 0, g = 0, b = 0, n = 0;
                        for (var j = 0; j < np.Neighbors.Length; j++)
                        {
                            var nnp = np.Neighbors[j];
                            if (!nnp.Empty)
                            {
                                r += nnp.Color.R;
                                g += nnp.Color.G;
                                b += nnp.Color.B;
                                n++;
                            }
                        }
                        var avg = new RGB
                        {
                            R = (byte)(r / n),
                            G = (byte)(g / n),
                            B = (byte)(b / n)
                        };
                        if (np.QueueIndex == -1)
                            queue.Add(np);
                        queue.Data[np.QueueIndex] = avg;
                    }
                }
            }
        }

        /// <summary>
        /// Almost like <see cref="AverageNeighborAlgorithm"/> (also lots of copypaste). The difference is that it uses a squared difference. This
        /// leaves lots of gaps between differently colored branches, so it grows like a coral, which is very cool. The downside is that the queue
        /// gets very big. It also has a blurred effect, but with some more rough edges.
        /// </summary>
        class AverageNeighborSqAlgorithm : AlgorithmBase
        {
            /// <summary>
            /// Let's say we have three neighbors: A, B, C. We want to add a fourth color, D. The squared difference for one channel is
            /// (D-A)^2+(D-B)^2+(D-C)^2. But for this, we need to read every neighbor. We need to make it quick by using precalculation. So we
            /// rearrange the formula like this: 3*D^2+((A^2+B^2+C^2)-2D(A+B+C)). Now we can precalculate the values which depend on the neighbors.
            /// This structure contains these values.
            /// </summary>
            struct AvgInfo
            {
                // sum of neighbors
                public int R, G, B;
                // squared sum of neighbors
                public int RSq, GSq, BSq;
                // number of neighbors
                public int Num;
            }

            PixelQueue<AvgInfo> queue = new PixelQueue<AvgInfo>();

            public override PixelQueue Queue
            {
                get { return queue; }
            }

            protected override Pixel placeImpl(RGB c)
            {
                // find the best pixel with parallel processing
                var q = queue.Pixels;
                var best = Partitioner.Create(0, queue.UsedUntil, Math.Max(256, queue.UsedUntil / Threads)).AsParallel().Min(range =>
                {
                    var bestdiff = int.MaxValue;
                    Pixel bestpixel = null;
                    for (var i = range.Item1; i < range.Item2; i++)
                    {
                        var qp = q[i];
                        if (qp != null)
                        {
                            var avg = queue.Data[qp.QueueIndex];
                            var cr = (int)c.R;
                            var cg = (int)c.G;
                            var cb = (int)c.B;
                            var rd = cr * cr * avg.Num + (avg.RSq - 2 * cr * avg.R);
                            var gd = cg * cg * avg.Num + (avg.GSq - 2 * cg * avg.G);
                            var bd = cb * cb * avg.Num + (avg.BSq - 2 * cb * avg.B);
                            var diff = rd + gd + bd;
                            // we have to use the same comparison as PixelWithValue!
                            if (diff < bestdiff || (diff == bestdiff && qp.Weight < bestpixel.Weight))
                            {
                                bestdiff = diff;
                                bestpixel = qp;
                            }
                        }
                    }
                    return new PixelWithValue
                    {
                        Pixel = bestpixel,
                        Value = bestdiff
                    };
                }).Pixel;

                // found the pixel, return it
                queue.Remove(best);
                return best;
            }

            protected override void changeQueue(Pixel p)
            {
                // recalculate the neighbors
                for (var i = 0; i < p.Neighbors.Length; i++)
                {
                    var np = p.Neighbors[i];
                    if (np.Empty)
                    {
                        int r = 0, g = 0, b = 0, rsq = 0, gsq = 0, bsq = 0, n = 0;
                        for (var j = 0; j < np.Neighbors.Length; j++)
                        {
                            var nnp = np.Neighbors[j];
                            if (!nnp.Empty)
                            {
                                var nr = (int)nnp.Color.R;
                                var ng = (int)nnp.Color.G;
                                var nb = (int)nnp.Color.B;
                                r += nr;
                                g += ng;
                                b += nb;
                                n++;
                                rsq += nr * nr;
                                gsq += ng * ng;
                                bsq += nb * nb;
                            }
                        }
                        var avg = new AvgInfo
                        {
                            R = r,
                            G = g,
                            B = b,
                            RSq = rsq,
                            GSq = gsq,
                            BSq = bsq,
                            Num = n
                        };
                        if (np.QueueIndex == -1)
                            queue.Add(np);
                        queue.Data[np.QueueIndex] = avg;
                    }
                }
            }
        }
        #endregion

        #region main
        /// <summary>
        /// Holds the big image.
        /// </summary>
        static Pixel[] Image;

        /// <summary>
        /// Number of threads, equal to the number of (logical) processors.
        /// </summary>
        static int Threads = Environment.ProcessorCount;

        static void Main(string[] args)
        {
            Console.WriteLine("RGB image generator by Jozsef Fejes");
            Console.WriteLine("Check out my blog for more information: http://joco.name/");
            Console.WriteLine();

            // parse command-line arguments
            if (!parseArgs(args))
            {
                printArgsHelp();
                Console.WriteLine("Press ENTER to exit");
                Console.ReadLine();
                return;
            }

            // we're just about to start
            System.Diagnostics.Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
            var start = DateTime.Now;
            Console.WriteLine("Running the precalculations and eating up your memory...");

            // create every color once and randomize their order
            var colors = new List<RGB>(Width * Height);
            for (var r = 0; r < NumColors; r++)
                for (var g = 0; g < NumColors; g++)
                    for (var b = 0; b < NumColors; b++)
                        colors.Add(new RGB
                        {
                            R = (byte)(r * 255 / (NumColors - 1)),
                            G = (byte)(g * 255 / (NumColors - 1)),
                            B = (byte)(b * 255 / (NumColors - 1))
                        });
            colors.Sort(Sorter);
            Debug.Assert(colors.Count == Width * Height);

            // create the pixels and their unique weights
            Image = new Pixel[Width * Height];
            var weights = new HashSet<int>();
            for (var y = 0; y < Height; y++)
                for (var x = 0; x < Width; x++)
                {
                    int weight;
                    do
                    {
                        weight = RndGen.Next();
                    } while (!weights.Add(weight));
                    Image[y * Width + x] = new Pixel
                    {
#if DEBUG
                        X = x,
                        Y = y,
#endif
                        Empty = true,
                        QueueIndex = -1,
                        Weight = weight
                    };
                }

            // precalculate the neighbors of every pixel
            Debug.Assert(NeighX.Length == NeighY.Length);
            for (var y = 0; y < Height; y++)
                for (var x = 0; x < Width; x++)
                {
                    Image[y * Width + x].Neighbors = Enumerable.Range(0, NeighX.Length).Select(n =>
                    {
                        var y2 = y + NeighY[n];
                        if (y2 < 0 || y2 == Height)
                            return null;
                        var x2 = x + NeighX[n];
                        if (x2 < 0 || x2 == Width)
                            return null;
                        return Image[y2 * Width + x2];
                    }).Where(p => p != null).ToArray();
                }

            // calculate the saving checkpoints in advance
            var checkpoints = Enumerable.Range(1, NumFrames).ToDictionary(i => (long)i * colors.Count / NumFrames - 1, i => i - 1);
            Thread pngthread = null;

            // loop through all colors that we want to place
            for (var i = 0; i < colors.Count; i++)
            {
                // give progress report to the impatient user
                if (i % 1024 == 0)
                {
                    Algorithm.Queue.Compress();
                    Console.WriteLine("{0:P}, queue size {1}", (double)i / Width / Height, Algorithm.Queue.Count);
                }

                // run the algorithm
                Algorithm.Place(colors[i]);

                // save a checkpoint if needed
                int chkidx;
                if (checkpoints.TryGetValue(i, out chkidx))
                {
                    // create the image
                    var img = new Bitmap(Width, Height, PixelFormat.Format24bppRgb);
                    var idata = img.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.WriteOnly, PixelFormat.Format24bppRgb);
                    var ibytes = new byte[idata.Stride * idata.Height];
                    for (var y = 0; y < Height; y++)
                    {
                        for (var x = 0; x < Width; x++)
                        {
                            var c = Image[y * Width + x].Color;
                            ibytes[y * idata.Stride + x * 3 + 2] = c.R;
                            ibytes[y * idata.Stride + x * 3 + 1] = c.G;
                            ibytes[y * idata.Stride + x * 3 + 0] = c.B;
                        }
                    }
                    Marshal.Copy(ibytes, 0, idata.Scan0, ibytes.Length);
                    img.UnlockBits(idata);

                    // png compression uses only one processor, so push it into the background, limiting to one thread at a time is more than enough
                    if (pngthread != null)
                        pngthread.Join();
                    pngthread = new Thread(new ThreadStart(delegate
                    {
                        img.Save(string.Format("result{0:D5}.png", chkidx), ImageFormat.Png);
                        img.Dispose();
                    }));
                    pngthread.Start();
                }
            }
            Debug.Assert(Algorithm.Queue.Count == 0);

            // wait for the final image
            pngthread.Join();

            // check the number of colors to be sure
            var img2 = (Bitmap)Bitmap.FromFile(string.Format("result{0:D5}.png", NumFrames - 1));
            var ch = new HashSet<RGB>();
            for (var y = 0; y < img2.Height; y++)
                for (var x = 0; x < img2.Width; x++)
                {
                    var pix = img2.GetPixel(x, y);
                    if (!ch.Add(new RGB
                    {
                        R = pix.R,
                        G = pix.G,
                        B = pix.B
                    }))
                    {
                        Console.WriteLine("Color {0}/{1}/{2} is added more than once!!!!!!", pix.R, pix.G, pix.B);
                    }
                }
            img2.Dispose();

            // we're done!
            Console.WriteLine("All done! It took this long: {0}", DateTime.Now.Subtract(start));
            Console.WriteLine("Press ENTER to exit");
            Console.ReadLine();
        }
        #endregion
    }
}
