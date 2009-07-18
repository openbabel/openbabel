#if !MONO
using OpenBabel;

namespace DotBabel
{
    using System.Windows.Media.Media3D;

    public static class MathExtensions
    {
        public static Vector3D ToCSVector(this OBVector3 vec)
        {
            return new Vector3D(vec.x(), vec.x(), vec.z());
        }

        public static Point3D ToCSPoint(this OBVector3 vec)
        {
            return new Point3D(vec.x(), vec.x(), vec.z());
        }
    }
}
#endif