package Utilities;

public class Hermite extends CurveBase{

    @Override
    public void setCoefs(Point3D p0, Point3D p1, Point3D p2, Point3D p3) {
        Point3D tangent1 = new Point3D(p1.X-p0.X, p1.Y-p0.Y,p1.Z-p0.Z);
        Point3D tangent2 = new Point3D(p2.X-p3.X, p2.Y-p3.Y,p2.Z-p3.Z);


        c[3][0] = (2 *p0.X) - (2*p3.X) + tangent1.X + tangent2.X;
        c[2][0] = (-3 * p0.X) + (3*p3.X) - 2 * tangent1.X - tangent2.X;
        c[1][0] = tangent1.X;
        c[0][0] = p0.X;
        c[3][1] = (2 *p0.Y) - (2*p3.Y) + tangent1.Y + tangent2.Y;
        c[2][1] = (-3 * p0.Y) + (3*p3.Y) - 2 * tangent1.Y - tangent2.Y;
        c[1][1] = tangent1.Y;
        c[0][1] = p0.Y;
        c[3][2] = (2 *p0.Z) - (2*p3.Z) + tangent1.Z + tangent2.Z;
        c[2][2] = (-3 * p0.Z) + (3*p3.Z) - 2 * tangent1.Z - tangent2.Z;
        c[1][2] = tangent1.Z;
        c[0][2] = p0.Z;
    }
}
