package Utilities;

public class Bezier extends CurveBase{

    @Override
    public void setCoefs(Point3D p0, Point3D p1, Point3D p2, Point3D p3) {
        c[3][0] = -p0.X + (3*p1.X) - (3 *p2.X) +p3.X;
        c[2][0] = (3 * p0.X) - (6*p1.X) + (3 * p2.X);
        c[1][0] = -(3 * p0.X) + (3 * p1.X);
        c[0][0] = p0.X;
        c[3][1] = -p0.Y + (3*p1.Y) - (3 *p2.Y) +p3.Y;
        c[2][1] = (3 * p0.Y) - (6*p1.Y) + (3 * p2.Y);
        c[1][1] = -(3 * p0.Y) + (3 * p1.Y);
        c[0][1] = p0.Y;
        c[3][2] = -p0.Z + (3*p1.Z) - (3 *p2.Z) +p3.Z;
        c[2][2] = (3 * p0.Z) - (6*p1.Z) + (3 * p2.Z);
        c[1][2] = -(3 * p0.Z) + (3 * p1.Z);
        c[0][2] = p0.Z;
    }
}
