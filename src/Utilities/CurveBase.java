package Utilities;


import Utilities.ImageUtilities;
import LineBase.LineBase;
import LineBase.Lines;


import java.io.IOException;
import java.util.ArrayList;

public abstract class CurveBase {

    // -- inner class to represent a 3D point
    public static class Point3D {
        public double X, Y, Z;

        public Point3D() {
            this.X = 0;
            this.Y = 0;
            this.Z = 0;
        }

        public Point3D(double X, double Y, double Z) {
            this.X = X;
            this.Y = Y;
            this.Z = Z;
        }
    }

    double c[][] = new double[4][3];

    // -- four control points are provided in order. P0 and P3 are the end points of the
    //    curve, P1 is the starting tangent for Hermite, the first control point for Bezier.
    //    P2 is the ending tangent for Hermite, the second control point for Bezier.
    //    For a Hermite curve the function will need to create the two derivatives (P1 - P0) and (P2 - P3).
    //    For a Bezier curve use the points as given in the order given.
    //    For a Lagrange curve don't use this function at all as it requires a list of points.
    public abstract void setCoefs(Point3D p0, Point3D p1, Point3D p2, Point3D p3);

    // -- The argument is the curve parameter 0 <= u <= 1.0
    //    Used for Hermite and Bezier, not used for Lagrange
    public Point3D GeneratePoint(double u) {
        //                        2        3       3    2
        // -- c0 + c1 * u + c2 * u + c3 * u  or  au + bu + cu + D in "standard form"
        Point3D p = new Point3D((((c[3][0] * u) + c[2][0]) * u + c[1][0]) * u + c[0][0],  // -- X
                (((c[3][1] * u) + c[2][1]) * u + c[1][1]) * u + c[0][1],  // -- Y
                (((c[3][2] * u) + c[2][2]) * u + c[1][2]) * u + c[0][2]); // -- Z
        return p;
    }

    public static void main(String[] args) {
        int framebuffer[][][] = new int[3][500][500];

        { // -- render a Bezier curve into the frame buffer
            Bezier bezier = new Bezier();

            ArrayList<Point3D> points = new ArrayList<Point3D>();

            // -- scan convert first section of curve into the point list
            Point3D p0 = new Point3D(25, 25, 0);
            Point3D p1 = new Point3D(50, 200, 0);
            Point3D p2 = new Point3D(200, 200, 0);
            Point3D p3 = new Point3D(225, 25, 0);
            // -- set Bezier curve parameters
            //    start, first knot, second knot, end points
            bezier.setCoefs(p0, p1, p2, p3);
            // -- generate Bezier curve points
            for (double u = 0; u <= 1; u += 0.05) {
                points.add(bezier.GeneratePoint(u));
            }


            // -- render the curve into the frame buffer using piecewise linear Bresenham scan conversion
            Lines lb = new Lines();
            Point3D point0 = null;
            Point3D point1 = null;
            LineBase.RGBColor c0 = new LineBase.RGBColor(255, 0, 0);
            LineBase.RGBColor c1 = new LineBase.RGBColor(255, 255, 0);
            float redstep = ((float) c1.R - (float) c0.R) / points.size();
            float grnstep = ((float) c1.G - (float) c0.G) / points.size();
            float blustep = ((float) c1.B - (float) c0.B) / points.size();
            float r = c0.R, g = c0.G, b = c0.B;
            for (int i = 0; i < points.size() - 1; ++i) {
                point0 = points.get(i);
                point1 = points.get(i + 1);
                lb.BresenhamFormRGB((int) point0.X, (int) point0.Y, (int) point1.X, (int) point1.Y, framebuffer,
                        new LineBase.RGBColor((int) r, (int) g, (int) b), new LineBase.RGBColor((int) r, (int) g, (int) b));
                r += redstep;
                g += grnstep;
                b += blustep;
            }
            // -- close out the curve (floating point calculations don't reach the end point)
            lb.BresenhamFormRGB((int) point1.X, (int) point1.Y, (int) p3.X, (int) p3.Y, framebuffer,
                    new LineBase.RGBColor((int) r, (int) g, (int) b), new LineBase.RGBColor((int) c1.R, (int) c1.G, (int) c1.B));
        }
        { // -- render a Hermite curve into the frame buffer
            Hermite hermite = new Hermite();

            ArrayList<Point3D> points = new ArrayList<Point3D>();

            // -- scan convert first section of curve into the point list
            Point3D p0 = new Point3D(76, 304, 0);
            Point3D p1 = new Point3D(443, 299, 0);
            Point3D p2 = new Point3D(443, 289, 0);
            Point3D p3 = new Point3D(416, 188, 0);
            // -- set Hermite curve parameters
            //       start point, endpoint, start derivative, end derivative
            //       if specified as points, derivatives need to subtract start and end
            //       points respectively to make them vectors
            hermite.setCoefs(p0, p1, p2, p3);
            // -- generate Hermite curve points
            for (double u = 0; u <= 1; u += 0.05) {
                points.add(hermite.GeneratePoint(u));
            }

            // -- scan convert next section of curve into the point list
            //    ensure G1 continuity by setting the starting point and tangent to the
            //    ending point and tangent of the previous segment
            p0 = p3;
            p1 = p2;
            p2 = new Point3D(310, 73, 0);
            p3 = new Point3D(120, 185, 0);

            // -- set Hermite curve parameters
            //       start point, endpoint, start derivative, end derivative
            //       if specified as points, derivatives need to subtract start and end
            //       points respectively to make them vectors
            hermite.setCoefs(p0, p1, p2, p3);
            // -- generate Hermite curve points
            for (double u = 0; u <= 1; u += 0.05) {
                points.add(hermite.GeneratePoint(u));
            }


            // -- render the curve into the frame buffer using piecewise linear Bresenham scan conversion
            Point3D point0 = null;
            Point3D point1 = null;
            LineBase.RGBColor c0 = new LineBase.RGBColor(0, 255, 0);
            LineBase.RGBColor c1 = new LineBase.RGBColor(255, 0, 255);
            float redstep = ((float)c1.R - (float)c0.R) / points.size();
            float grnstep = ((float)c1.G - (float)c0.G) / points.size();
            float blustep = ((float)c1.B - (float)c0.B) / points.size();
            float r = c0.R;
            float g = c0.G;
            float b = c0.B;
            LineBase lb = new Lines();
            for (int i = 0; i < points.size() - 1; ++i) {
                point0 = points.get(i);
                point1 = points.get(i + 1);
                lb.BresenhamFormRGB((int)point0.X, (int)point0.Y, (int)point1.X, (int)point1.Y, framebuffer,
                        new LineBase.RGBColor((int)r, (int)g, (int)b), new LineBase.RGBColor((int)r, (int)g, (int)b));
                r += redstep;
                g += grnstep;
                b += blustep;
            }
            // -- close out the curve (floating point calculations don't reach the end point)
            lb.BresenhamFormRGB((int)point1.X, (int)point1.Y, (int)p3.X, (int)p3.Y, framebuffer,
                    new LineBase.RGBColor((int)r, (int)g, (int)b), c1);
        }

        try {
            ImageUtilities.ImageWriteRGB(framebuffer, "curves.png");
        } catch (IOException e) {
            System.out.println(e);
        }
    }


}