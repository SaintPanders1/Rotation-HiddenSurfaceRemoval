package Cube;

import LineBase.LineBase;
import LineBase.Lines;
import Utilities.CurveBase;
import Utilities.CurveBase.Point3D;
import Utilities.Bezier;
import Utilities.Hermite;
import Utilities.ImageUtilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;

public class CubeMain {

    //faster version of Matrix Multiplication
    public static double[][] WinogradMethod(double[][] matrix1, double[][] matrix2) throws IllegalArgumentException
    {
        if (matrix1[0].length != matrix2.length)
        {
            throw new IllegalArgumentException("Matrices are not compatible");
        }
        double[][] result = new double[matrix1.length][matrix2[0].length];
        int b = matrix1[0].length;
        int d = b/2;
        double[] rowFactor = new double[matrix1.length];
        double[] columnFactor = new double[matrix2[0].length];
        for(int i = 0; i < matrix1.length; i++)
        {
            rowFactor[i] = (matrix1[i][0] * matrix1[i][1]);
            for(int j = 1; j < d; j++)
            {
                rowFactor[i] += matrix1[i][2 * j] * matrix1[i][2 * j + 1];
            }
        }

        for(int i = 0; i < matrix2[0].length; i++)
        {
            columnFactor[i] = matrix2[0][i] * matrix2[1][i];
            for(int j =1; j < d; j++)
            {
                columnFactor[i] += matrix2[2*j][i] * matrix2[2*j+1][i];
            }
        }

        for(int i = 0; i < rowFactor.length; i++)
        {
            for(int j = 0; j < columnFactor.length; j++)
            {
                result[i][j] = -rowFactor[i] - columnFactor[j];
                for(int k = 0; k < d; k++)
                {
                    result[i][j] = result[i][j] + (matrix1[i][2*k] + matrix2[2*k+1][j]) * (matrix1[i][2*k+1] + matrix2[2*k][j]);
                }
            }
        }

        if(2 * (matrix1[0].length /2) != matrix1[0].length)
        {
            for(int i = 0; i < matrix1.length; i++)
            {
                for(int j = 0; j < matrix2[0].length; j++)
                {
                    result[i][j] = result[i][j] + matrix1[i][b] * matrix2[b][j];
                }
            }
        }
        return result;
    }

//    public static double[][] MatrixMultiplication(double[][] matrix1, double[][] matrix2)throws IllegalArgumentException
//    {
//        if(matrix1[0].length != matrix2.length)
//        {
//            throw new IllegalArgumentException("Matrices do not match");
//        }
//        double[][] result = new double[matrix1.length][matrix2[0].length];
//        for(int i = 0; i < matrix2.length; i++)
//        {
//            for(int j = 0; j <matrix2[0].length; j++ )
//            {
//                for(int k = 0; k < matrix1[0].length; k++)
//                {
//                    result[i][j] = result[i][j] + matrix1[j][k] * matrix2[i][j];
//                }
//            }
//        }
//        return result;
//    }

    public static double[][] scale(double x, double y, double z)
    {
        return new double[][]  {{x,0.0,0.0,0.0},{0.0,y,0.0,0.0}, {0.0,0.0,z,0.0}, {0.0,0.0,0.0,1.0}};
    }
    public static double[][] translate(double x, double y, double z)
    {
        return new double[][] {{1.0,0.0,0.0,x},{0.0,1.0,0.0,y}, {0.0,0.0,1.0,z}, {0.0,0.0,0.0,1.0}};
    }
    public static double[][] rotateX(double angle)
    {
        angle = Math.toRadians(angle);
        return new double[][] {{1.0,0.0,0.0,0.0},{0.0,Math.cos(angle),-Math.sin(angle),0.0}, {0.0,Math.sin(angle),Math.cos(angle),0.0}, {0.0,0.0,0.0,1.0}};
    }
    public static double[][] rotateY(double angle)
    {
        angle = Math.toRadians(angle);
        return new double[][] {{Math.cos(angle),0.0,Math.sin(angle),0.0},{0.0,1.0,0.0,0.0}, {-Math.sin(angle),0.0,Math.cos(angle),0.0}, {0.0,0.0,0.0,1.0}};
    }
    public static double[][] rotateZ(double angle)
    {
        angle = Math.toRadians(angle);
        return new double[][] {{Math.cos(angle),-Math.sin(angle),0.0,0.0},{Math.sin(angle),Math.cos(angle),0.0,0.0}, {0.0,0.0,1.0,0.0}, {0.0,0.0,0.0,1.0}};
    }

    //public static double[][] ArbitraryAxisRotation(int angle, )

    public static double vectorMag(double x, double y, double z)
    {
        return Math.sqrt(Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2));
    }

    public static double[][] rotation(double[][] masterMatrix, Point3D angle, double angleDegrees)
    {
        Point3D centroid = new Point3D(masterMatrix[0][masterMatrix[0].length-1],masterMatrix[1][masterMatrix[1].length-1],masterMatrix[2][masterMatrix[2].length-1]);
        Point3D alpha = new Point3D(angle.X / vectorMag(angle.X,angle.Y,angle.Z), angle.Y / vectorMag(angle.X,angle.Y,angle.Z), angle.Z / vectorMag(angle.X,angle.Y,angle.Z));
        double[][] translateCenter = translate(-centroid.X, -centroid.Y, -centroid.Z);
        double d = Math.sqrt(Math.pow(alpha.Y,2) + Math.pow(alpha.Z,2));
        double[][] Rx = {
                {1.0, 0, 0, 0},
                {0, alpha.Z / d, -alpha.Y / d, 0},
                {0, alpha.Y / d, alpha.Z / d, 0},
                {0, 0, 0, 1}};
        double[][] Ry = {
                {d, 0, -alpha.X, 0},
                {0, 1, 0, 0},
                {alpha.X, 0, d, 0},
                {0, 0, 0, 1}};
        double[][] notRx = {{1.0, 0, 0, 0},
                {0, alpha.Z / d, alpha.Y / d, 0},
                {0, -alpha.Y / d, alpha.Z / d, 0},
                {0, 0, 0, 1}};
        double[][] notRy = {
                {d, 0, alpha.X, 0},
                {0, 1, 0, 0},
                {-alpha.X, 0, d, 0},
                {0, 0, 0, 1}};
        double[][] translateCenterBack = translate(centroid.X, centroid.Y, centroid.Z);
        double[][] FinalAugmentMatrix = WinogradMethod(translateCenterBack,
                WinogradMethod(notRx,
                        WinogradMethod(notRy,
                                WinogradMethod(rotateZ(angleDegrees),
                                        WinogradMethod(Ry,
                                                WinogradMethod(Rx,translateCenter))))));
        return WinogradMethod(FinalAugmentMatrix, masterMatrix);
    }
    public static Point3D VectorSubtraction(Point3D point1, Point3D point2)
    {
        return new Point3D(point1.X-point2.X,point1.Y-point2.Y,point1.Z-point2.Z );
    }

    public static Point3D crossProduct(Point3D point1, Point3D point2, Point3D fixedPoint)
    {
        Point3D a = new Point3D(point1.X - fixedPoint.X, point1.Y - fixedPoint.Y, point1.Z - fixedPoint.Z);
        Point3D b = new Point3D(point2.X - fixedPoint.X, point2.Y - fixedPoint.Y, point2.Z - fixedPoint.Z);

        return new Point3D((a.Y*b.Z) - (a.Z*b.Y),(a.Z*b.X) - (a.X*b.Z),(a.X*b.Y) - (a.Y*b.X));
    }

    public static double dotProduct(Point3D vector1, Point3D vector2, Point3D anchor1, Point3D anchor2)
    {
        vector1 = VectorSubtraction(vector1,anchor1);
        vector2 = VectorSubtraction(vector2,anchor2);

        return (vector1.X*vector2.X)+(vector1.Y*vector2.Y)+(vector1.Z*vector2.Z);
    }

    public static double angleBetween(Point3D vector1, Point3D vector2, Point3D anchor1, Point3D anchor2)
    {
        return Math.toDegrees(Math.acos(dotProduct(vector1,vector2, anchor1, anchor2)/((vectorMag(vector1.X,vector1.Y,vector1.Z)) * vectorMag(vector2.X,vector2.Y,vector2.Z))));
    }

    public static void renderLine(double[][] vertices, int[][] faces, int[][]surfaceNormals, LineBase.RGBColor color, int[][][]fb)
    {
        LineBase lb = new Lines();
        for(int i = 0; i < faces[0].length;i++)
        {
            Point3D vertex1 = new Point3D(vertices[0][surfaceNormals[0][i]],vertices[1][surfaceNormals[0][i]],vertices[2][surfaceNormals[0][i]]);
            Point3D anchor1 = new Point3D(vertices[0][surfaceNormals[1][i]], vertices[1][surfaceNormals[1][i]], vertices[2][surfaceNormals[1][i]]);
            Point3D vertex2 = new Point3D(0,0,-1);
            Point3D anchor2 = new Point3D(0,0,0);
            if(angleBetween(vertex1,vertex2,anchor1,anchor2) < 90) {
                for (int j = 0; j < faces.length - 1; j++) {
                    lb.BresenhamFormRGB((int) vertices[0][faces[j][i]], (int) vertices[1][faces[j][i]], (int) vertices[0][faces[j + 1][i]], (int) vertices[1][faces[j + 1][i]], fb, color, color);
                }
            }
        }
    }
    public static void renderLine(double[][] vertices, int[][] faces, int[][]surfaceNormals, LineBase.RGBColor[] color, int[][][]fb)
    {
        LineBase lb = new Lines();


        for(int i = 0; i < faces[0].length;i++)
        {
            Point3D vertex1 = new Point3D(vertices[0][surfaceNormals[0][i]],vertices[1][surfaceNormals[0][i]],vertices[2][surfaceNormals[0][i]]);
            Point3D anchor1 = new Point3D(vertices[0][surfaceNormals[1][i]], vertices[1][surfaceNormals[1][i]], vertices[2][surfaceNormals[1][i]]);
            Point3D vertex2 = new Point3D(0,0,-1);
            Point3D anchor2 = new Point3D(0,0,0);
            if(angleBetween(vertex1,vertex2,anchor1,anchor2) < 90) {
                for (int j = 0; j < faces.length - 1; j++) {

                    lb.BresenhamFormRGB((int) vertices[0][faces[j][i]], (int) vertices[1][faces[j][i]], (int) vertices[0][faces[j + 1][i]], (int) vertices[1][faces[j + 1][i]], fb, color[faces[j][i]], color[faces[j + 1][i]]);
                }
            }
        }
    }

    public static void main(String [] args) throws IOException {

        // Need to change how I input values into my master Matrix
        File file = new File("C:\\Users\\andre\\IdeaProjects\\Rotation-HiddenSurfaceRemoval\\src\\Cube\\cube.txt");
        Scanner scan = new Scanner(file);
        String line = scan.nextLine();
        int numVertices = Integer.parseInt(line.split(" ")[1]);
        double[][] vertices = new double[4][numVertices];
        LineBase.RGBColor[] colorz = new LineBase.RGBColor[numVertices];
        for(int i = 0; i < numVertices; i++)
        {
            String[] in = scan.nextLine().split(", ");
            vertices[0][i] = Double.parseDouble(in[0]);
            vertices[1][i] = Double.parseDouble(in[1]);
            vertices[2][i] = Double.parseDouble(in[2]);
            vertices[3][i] = 1.0;
            colorz[i] = new LineBase.RGBColor(Integer.parseInt(in[3]),Integer.parseInt(in[4]),Integer.parseInt(in[5]));
        }
        line = scan.nextLine();
        int numEdges = Integer.parseInt(line.split(" ")[1]);
        //Integer[][] edges = new Integer[2][numEdges];
        for(int i = 0; i < numEdges; i++)
        {
            String[] in = scan.nextLine().split(", ");
            //edges[0][i] = Integer.parseInt(in[0]);
            //edges[1][i] = Integer.parseInt(in[1]);

        }
        line = scan.nextLine();;
        int numFaces = Integer.parseInt(line.split(" ")[1]);
        int[][] faces = new int[5][numFaces];
        for(int i = 0; i <numFaces; i++ )
        {
            line = scan.nextLine();
            String[] corners = line.split(", ");
            faces[0][i] = Integer.parseInt(corners[1]);
            faces[1][i] = Integer.parseInt(corners[2]);
            faces[2][i] = Integer.parseInt(corners[3]);
            faces[3][i] = Integer.parseInt(corners[4]);
            faces[4][i] = Integer.parseInt(corners[5]);
        }

        double[][] surfaceCentroids = new double[4][numFaces];
        double[][] surfaceNormals = new double[4][numFaces];
        for(int i = 0; i < numFaces;i++)
        {
            Point3D surfaceNormal = crossProduct(new Point3D(vertices[0][faces[0][i]],vertices[1][faces[0][i]], vertices[2][faces[0][i]]),
                    new Point3D(vertices[0][faces[2][i]], vertices[1][faces[2][i]], vertices[2][faces[2][i]]),
                    new Point3D(vertices[0][faces[1][i]],vertices[1][faces[1][i]], vertices[2][faces[1][i]]));

            surfaceNormals[0][i] = surfaceNormal.X;
            surfaceNormals[1][i] = surfaceNormal.Y;
            surfaceNormals[2][i] = surfaceNormal.Z;
            surfaceNormals[3][i] = 1;

            double sumx = 0;
            double sumy = 0;
            double sumz = 0;
            for(int j = 0; j < faces.length -1; j++)
            {
                sumx += vertices[0][faces[j][i]];
                sumy += vertices[1][faces[j][i]];
                sumz += vertices[2][faces[j][i]];
            }

            surfaceCentroids[0][i] = sumx/(faces.length-1);
            surfaceCentroids[1][i] = sumy/(faces.length-1);
            surfaceCentroids[2][i] = sumz/(faces.length-1);
            surfaceCentroids[3][i] = 1;
        }
        double[][] everythingCube = new double[4][numVertices + (numFaces*2) + 1];
        for(int i = 0; i < vertices[0].length; i++)
        {
            everythingCube[0][i] = vertices[0][i];
            everythingCube[1][i] = vertices[1][i];
            everythingCube[2][i] = vertices[2][i];
            everythingCube[3][i] = vertices[3][i];
        }
        for (int i = vertices[0].length; i < surfaceNormals[0].length+vertices[0].length; i++)
        {
            everythingCube[0][i] = surfaceNormals[0][i-vertices[0].length];
            everythingCube[1][i] = surfaceNormals[1][i-vertices[0].length];
            everythingCube[2][i] = surfaceNormals[2][i-vertices[0].length];
            everythingCube[3][i] = surfaceNormals[3][i-vertices[0].length];
        }
        for(int i = surfaceNormals[0].length+vertices[0].length; i < everythingCube[0].length-1; i++)
        {
            everythingCube[0][i] = surfaceCentroids[0][i-(surfaceNormals[0].length+vertices[0].length)];
            everythingCube[1][i] = surfaceCentroids[1][i-(surfaceNormals[0].length+vertices[0].length)];
            everythingCube[2][i] = surfaceCentroids[2][i-(surfaceNormals[0].length+vertices[0].length)];
            everythingCube[3][i] = surfaceCentroids[3][i-(surfaceNormals[0].length+vertices[0].length)];
        }


        double centerX = 0;
        double centerY = 0;
        double centerZ = 0;
        for(int i = 0; i < vertices[0].length; i++)
        {
            centerX += vertices[0][i];
            centerY += vertices[1][i];
            centerZ += vertices[2][i];
        }
        everythingCube[0][everythingCube[0].length-1] = centerX / vertices[0].length;
        everythingCube[1][everythingCube[1].length-1] = centerY / vertices[0].length;
        everythingCube[2][everythingCube[2].length-1] = centerZ / vertices[0].length;
        everythingCube[3][everythingCube[3].length-1] = 1;
        for(int i = 0; i < everythingCube.length; i++)
        {
            for(int j = 0; j < everythingCube[0].length; j++)
            {
                System.out.print(everythingCube[i][j] + ", ");
            }
            System.out.println();
        }


//        for(int i = 0; i < edges.length; i++)
//        {
//            for(int j = 0; j < edges[i].length; j++)
//            {
//                System.out.print(edges[i][j]);
//
//            }
//        }
//        for(int i = 0; i < vertices.length; i++)
//        {
//            for(int j = 0; j < vertices[i].length; j++)
//            {
//                System.out.print(vertices[i][j] + ", ");
//            }
//            System.out.println();
//        }
        int[][] normalsCentroids= new int[2][numFaces];
        for(int i = 0; i < faces[0].length; i++)
        {
            normalsCentroids[0][i] = i + vertices[0].length;
            normalsCentroids[1][i] = i + vertices[0].length + faces[0].length;
        }

        int[][][] framebuffer = new int[3][512][512];
        double frambufferlength = (double)framebuffer[0].length/2;
        System.out.println(frambufferlength);
        double[][] FinalMatrix = rotation(everythingCube,new Point3D(1,1,1), 57);
        double[][] temp = WinogradMethod(scale(100,100,100), FinalMatrix);
        double[][] result = WinogradMethod(translate((double)framebuffer[0].length/2, (double)framebuffer[0].length/2,0 ), temp);
        renderLine(result, faces, normalsCentroids, colorz, framebuffer);
        LineBase.ImageWriteRGB(framebuffer,"cube.png");
        renderLine(result,normalsCentroids, normalsCentroids, new LineBase.RGBColor(255,255,255), framebuffer);
        LineBase.ImageWriteRGB(framebuffer,"cube1.png");

    }



}

