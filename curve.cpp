#ifdef __APPLE__
#include <OpenGL/gl.h>
/* Just in case we need these later
// References:
// http://alumni.cs.ucsb.edu/~wombatty/tutorials/opengl_mac_osx.html
// # include <OpenGL/gl.h>
// # include <OpenGL/glu.h>
*/
#else
#include <GL/gl.h>
#endif

#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx(const Vector3f &lhs, const Vector3f &rhs)
    {
        const float eps = 1e-8f;
        return (lhs - rhs).absSquared() < eps;
    }

}

Vector3f prevBinorm;

Curve evalBezier(const vector<Vector3f> &P, unsigned steps)
{
    // Check
    if (P.size() < 4 || P.size() % 3 != 1)
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit(0);
    }

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    Curve B((steps + 1) * (P.size() - 3));
    cout << "Number of Control Points : " << P.size() << endl;
    for (int i = 0; i < P.size(); ++i)
    {
        cout << P[i][0] << P[i][1] << P[i][2] << endl;
    }

    //V = (1−t)^3P1 + 3(1−t)^2tP2 +3(1−t)t^2P3 + t^3P4
    //T = -3(1-t)^2P1 + 3(1-t)^2P2 - 6t(1-t)P2 - 3t^2P3 + 6t(1-t)P3 + 3t^2P4
    // N = prevBxT.normalize()
    // B = TxN.normalize()
    int Index = 0;
    for (int j = 0; j < P.size() - 3; j += 3)
    {
        Vector3f Binorm;

        for (int i = 0; i <= steps; i++)
        {
            if (prevBinorm.x() == 0 && prevBinorm.y() == 0 && prevBinorm.z() == 0)
            {
                Binorm = Vector3f(0.0f, 0.0f, 1.0f);
            }
            else
            {
                Binorm = prevBinorm;
            }
            float t = float(i) / steps;
            float x = float(pow(1 - t, 3) * P[j][0] + 3 * pow(1 - t, 2) * (t)*P[j + 1][0] + 3 * (1 - t) * pow(t, 2) * P[j + 2][0] + pow(t, 3) * P[j + 3][0]);
            float y = float(pow(1 - t, 3) * P[j][1] + 3 * pow(1 - t, 2) * (t)*P[j + 1][1] + 3 * (1 - t) * pow(t, 2) * P[j + 2][1] + pow(t, 3) * P[j + 3][1]);
            float z = float(pow(1 - t, 3) * P[j][2] + 3 * pow(1 - t, 2) * (t)*P[j + 1][2] + 3 * (1 - t) * pow(t, 2) * P[j + 2][2] + pow(t, 3) * P[j + 3][2]);
            B[Index].V = Vector3f(x, y, z);
            x = -3 * pow((1 - t), 2) * P[j][0] + 3 * pow((1 - t), 2) * P[j + 1][0] - 6 * t * (1 - t) * P[j + 1][0] - 3 * pow(t, 2) * P[j + 2][0] + 6 * t * (1 - t) * P[j + 2][0] + 3 * pow(t, 2) * P[j + 3][0];
            y = -3 * pow((1 - t), 2) * P[j][1] + 3 * pow((1 - t), 2) * P[j + 1][1] - 6 * t * (1 - t) * P[j + 1][1] - 3 * pow(t, 2) * P[j + 2][1] + 6 * t * (1 - t) * P[j + 2][1] + 3 * pow(t, 2) * P[j + 3][1];
            z = -3 * pow((1 - t), 2) * P[j][2] + 3 * pow((1 - t), 2) * P[j + 1][2] - 6 * t * (1 - t) * P[j + 1][2] - 3 * pow(t, 2) * P[j + 2][2] + 6 * t * (1 - t) * P[j + 2][2] + 3 * pow(t, 2) * P[j + 3][2];
            B[Index].T = Vector3f(x, y, z).normalized();
            x = 6 * (1 - t) * P[j][0] - 6 * (1 - t) * P[j + 1][0] - 6 * P[j + 1][0] + 12 * t * P[j + 1][0] - 6 * t * P[2][0] + 6 * P[j + 2][0] - 12 * t * P[j + 2][0] + 6 * t * P[j + 3][0];
            y = 6 * (1 - t) * P[j][1] - 6 * (1 - t) * P[j + 1][1] - 6 * P[j + 1][1] + 12 * t * P[j + 1][1] - 6 * t * P[2][1] + 6 * P[j + 2][1] - 12 * t * P[j + 2][1] + 6 * t * P[j + 3][1];
            z = 6 * (1 - t) * P[j][2] - 6 * (1 - t) * P[j + 1][2] - 6 * P[j + 1][2] + 12 * t * P[j + 1][2] - 6 * t * P[2][2] + 6 * P[j + 2][2] - 12 * t * P[j + 2][2] + 6 * t * P[j + 3][2];

            B[Index].N = Vector3f::cross(Binorm, B[i].T).normalized();
            B[Index].B = Vector3f::cross(B[i].T, B[i].N).normalized();
            prevBinorm = B[Index].B;
            Index++;
        }
    }

    cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
    for (unsigned i = 0; i < P.size(); ++i)
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning empty curve." << endl;

    // Right now this will just return this empty curve.
    return B;
}

Curve evalBspline(const vector<Vector3f> &P, unsigned steps)
{
    // Check
    if (P.size() < 4)
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit(0);
    }
    cout << "Number of control points : " << P.size() << endl;
    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    // Convert B spline control points into Bezier control points.
    Curve BFinal((steps + 1) * (P.size() - 3));
    Curve B(steps + 1);

    Matrix4f Bspline(1.0f, -3.0f, 3.0f, -1.0f,
                     4.0f, 0.0f, -6.0f, 3.0f,
                     1.0f, 3.0f, 3.0f, -3.0f,
                     0.0f, 0.0f, 0.0f, 1.0f);
    Bspline /= 6;
    Matrix4f BezierInverse(1.0f, 1.0f, 1.0f, 1.0f,
                           0.0f, 0.333333333f, 0.666666667f, 1.0f,
                           0.0f, 0.0f, 0.333333333f, 1.0f,
                           0.0f, 0.0f, 0.0f, 1.0f);
    Matrix4f Product;
    Product = Bspline * BezierInverse;
    int Index = 0;
    for (int i = 0; i < P.size() - 3; i++)
    {
        vector<Vector3f> NewControlPoints;
        float x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
        x1 = P[i].x() * Product(0, 0) + P[i + 1].x() * Product(1, 0) + P[i + 2].x() * Product(2, 0) + P[i + 3].x() * Product(3, 0);
        x2 = P[i].x() * Product(0, 1) + P[i + 1].x() * Product(1, 1) + P[i + 2].x() * Product(2, 1) + P[i + 3].x() * Product(3, 1);
        x3 = P[i].x() * Product(0, 2) + P[i + 1].x() * Product(1, 2) + P[i + 2].x() * Product(2, 2) + P[i + 3].x() * Product(3, 2);
        x4 = P[i].x() * Product(0, 3) + P[i + 1].x() * Product(1, 3) + P[i + 2].x() * Product(2, 3) + P[i + 3].x() * Product(3, 3);

        y1 = P[i].y() * Product(0, 0) + P[i + 1].y() * Product(1, 0) + P[i + 2].y() * Product(2, 0) + P[i + 3].y() * Product(3, 0);
        y2 = P[i].y() * Product(0, 1) + P[i + 1].y() * Product(1, 1) + P[i + 2].y() * Product(2, 1) + P[i + 3].y() * Product(3, 1);
        y3 = P[i].y() * Product(0, 2) + P[i + 1].y() * Product(1, 2) + P[i + 2].y() * Product(2, 2) + P[i + 3].y() * Product(3, 2);
        y4 = P[i].y() * Product(0, 3) + P[i + 1].y() * Product(1, 3) + P[i + 2].y() * Product(2, 3) + P[i + 3].y() * Product(3, 3);

        z1 = P[i].z() * Product(0, 0) + P[i + 1].z() * Product(1, 0) + P[i + 2].z() * Product(2, 0) + P[i + 3].z() * Product(3, 0);
        z2 = P[i].z() * Product(0, 1) + P[i + 1].z() * Product(1, 1) + P[i + 2].z() * Product(2, 1) + P[i + 3].z() * Product(3, 1);
        z3 = P[i].z() * Product(0, 2) + P[i + 1].z() * Product(1, 2) + P[i + 2].z() * Product(2, 2) + P[i + 3].z() * Product(3, 2);
        z4 = P[i].z() * Product(0, 3) + P[i + 1].z() * Product(1, 3) + P[i + 2].z() * Product(2, 3) + P[i + 3].z() * Product(3, 3);

        NewControlPoints.push_back(Vector3f(x1, y1, z1));
        NewControlPoints.push_back(Vector3f(x2, y2, z2));
        NewControlPoints.push_back(Vector3f(x3, y3, z3));
        NewControlPoints.push_back(Vector3f(x4, y4, z4));

        B = evalBezier(NewControlPoints, steps);
        for (int j = 0; j <= steps; j++)
        {
            BFinal[Index] = B[j];
            Index++;
        }
    }

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
    for (unsigned i = 0; i < P.size(); ++i)
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning empty curve." << endl;

    // Return an empty curve right now.
    return BFinal;
}

Curve evalcatMullRomspline(const std::vector<Vector3f> &P, unsigned steps)
{
    // Check
    if (P.size() < 4 || P.size() % 3 != 1)
    {
        cerr << "evalcatMullRomspline must be called with 3n+1 control points." << endl;
        exit(0);
    }

    // Convert CatMullRomspline control points into Bezier control points.
    Curve BFinal((steps + 1) * (P.size() - 3));
    Curve B(steps + 1);

    Matrix4f CatMullRomspline(0.0f, -1.0f, 2.0f, -1.0f,
                              2.0f, 0.0f, -5.0f, 3.0f,
                              0.0f, 1.0f, 4.0f, -3.0f,
                              0.0f, 0.0f, -1.0f, 1.0f);
    CatMullRomspline /= 2;

    Matrix4f BezierInverse(1.0f, 1.0f, 1.0f, 1.0f,
                           0.0f, 0.333333333f, 0.666666667f, 1.0f,
                           0.0f, 0.0f, 0.333333333f, 1.0f,
                           0.0f, 0.0f, 0.0f, 1.0f);
    Matrix4f Product;
    Product = CatMullRomspline * BezierInverse;
    int Index = 0;
    for (int i = 0; i < P.size() - 3; i++)
    {
        vector<Vector3f> NewControlPoints;
        float x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
        x1 = P[i].x() * Product(0, 0) + P[i + 1].x() * Product(1, 0) + P[i + 2].x() * Product(2, 0) + P[i + 3].x() * Product(3, 0);
        x2 = P[i].x() * Product(0, 1) + P[i + 1].x() * Product(1, 1) + P[i + 2].x() * Product(2, 1) + P[i + 3].x() * Product(3, 1);
        x3 = P[i].x() * Product(0, 2) + P[i + 1].x() * Product(1, 2) + P[i + 2].x() * Product(2, 2) + P[i + 3].x() * Product(3, 2);
        x4 = P[i].x() * Product(0, 3) + P[i + 1].x() * Product(1, 3) + P[i + 2].x() * Product(2, 3) + P[i + 3].x() * Product(3, 3);

        y1 = P[i].y() * Product(0, 0) + P[i + 1].y() * Product(1, 0) + P[i + 2].y() * Product(2, 0) + P[i + 3].y() * Product(3, 0);
        y2 = P[i].y() * Product(0, 1) + P[i + 1].y() * Product(1, 1) + P[i + 2].y() * Product(2, 1) + P[i + 3].y() * Product(3, 1);
        y3 = P[i].y() * Product(0, 2) + P[i + 1].y() * Product(1, 2) + P[i + 2].y() * Product(2, 2) + P[i + 3].y() * Product(3, 2);
        y4 = P[i].y() * Product(0, 3) + P[i + 1].y() * Product(1, 3) + P[i + 2].y() * Product(2, 3) + P[i + 3].y() * Product(3, 3);

        z1 = P[i].z() * Product(0, 0) + P[i + 1].z() * Product(1, 0) + P[i + 2].z() * Product(2, 0) + P[i + 3].z() * Product(3, 0);
        z2 = P[i].z() * Product(0, 1) + P[i + 1].z() * Product(1, 1) + P[i + 2].z() * Product(2, 1) + P[i + 3].z() * Product(3, 1);
        z3 = P[i].z() * Product(0, 2) + P[i + 1].z() * Product(1, 2) + P[i + 2].z() * Product(2, 2) + P[i + 3].z() * Product(3, 2);
        z4 = P[i].z() * Product(0, 3) + P[i + 1].z() * Product(1, 3) + P[i + 2].z() * Product(2, 3) + P[i + 3].z() * Product(3, 3);

        NewControlPoints.push_back(Vector3f(x1, y1, z1));
        NewControlPoints.push_back(Vector3f(x2, y2, z2));
        NewControlPoints.push_back(Vector3f(x3, y3, z3));
        NewControlPoints.push_back(Vector3f(x4, y4, z4));

        B = evalBezier(NewControlPoints, steps);
        for (int j = 0; j <= steps; j++)
        {
            BFinal[Index] = B[j];
            Index++;
        }
    }

    cerr << "\t>>> evalcatMullRomspline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
    for (unsigned i = 0; i < P.size(); ++i)
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning empty curve." << endl;

    // Return an empty curve right now.
    return BFinal;
}

Curve evalCircle(float radius, unsigned steps)
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).

    // Preallocate a curve with steps+1 CurvePoints
    Curve R(steps + 1);

    // Fill it in counterclockwise
    for (unsigned i = 0; i <= steps; ++i)
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float(i) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f(cos(t), sin(t), 0);

        // Tangent vector is first derivative
        R[i].T = Vector3f(-sin(t), cos(t), 0);

        // Normal vector is second derivative
        R[i].N = Vector3f(-cos(t), -sin(t), 0);

        // Finally, binormal is facing up.
        R[i].B = Vector3f(0, 0, 1);
    }

    return R;
}

void drawCurve(const Curve &curve, float framesize)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    // Setup for line drawing
    glDisable(GL_LIGHTING);
    glColor4f(1, 1, 1, 1);
    glLineWidth(1);

    // Draw curve
    glBegin(GL_LINE_STRIP);
    for (unsigned i = 0; i < curve.size(); ++i)
    {
        glVertex(curve[i].V);
    }
    glEnd();

    glLineWidth(1);

    // Draw coordinate frames if framesize nonzero
    if (framesize != 0.0f)
    {
        Matrix4f M;

        for (unsigned i = 0; i < curve.size(); ++i)
        {
            M.setCol(0, Vector4f(curve[i].N, 0));
            M.setCol(1, Vector4f(curve[i].B, 0));
            M.setCol(2, Vector4f(curve[i].T, 0));
            M.setCol(3, Vector4f(curve[i].V, 1));

            glPushMatrix();
            glMultMatrixf(M);
            glScaled(framesize, framesize, framesize);
            glBegin(GL_LINES);
            glColor3f(1, 0, 0);
            glVertex3d(0, 0, 0);
            glVertex3d(1, 0, 0);
            glColor3f(0, 1, 0);
            glVertex3d(0, 0, 0);
            glVertex3d(0, 1, 0);
            glColor3f(0, 0, 1);
            glVertex3d(0, 0, 0);
            glVertex3d(0, 0, 1);
            glEnd();
            glPopMatrix();
        }
    }

    // Pop state
    glPopAttrib();
}
