#include "surf.h"
#include "extra.h"
using namespace std;

namespace
{

    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i = 0; i < profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0)
                return false;

        return true;
    }
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.

    // Matrix_M * N * Matrix_T * V = 0
    // Matrix_T -> Rotation Matrix around Y axis
    // Matrix_M -> -1 * inverse(transpose(getSubMatrix(Matrix_T)))

    for (unsigned i = 0; i < profile.size(); i++)
    {
        for (unsigned j = 0; j <= steps; j++)
        {
            // Calculate t
            float t = 2.0f * 3.14159265358979323846264338327950288 * float(j) / steps; // 2*PI*r Circumference equally divided

            // Vertex Rotation Matrix about Y axis
            Matrix4f Matrix_T = Matrix4f::rotateY(t);

            // Normal Rotation Matrix
            Matrix3f Matrix_M = Matrix_T.getSubmatrix3x3(0, 0).transposed().inverse();

            // Vertex and Normal Calculation
            Vector4f VertexCalculation = Matrix_T * Vector4f(profile[i].V.x(), profile[i].V.y(), profile[i].V.z(), 1.0f);
            Vector3f finalVertex = Vector3f(VertexCalculation.x(), VertexCalculation.y(), VertexCalculation.z());
            Vector3f finalNormal = Matrix_M * profile[i].N; // Multiplying with -1 to force it to be left of Y axis

            // Push Calculated Vertex and Normal into surface
            surface.VV.push_back(finalVertex);
            surface.VN.push_back(-1 * finalNormal);
        }
    }

    // Generate Triangles according to Number of Surface Vertices
    for (unsigned i = 0; i < surface.VV.size() - (steps + 1); i++)
    {
        Tup3u T1, T2; // 3 Points per tuple to save vertices or triangle

        if ((i + 1) % (steps + 1))
        {
            T1 = Tup3u(i + 1, i, i + (steps + 1));
            T2 = Tup3u(i + 1, (i + 1) + steps, (i + 1) + (steps + 1));
        }

        surface.VF.push_back(T1);
        surface.VF.push_back(T2);
    }

    cerr << "\t>>> makeSurfRev called (but not implemented).\n\t>>> Returning empty surface." << endl;

    return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep)
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.
    for (unsigned i = 0; i < profile.size(); i++)
    {
        for (unsigned j = 0; j < sweep.size(); j++)
        {
            // Create Sweep Matrix
            Matrix4f Matrix_Sweep(sweep[j].N[0], sweep[j].B[0], sweep[j].T[0], sweep[j].V[0],
                                  sweep[j].N[1], sweep[j].B[1], sweep[j].T[1], sweep[j].V[1],
                                  sweep[j].N[2], sweep[j].B[2], sweep[j].T[2], sweep[j].V[2],
                                  0.0f, 0.0f, 0.0f, 1.0f);

            // Normal Rotation Matrix
            Matrix3f Matrix_M = Matrix_Sweep.getSubmatrix3x3(0, 0).transposed().inverse();

            // Vertex and Normal Calculation
            Vector4f VertexCalculation = Matrix_Sweep * Vector4f(profile[i].V.x(), profile[i].V.y(), profile[i].V.z(), 1.0f);
            Vector3f finalVertex = Vector3f(VertexCalculation.x(), VertexCalculation.y(), VertexCalculation.z());
            Vector3f finalNormal = Matrix_M * profile[i].N; // Multiplying with -1 to force it to be left of Y axis

            // Push Calculated Vertex and Normal into surface
            surface.VV.push_back(finalVertex);
            surface.VN.push_back(-1 * finalNormal);
        }
    }

    //Calculate faces once all the vertices are added
    for (unsigned i = 0; i < surface.VV.size() - sweep.size(); i++)
    {
        Tup3u T1, T2; // 3 Points per tuple to save vertices or triangle

        if ((i + 1) % sweep.size())
        {
            T1 = Tup3u(i + 1, i, i + sweep.size());
            T2 = Tup3u(i + 1, i + sweep.size(), (i + 1) + sweep.size());
        }

        surface.VF.push_back(T1);
        surface.VF.push_back(T2);
    }

    cerr << "\t>>> makeGenCyl called (but not implemented).\n\t>>> Returning empty surface." << endl;

    return surface;
}

void drawSurface(const Surface &surface, bool shaded)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (shaded)
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    else
    {
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glColor4f(0.4f, 0.4f, 0.4f, 1.f);
        glLineWidth(1);
    }

    glBegin(GL_TRIANGLES);
    for (unsigned i = 0; i < surface.VF.size(); i++)
    {
        glNormal(surface.VN[surface.VF[i][0]]);
        glVertex(surface.VV[surface.VF[i][0]]);
        glNormal(surface.VN[surface.VF[i][1]]);
        glVertex(surface.VV[surface.VF[i][1]]);
        glNormal(surface.VN[surface.VF[i][2]]);
        glVertex(surface.VV[surface.VF[i][2]]);
    }
    glEnd();

    glPopAttrib();
}

void drawNormals(const Surface &surface, float len)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    glColor4f(0, 1, 1, 1);
    glLineWidth(1);

    glBegin(GL_LINES);
    for (unsigned i = 0; i < surface.VV.size(); i++)
    {
        glVertex(surface.VV[i]);
        glVertex(surface.VV[i] + surface.VN[i] * len);
    }
    glEnd();

    glPopAttrib();
}

void outputObjFile(ostream &out, const Surface &surface)
{

    for (unsigned i = 0; i < surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (unsigned i = 0; i < surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;

    for (unsigned i = 0; i < surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j = 0; j < 3; j++)
        {
            unsigned a = surface.VF[i][j] + 1;
            out << a << "/"
                << "1"
                << "/" << a << " ";
        }
        out << endl;
    }
}
