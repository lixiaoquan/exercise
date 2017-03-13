/* Follow code in Camera rays. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <iomanip> 
#include <limits> 
#include <vector> 
#include <random> 
#include "stdint.h"
#include "string.h"

#include "geometry.h"
#include "vertexdata.h"

using namespace std;

float kEpsilon = 1e-8;

#define CULLING 1
#ifndef M_PI
#define M_PI 3.14159265
#endif

uint32_t width = 640;
uint32_t height = 480;
float angleOfView  = 51.52;

float clamp(const float &a, const float &b, const float &c)
{
    return max(a, min(b, c));
}


bool rayTriangleIntersect(
    Vec3f &orig,
    Vec3f &dir,
    Vec3f &v0,
    Vec3f &v1,
    Vec3f &v2,
    float &t,
    float &u,
    float &v
    )
{
    Vec3f v0v1 = v1 - v0;
    Vec3f v0v2 = v2 - v0;

    Vec3f pvec = dir.crossProduct(v0v2);
    float det = v0v1.dot(pvec);

#ifdef CULLING
    if (det < kEpsilon)
    {
        return false;
    }
#else
    if (fabs(det) < kEpsilon)
    {
        return false;
    }
#endif

  //  cerr << "after culling" << endl;

    float invDet = 1/ det;

    Vec3f tvec = orig - v0;

    u = tvec.dot(pvec) * invDet;

    if (u < 0 || u > 1)
    {
        return false;
    }

    Vec3f qvec = tvec.crossProduct(v0v1);

    v = dir.dot(qvec) * invDet;

    if (v < 0 || u + v > 1) return false;

    t = v0v2.dot(qvec) * invDet;

   // cerr << true << endl;
    return true;
}


const float kInfinity = std::numeric_limits<float>::max();

void render (
    )
{
    Vec3f v0(-1, -1, -5);
    Vec3f v1(1, -1, -5);
    Vec3f v2(0, 1, -5);
#ifdef CULLING
    Vec3f vertices[6] = {{0, 1, -5}, {1, -1, -5}, {-1, -1, -5}, {-1, -1, -5}, {2, -1, -6}, {0, 2, -6}};
#else
    Vec3f vertices[6] = {{-1, -1, -5}, {1, -1, -5}, {0, 1, -5}, {-1, -1, -5}, {2, -1, -6}, {0, 2, -6}};
#endif

    Vec3f cols[3] = {{0.6, 0.4, 0.1}, {0.1, 0.5, 0.3}, {0.1, 0.3, 0.7}};

    Mat44f cameraToWorld = Mat44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

    float imageAspectRatio = width / (float)height;
    float scale = tan(angleOfView / 2 * M_PI / 180);


    Vec3f *buffer = new Vec3f[width * height];
    Vec3f *pix = buffer;

    memset(buffer, 0, width * height);

    Vec3f orig(0);

    //cameraToWorld.multiple(Vec3f(0), orig);

    cerr << "orig: " << orig << endl;

    for (uint32_t j = 0; j < height; j++)
    {
        for (uint32_t i = 0; i < width; i++)
        {
            /* NDC */

            float x = (2 * (i + 0.5) / (float)width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)height) * scale;

            Vec3f dir(x, y, -1);

            /* NOTE this is a dir transform. */
//            cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();

            float t, u, v;
            float tNear;

            tNear = kInfinity;
            for (int k = 0; k < 2; k++)
            {
                v0 = vertices[k * 3];
                v1 = vertices[k * 3 + 1];
                v2 = vertices[k * 3 + 2];

                if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v))
                {
                    if (t < tNear)
                    {
                        tNear = t;
                        *pix = cols[0] * u + cols[1] * v + cols[2] * (1 - u - v);
                    }
                }
            }

            pix++;

        }
    }


    // save to file
    std::ofstream ofs;
    ofs.open("./raytri.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";


     for (uint32_t i = 0; i < height * width; ++i) {
     char r = (char)(255 * clamp(0, 1, buffer[i].x));
     char g = (char)(255 * clamp(0, 1, buffer[i].y));
     char b = (char)(255 * clamp(0, 1, buffer[i].z));
     ofs << r << g << b;
     } 
    ofs.write((char*)buffer, width * height);
    ofs.close();

    delete [] buffer;
}

int main()
{
    render();
}

