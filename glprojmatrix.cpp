/* Follow code in Projective projection matrix. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <iomanip> 
#include "stdint.h"
#include "string.h"

#include "geometry.h"
#include "vertexdata.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif

void glFrustum(
    Mat44f &M,
    float &b,
    float &t,
    float &l,
    float &r,
    float &near,
    float &far
    )
{
    M[0][0] = 2 * near / (r - l);
    M[0][1] = M[0][1] = M[0][3] = 0;

    M[1][1] = 2 * near / (t - b);
    M[1][0] = M[1][2] = M[1][3] = 0;

    M[2][0] = (r + l) / (r - l);
    M[2][1] = (t + b) / (t - b);
    M[2][2] = - (far + near) / (far - near);
    M[2][3] = -1;

    M[3][0] = 0;
    M[3][1] = 0;
    M[3][2] = - 2 * far * near / (far -near);
    M[3][3] = 0;
}

void gluPerspective(
    float &angleOfView,
    float imageAspectRatio,
    float &near, float &far,
    float &b, float &t, float &l, float &r
    )
{
    t = near * tan(angleOfView * 0.5 * M_PI / 180);
    r = imageAspectRatio * t;
    l = -r;
    b = -t;
}

int main()
{
    Mat44f worldToCamera;
    Mat44f Mproj;


    worldToCamera[3][1] = -10;
    worldToCamera[3][2] = -20;

    float angleOfView  = 90;
    float near = 0.1;
    float far = 100;
    float b, t, l, r=0;

    uint32_t imageWidth = 512;
    uint32_t imageHeight = 512;

    gluPerspective(angleOfView, imageWidth / (float)imageHeight, near, far, b, t, l, r);

    glFrustum(Mproj, b, t, l, r, near, far);

    uint8_t *buffer = new uint8_t[imageWidth * imageHeight];

    memset(buffer, 0, imageWidth * imageHeight);

    std::cerr << Mproj << std::endl;
    std::cerr << worldToCamera << std::endl;

    for (int i = 0; i < numVertices; i++)
    {
        Vec3f world = vertices[i];
        Vec3f camera;
        Vec3f proj;

        worldToCamera.multiple(world, camera);
        Mproj.multiple(camera, proj);

        if (proj.x < -1 || proj.x > 1 || proj.y < -1 || proj.y > 1)
        {
            continue;
        }

        Vec2f ndc;

        /* [0, 1] */
        ndc.x = (proj.x + 1) / 2;
        ndc.y = 1 - (proj.y + 1) / 2;

        uint32_t x = std::min(imageWidth - 1, (uint32_t)(ndc.x * imageWidth));
        uint32_t y = std::min(imageHeight - 1, (uint32_t)(ndc.y * imageHeight));

        buffer[y * imageWidth + x] = 255;
    }


    // save to file
    std::ofstream ofs;
    ofs.open("./projmatrix.ppm");
    ofs << "P5\n" << imageWidth << " " << imageHeight << "\n255\n";
    ofs.write((char*)buffer, imageWidth * imageHeight);
    ofs.close();

    delete [] buffer;

    return 0;
}
