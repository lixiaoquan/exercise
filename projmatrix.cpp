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

#define M_PI 3.14159265

void setProjectionMatrix(
    Mat44f &M,
    float &angleOfView,
    float &far,
    float &near
    )
{
    float scale = 1 / tan(angleOfView * 0.5 * M_PI / 180);

    M[0][0] = M[1][1] = scale;

    M[2][2] = - far / (far -near);
    M[3][2] = - far*near / (far -near);

    M[2][3] = -1;
    M[3][3] = 0;
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

    setProjectionMatrix(Mproj, angleOfView, far, near);

#if 0
    float canvasWidth = 1;
    float canvasHeight = 1;
#else
    float canvasWidth = 2;
    float canvasHeight = 2;
#endif

    uint32_t imageWidth = 512;
    uint32_t imageHeight = 512;

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
