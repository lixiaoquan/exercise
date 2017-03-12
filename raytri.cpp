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

#ifndef M_PI
#define M_PI 3.14159265
#endif

uint32_t width = 640;
uint32_t height = 480;
float angleOfView  = 51.52;

void render (
    )
{
    Mat44f cameraToWorld = Mat44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

    float imageAspectRatio = width / (float)height;
    float scale = tan(angleOfView / 2 * M_PI / 180);


    Vec3f *buffer = new Vec3f[width * height];

    memset(buffer, 0, width * height);

    Vec3f orig;

    cameraToWorld.multiple(Vec3f(0), orig);

    cerr << "orig: " << orig << endl;

    for (uint32_t j = 0; j < height; j++)
    {
        for (uint32_t i = 0; i < width; i++)
        {
            /* NDC */

            float x = (2 * (i + 0.5) / (float)width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)height) * scale;

            Vec3f dir;

            /* NOTE this is a dir transform. */
            cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();

        }
    }


    // save to file
    std::ofstream ofs;
    ofs.open("./simpleshapes.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";


#if 0
     for (uint32_t i = 0; i < height * width; ++i) {
     char r = (char)(255 * clamp(0, 1, buffer[i].x));
     char g = (char)(255 * clamp(0, 1, buffer[i].y));
     char b = (char)(255 * clamp(0, 1, buffer[i].z));
     ofs << r << g << b;
     } 
#endif
    ofs.write((char*)buffer, width * height);
    ofs.close();

    delete [] buffer;
}

int main()
{
    render();
}

