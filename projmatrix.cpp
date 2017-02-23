/* Follow code in Projective projection matrix. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <iomanip> 
#include "stdint.h"

#include "geometry.h"
#include "vertexdata.h"

int main()
{
    Mat44f worldToCamera;

    worldToCamera[3][1] = -10;
    worldToCamera[3][2] = -20;

    float angleOfView  = 90;
    float near = 0.1;
    float far = 100;

#if 0
    float canvasWidth = 1;
    float canvasHeight = 1;
#else
    float canvasWidth = 2;
    float canvasHeight = 2;
#endif

    uint32_t imageWidth = 512;
    uint32_t imageHeight = 512;

    Mat44f cameraToWorld(0.871214, 0, -0.490904, 0, -0.192902, 0.919559, -0.342346, 0, 0.451415, 0.392953, 0.801132, 0, 14.777467, 29.361945, 27.993464, 1);

    worldToCamera = cameraToWorld.inverse();

    for (int i = 0; i < numVertices; i++)
    {

    }

    return 0;
}
