/* Follow code in Virtual Camera. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <iomanip> 
#include "stdint.h"

#include "geometry.h"

/* Vertices data from lesson page. */
const Vec3f verts[146] = { { -2.5703, 0.78053, -2.4e-05}, { -0.89264, 0.022582, 0.018577}, { 1.6878, -0.017131, 0.022032}, { 3.4659, 0.025667, 0.018577}, { -2.5703, 0.78969, -0.001202}, { -0.89264, 0.25121, 0.93573}, { 1.6878, 0.25121, 1.1097}, { 3.5031, 0.25293, 0.93573}, { -2.5703, 1.0558, -0.001347}, { -0.89264, 1.0558, 1.0487}, { 1.6878, 1.0558, 1.2437}, { 3.6342, 1.0527, 1.0487}, { -2.5703, 1.0558, 0}, { -0.89264, 1.0558, 0}, { 1.6878, 1.0558, 0}, { 3.6342, 1.0527, 0}, { -2.5703, 1.0558, 0.001347}, { -0.89264, 1.0558, -1.0487}, { 1.6878, 1.0558, -1.2437}, { 3.6342, 1.0527, -1.0487}, { -2.5703, 0.78969, 0.001202}, { -0.89264, 0.25121, -0.93573}, { 1.6878, 0.25121, -1.1097}, { 3.5031, 0.25293, -0.93573}, { 3.5031, 0.25293, 0}, { -2.5703, 0.78969, 0}, { 1.1091, 1.2179, 0}, { 1.145, 6.617, 0}, { 4.0878, 1.2383, 0}, { -2.5693, 1.1771, -0.081683}, { 0.98353, 6.4948, -0.081683}, { -0.72112, 1.1364, -0.081683}, { 0.9297, 6.454, 0}, { -0.7929, 1.279, 0}, { 0.91176, 1.2994, 0} }; 

const int numTris = 51;

const uint32_t tris[numTris * 3] = {
4, 0, 5, 0, 1, 5, 1, 2, 5, 5, 2, 6, 3, 7, 2,
2, 7, 6, 5, 9, 4, 4, 9, 8, 5, 6, 9, 9, 6, 10,
7, 11, 6, 6, 11, 10, 9, 13, 8, 8, 13, 12, 10, 14, 9,
9, 14, 13, 10, 11, 14, 14, 11, 15, 17, 16, 13, 12, 13, 16,
13, 14, 17, 17, 14, 18, 15, 19, 14, 14, 19, 18, 16, 17, 20,
20, 17, 21, 18, 22, 17, 17, 22, 21, 18, 19, 22, 22, 19, 23,
20, 21, 0, 21, 1, 0, 22, 2, 21, 21, 2, 1, 22, 23, 2,
2, 23, 3, 3, 23, 24, 3, 24, 7, 24, 23, 15, 15, 23, 19,
24, 15, 7, 7, 15, 11, 0, 25, 20, 0, 4, 25, 20, 25, 16,
16, 25, 12, 25, 4, 12, 12, 4, 8, 26, 27, 28, 29, 30, 31,
32, 34, 33
};

bool computerPixelCoordinate(
    Vec3f &WorldCord,
    Vec2i &PixelCorrdinate,
    Mat44f &worldToCamera,
    float &l,
    float &b,
    float &r,
    float &t,
    float &zNear,
    uint32_t &ImageWidth,
    uint32_t &ImageHeight
    )
{
    /* Convert to camera space. */
    Vec3f cameraCord;
    Vec2f screenCord;
    Vec2f ndcCord;

    std::cerr << "World Coord: " << WorldCord << std::endl;

    worldToCamera.multiple(WorldCord, cameraCord);

    std::cerr << "Camera Coord:" << cameraCord << std::endl;

    screenCord.x = - cameraCord.x / cameraCord.z * zNear;
    screenCord.y = - cameraCord.y / cameraCord.z * zNear;

    std::cerr << "Screen Coord:" << screenCord << std::endl;

    /* NDC. */
    ndcCord.x = (screenCord.x + r) / (2 * r);
    ndcCord.y = (screenCord.y + t) / (2 * t);

    std::cerr << "NDC Coord:" << ndcCord << std::endl;

    /* Raster. */
    PixelCorrdinate.x = ndcCord.x * ImageWidth;
    PixelCorrdinate.y = (1 - ndcCord.y) * ImageHeight;

    if (screenCord.x >= l && screenCord.x <= r && screenCord.y >= b && screenCord.y <= t )
        return true;

    return false;
}

int main()
{
    Mat44f worldToCamera;

    float filmApertureWidth = 0.825;
    float filmApertureHeight = 0.446;

    float focalLength = 35; // mm

    float zNear  = 0.1;
    float zFar   = 1000;

    float canvasWidth = 2;
    float canvasHeight = 2;

    float right;
    float left;
    float top;
    float bottom;

    uint32_t imageWidth = 512;
    uint32_t imageHeight = 512;

    top = (filmApertureHeight * 25.4 / 2 / focalLength) * zNear;

    bottom = -top;

    right = (filmApertureWidth * 25.4 / 2 / focalLength) * zNear;

    left = - right;

    std::cerr << left << " " << bottom << " " << right <<  " " << top << std::endl;

    Mat44f cameraToWorld(-0.95424, 0, 0.299041, 0, 0.0861242, 0.95763, 0.274823, 0, -0.28637, 0.288002, -0.913809, 0, -3.734612, 7.610426, -14.152769, 1);

    worldToCamera = cameraToWorld.inverse();

    std::ofstream ofs; 
    ofs.open("./pinhole.svg"); 
    ofs << "<svg version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns=\"http://www.w3.org/2000/svg\" height=\"512\" width=\"512\">" << std::endl; 
    
    for (int i = 0; i < numTris; i++)
    {
        Vec3f v0World = verts[tris[i * 3]];
        Vec3f v1World = verts[tris[i * 3 + 1]];
        Vec3f v2World = verts[tris[i * 3 + 2]];
        Vec2i v0Raster;
        Vec2i v1Raster;
        Vec2i v2Raster;
        bool v0Visiable;
        bool v1Visiable;
        bool v2Visiable;

        v0Visiable = computerPixelCoordinate(v0World, v0Raster, worldToCamera, left, bottom, right, top, zNear, imageWidth, imageHeight);
        v1Visiable = computerPixelCoordinate(v1World, v1Raster, worldToCamera, left, bottom, right, top, zNear, imageWidth, imageHeight);
        v2Visiable = computerPixelCoordinate(v2World, v2Raster, worldToCamera, left, bottom, right, top, zNear, imageWidth, imageHeight);
        std::cerr << v0Raster << ", " << v1Raster << ", " << v2Raster << std::endl;

        if (v0Visiable && v1Visiable)
        {
            ofs << "<line x1=\"" << v0Raster.x << "\" y1=\"" << v0Raster.y << "\" x2=\"" << v1Raster.x << "\" y2=\"" << v1Raster.y << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
        }
        else
        {
            ofs << "<line x1=\"" << v0Raster.x << "\" y1=\"" << v0Raster.y << "\" x2=\"" << v1Raster.x << "\" y2=\"" << v1Raster.y << "\" style=\"stroke:rgb(0,0,255);stroke-width:1\" />\n";
        }

        if (v1Visiable && v2Visiable)
        {
            ofs << "<line x1=\"" << v1Raster.x << "\" y1=\"" << v1Raster.y << "\" x2=\"" << v2Raster.x << "\" y2=\"" << v2Raster.y << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
        }
        else
        {
            ofs << "<line x1=\"" << v1Raster.x << "\" y1=\"" << v1Raster.y << "\" x2=\"" << v2Raster.x << "\" y2=\"" << v2Raster.y << "\" style=\"stroke:rgb(0,0,255);stroke-width:1\" />\n";
        }

        if (v0Visiable && v1Visiable)
        {
            ofs << "<line x1=\"" << v2Raster.x << "\" y1=\"" << v2Raster.y << "\" x2=\"" << v0Raster.x << "\" y2=\"" << v0Raster.y << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
        }
        else
        {
            ofs << "<line x1=\"" << v2Raster.x << "\" y1=\"" << v2Raster.y << "\" x2=\"" << v0Raster.x << "\" y2=\"" << v0Raster.y << "\" style=\"stroke:rgb(0,0,255);stroke-width:1\" />\n";
        }

    }

    ofs << "</svg>\n";
    ofs.close();


    return 0;
}
