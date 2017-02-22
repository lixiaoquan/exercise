/* Follow code in Rasterization. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <iomanip> 
#include "stdint.h"

#include "geometry.h"
#include "cow.h"

const int numTris = 3156;

void computerPixelCoordinate(
    Vec3f &WorldCord,
    Vec3f &PixelCorrdinate,
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

//    std::cerr << "World Coord: " << WorldCord << std::endl;

    worldToCamera.multiple(WorldCord, cameraCord);

//    std::cerr << "Camera Coord:" << cameraCord << std::endl;

    screenCord.x = - cameraCord.x / cameraCord.z * zNear;
    screenCord.y = - cameraCord.y / cameraCord.z * zNear;

//    std::cerr << "Screen Coord:" << screenCord << std::endl;

    /* NDC. [-1, 1]*/
    ndcCord.x = 2 * screenCord.x / (r - l) - (r + l) / (r - l);
    ndcCord.y = 2 * screenCord.y / (t - b) - (t + b) / (t - b);

//    std::cerr << "NDC Coord:" << ndcCord << std::endl;

    /* Raster. */
    PixelCorrdinate.x = (ndcCord.x + 1) / 2 * ImageWidth;
    PixelCorrdinate.y = (1 - ndcCord.y) / 2 * ImageHeight;

    PixelCorrdinate.z = - cameraCord.z;
}

float min3(float &a, float &b, float &c)
{
    return std::min(a, std::min(b, c));
}

float max3(float &a, float &b, float &c)
{
    return std::max(a, std::max(b, c));
}

float edgeFunction(Vec3f &a, Vec3f &b, Vec3f c)
{
    return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
}

#define PRINT 0

int main()
{

    float filmApertureWidth = 0.980;
    float filmApertureHeight = 0.735;

    float focalLength = 20; // mm

    float zNear  = 1;
    float zFar   = 1000;
    float zSmallest = 0;
    float zLargest = 1000;

    float right;
    float left;
    float top;
    float bottom;

    float xscale = 1;
    float yscale = 1;

    uint32_t imageWidth = 640;
    uint32_t imageHeight = 480;

    float filmRatio = filmApertureWidth / filmApertureHeight;
    float deviceRatio =  imageWidth / imageHeight;

#if 0
    if (filmRatio > deviceRatio)
    {
        xscale = deviceRatio / filmRatio;
    }
    else
    {
        yscale = deviceRatio / filmRatio;
    }
#endif


    top = (filmApertureHeight * 25.4 / 2 / focalLength) * zNear;

    top *= yscale;

    bottom = -top;

    right = (filmApertureWidth * 25.4 / 2 / focalLength) * zNear;

    right *= xscale;

    left = - right;

#if PRINT
    std::cerr << left << " " << bottom << " " << right <<  " " << top << " " << xscale << " " << yscale << std::endl;
#endif

    Vec3<uint8_t> *frameBuffer = new Vec3<uint8_t>[imageHeight * imageWidth];

    float *depthBuffer = new float[imageHeight * imageWidth];

    for (int i = 0; i < imageHeight * imageWidth; i++)
    {
        /* Reset framebuffer and depth. */
        frameBuffer[i] = Vec3<uint8_t>(255);

        depthBuffer[i] = zFar;
    }

    Mat44f worldToCamera = {0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295, 0.624695, 0, -1.63871, -5.747777, -40.400412, 1};

    /* Per triangle loop. */
    int n = numTris;
    //int n = 1;
    for (int i = 0; i < n; i++)
    {
        Vec3f v0World = vertices[nvertices[i * 3]];
        Vec3f v1World = vertices[nvertices[i * 3 + 1]];
        Vec3f v2World = vertices[nvertices[i * 3 + 2]];
        Vec3f v0Raster;
        Vec3f v1Raster;
        Vec3f v2Raster;

        bool v0Visiable;
        bool v1Visiable;
        bool v2Visiable;

        computerPixelCoordinate(v0World, v0Raster, worldToCamera, left, bottom, right, top, zNear, imageWidth, imageHeight);
        computerPixelCoordinate(v1World, v1Raster, worldToCamera, left, bottom, right, top, zNear, imageWidth, imageHeight);
        computerPixelCoordinate(v2World, v2Raster, worldToCamera, left, bottom, right, top, zNear, imageWidth, imageHeight);
//        std::cerr << v0Raster << ", " << v1Raster << ", " << v2Raster << std::endl;

        v0Raster.z = 1/v0Raster.z;
        v1Raster.z = 1/v1Raster.z;
        v2Raster.z = 1/v2Raster.z;

        Vec2f st0 = st[stindices[i * 3]];
        Vec2f st1 = st[stindices[i * 3 + 1]];
        Vec2f st2 = st[stindices[i * 3 + 2]];

        /* Prepare for interploation. */
        st0 *= v0Raster.z;
        st1 *= v1Raster.z;
        st2 *= v2Raster.z;

        /* Find bounding box. */
        float xmin = min3(v0Raster.x, v1Raster.x, v2Raster.x);
        float ymin = min3(v0Raster.y, v1Raster.y, v2Raster.y);
        float xmax = max3(v0Raster.x, v1Raster.x, v2Raster.x);
        float ymax = max3(v0Raster.y, v1Raster.y, v2Raster.y);

        if (xmin < 0 || ymin < 0 || xmax > imageWidth - 1 || ymax > imageHeight - 1)
        {
            continue;
        }

//        std::cerr << "Bounding box : " << xmin << ymin << xmax << ymax << std::endl;

        uint32_t xstart = std::max(int32_t(0), (int32_t)std::floor(xmin));
        uint32_t xend   = std::min(int32_t(imageWidth - 1), (int32_t)std::floor(xmax));
        uint32_t ystart = std::max(int32_t(0), (int32_t)std::floor(ymin));
        uint32_t yend   = std::min(int32_t(imageHeight - 1), (int32_t)std::floor(ymax));

        float area = edgeFunction(v0Raster, v1Raster, v2Raster);

//        printf("%f %f %f area=%f\n", v0Raster.z, v1Raster.z, v2Raster.z, area);

        for (uint32_t y = ystart; y <= yend; y++)
        {
            for (uint32_t x = xstart; x <= xend; x++)
            {
                Vec3f pos (x + 0.5, y + 0.5, 0);
                float w0 = edgeFunction(v1Raster, v2Raster, pos);
                float w1 = edgeFunction(v2Raster, v0Raster, pos);
                float w2 = edgeFunction(v0Raster, v1Raster, pos);

                if (w0 >= 0 && w1 >=0 && w2 >=0)
                {
                    w0 /= area;
                    w1 /= area;
                    w2 /= area;

                    float z = 1 / (w0 * v0Raster.z + w1 * v1Raster.z + w2 * v2Raster.z);

#if PRINT
                    printf("%f %f %f z = %f\n", w0, w1, w2, z);
#endif

                    if (z < depthBuffer[y * imageWidth + x])
                    {
                        /* Check z buffer. */
                        depthBuffer[y * imageWidth + x] = z;

                        if (zSmallest == 0.f)
                        {
                            zSmallest = z;
                        }

                        if (zLargest == 1000.f)
                        {
                            zLargest = z;
                        }

                        zSmallest = std::min(z, zSmallest);
                        zLargest = std::max(z, zLargest);


                        /* interploated attribute. */
                        Vec2f st = (st0 * w0 + st1 * w1 + st2 * w2) * z;

#if PRINT
                        std::cerr << st0 * w0 << " " << " " <<  st1 * w1 << " " << st2 * w2 << std::endl;
                        std::cerr << (st0 * w0 + st1 * w1 + st2 * w2) << std::endl;
                        std::cerr << st0 << " " << " " <<  st1 << " " << st2 << " " << " st: " << st << std::endl;
#endif

                        /* Shading. */

                        /* Get camera space coordinate. */
                        Vec3f v0Cam, v1Cam, v2Cam;

                        worldToCamera.multiple(v0World, v0Cam);
                        worldToCamera.multiple(v1World, v1Cam);
                        worldToCamera.multiple(v2World, v2Cam);

                        float px = (v0Cam.x / -v0Cam.z) * w0 + (v1Cam.x / -v1Cam.z) * w1 + (v2Cam.x / -v2Cam.z) * w2;
                        float py = (v0Cam.y / -v0Cam.z) * w0 + (v1Cam.y / -v1Cam.z) * w1 + (v2Cam.y / -v2Cam.z) * w2;

                        Vec3f pt(px * z, py * z, -z);

#if PRINT
                        std::cerr << "pt: " << pt << std::endl;
#endif

                        Vec3f n = (v1Cam - v0Cam).crossProduct(v2Cam - v0Cam);
                        n.normalize();

#if PRINT
                        std::cerr << "n: " << n << std::endl;
#endif

                        Vec3f viewDirection = -pt;
                        viewDirection.normalize();

#if PRINT
                        std::cerr << "viewDirection " << viewDirection << std::endl;
#endif

                        float nDotView = std::max(0.f, n.dot(viewDirection));

                        const int M = 10; 
                        float checker = (fmod(st.x * M, 1.0) > 0.5) ^ (fmod(st.y * M, 1.0) < 0.5);
                        float c = 0.3 * (1 - checker) + 0.7 * checker;
                        nDotView *= c;

#if PRINT
                        printf("%f c= %f\n", nDotView, c);
#endif

                        uint8_t color = (uint8_t)(nDotView * 255);

                        frameBuffer[y * imageWidth + x] = {color, color, color};

#if PRINT
                        printf("color= %d\n", frameBuffer[y * imageWidth + x].x);
#endif

                    }
                }
            }
        }
    }

    std::ofstream ofs;
    ofs.open("./output.ppm");
    ofs << "P6\n" << imageWidth << " " << imageHeight << "\n255\n";
    ofs.write((char*)frameBuffer, imageWidth * imageHeight * 3);
    ofs.close();

    printf("Depth range : %f %f\n", zSmallest, zLargest);
    for (int i = 0; i < imageHeight * imageWidth; i++)
    {
        /* Reset framebuffer and depth. */
        if (depthBuffer[i] > 255)
        {
            frameBuffer[i] = Vec3<uint8_t>(255);
        }
        else
        {
            frameBuffer[i] = Vec3<uint8_t>(((depthBuffer[i] - zSmallest)/ (zLargest - zSmallest)) * 255);
        }
    }

    ofs.open("./depth.ppm");
    ofs << "P6\n" << imageWidth << " " << imageHeight << "\n255\n";
    ofs.write((char*)frameBuffer, imageWidth * imageHeight * 3);
    ofs.close();


    delete [] frameBuffer;
    delete [] depthBuffer;

    return 0;
}
