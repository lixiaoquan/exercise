/* Follow code in Camera rays. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <iomanip> 
#include <limits> 
#include <vector> 
#include "stdint.h"
#include "string.h"

#include "geometry.h"
#include "vertexdata.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif

using namespace std;

float clamp(const float &a, const float &b, const float &c)
{
    return max(a, min(b, c));
}


class Object
{
public:
    Object() {}
    virtual ~Object() {}
    virtual bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const = 0;
};

class Light
{
public:
    Light() {}
    virtual ~Light() {}
};

class Sephere : public Object
{
public:
    Sephere(Vec3f& c, float &r): radius(r), radius2(r * r), center(r) {}

    bool instersect(const Vec3f &orij, const Vec3f &dir, float &t)
    {
        return false;
    }

private:
    float radius;
    float radius2;
    float center;
};


uint32_t width = 640;
uint32_t height = 480;
float angleOfView  = 90;

Vec3f castRay(
    Vec3f &orig,
    Vec3f &dir,
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Object>> &lights,
    uint32_t depth
    )
{
    Vec3f hitColor = (dir + Vec3f(1)) * 0.5;
    return hitColor;
}

void render (
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Object>> &lights
    )
{
    Mat44f cameraToWorld;

    float imageAspectRatio = width / (float)height;
    float scale = tan(angleOfView / 2 * M_PI / 180);

    const float kInfinity = std::numeric_limits<float>::max();

    Vec3f minWorld(kInfinity), maxWorld(-kInfinity);

    Vec3f *buffer = new Vec3f[width * height];

    memset(buffer, 0, width * height);

    Vec3f orig;

    cameraToWorld.multiple(Vec3f(0), orig);

    for (uint32_t j = 0; j < height; j++)
    {
        for (uint32_t i = 0; i < width; i++)
        {
            /* NDC */

            float x = (2 * (i + 0.5) / (float)width - 1) * imageAspectRatio * scale;
            float y = (2 * (j + 0.5) / (float)height - 1) * scale;

            Vec3f dir;

            cameraToWorld.multiple(Vec3f(x, y, -1), dir);
            dir.normalize();

            buffer[j * width + i] = castRay(orig, dir, objects, lights, 0);
        }
    }


    // save to file
    std::ofstream ofs;
    ofs.open("./camerarays.ppm");
    ofs << "P5\n" << width << " " << height << "\n255\n";


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
    vector<unique_ptr<Object>> objects;
    vector<unique_ptr<Object>> lights;
    render(objects, lights);
}

