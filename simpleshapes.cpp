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

#ifndef M_PI
#define M_PI 3.14159265
#endif


using namespace std;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);


const float kInfinity = std::numeric_limits<float>::max();

float clamp(const float &a, const float &b, const float &c)
{
    return max(a, min(b, c));
}

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

bool solveQuadratic(float a, float b, float c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;

    if (discr < 0)
    {
    //    printf("%s %d %f\n", __func__, __LINE__, discr);
        return false;
    }

    if (discr == 0)
    {
        x0 = x1 = -b / (2 * a);
    }
    else
    {
        float q = (b > 0) ?
                -0.5 * (b + sqrt(discr)) :
                -0.5 * (b - sqrt(discr)) ;

        x0 = q/a;
        x1 = c/q;
    }

    if (x0 > x1)
    {
        swap(x0, x1);
    }

    return true;
}


class Object
{
public:
    Object() : color(dis(gen), dis(gen), dis(gen)) {}
    virtual ~Object() {}
    virtual bool intersect(Vec3f &orig, Vec3f &dir, float &t) const = 0;
    virtual void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const = 0;

    Vec3f color;
};

class Light
{
public:
    Light() {}
    virtual ~Light() {}
};

class Sphere : public Object
{
public:
    Sphere(Vec3f& c, float &r): radius(r), radius2(r * r), center(c) {
        cerr << c << endl;
        cerr << r << endl;
    }

    bool intersect(Vec3f &orig, Vec3f &dir, float &t) const
    {
        Vec3f L = orig - center;
        float a = dir.dot(dir);
        float b = 2 * dir.dot(L);
        float c = L.dot(L) - radius2;
        float t0, t1;

   //     cerr << " L " << L <<"dir" << dir << endl;
   //     printf("%s %d %f %f %f \n", __func__, __LINE__, a, b, c);

   //     exit(0);

        if (!solveQuadratic(a, b, c, t0, t1))
        {
            return false;
        }

        if (t0 < 0)
        {
            t0 = t1;

            if (t0 < 0)
            {
                return false;
            }
        }

        t = t0;

        return true;
    }

    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const 
    { 
        Nhit = Phit - center; 
        Nhit.normalize(); 
        // In this particular case, the normal is simular to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of Phit.
        // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
        // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5; 
        tex.y = acosf(Nhit.y) / M_PI; 
    } 

private:
    float radius;
    float radius2;
    Vec3f center;
};


uint32_t width = 640;
uint32_t height = 480;
float angleOfView  = 51.52;

bool trace(
    Vec3f &orig,
    Vec3f &dir,
    vector<unique_ptr<Object>> &objects,
    float &tNear,
    const Object *&hitObject
    )
{
    tNear = kInfinity;

    vector<unique_ptr<Object>>::const_iterator iter = objects.begin();

    /* Check all objects. */
    for (; iter != objects.end(); iter++)
    {
        float t = kInfinity;

        if ((*iter)->intersect(orig, dir, t) && t < tNear)
        {
            hitObject = iter->get();
            tNear = t;
        }
    }

    return hitObject != nullptr;
}

/* Compute the color. */
Vec3f castRay(
    Vec3f &orig,
    Vec3f &dir,
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Object>> &lights,
    uint32_t depth
    )
{
    Vec3f hitColor = 0;
    const Object *hitObject = nullptr;
    float t;

    if (trace(orig, dir, objects, t, hitObject))
    {
        Vec3f Phit = orig + dir * t;
        Vec3f Nhit;
        Vec2f tex;
        hitObject->getSurfaceData(Phit, Nhit, tex);
        Vec3f dir1 = -dir;

        float scale = 4;

        float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5); 
        //hitColor = max(0.f, Nhit.dot(dir1)) * mix(hitObject->color, hitObject->color * 0.8, pattern);
        hitColor =  mix(hitObject->color, hitObject->color * 0.8, pattern) * max(0.f, Nhit.dot(dir1));
    }

    return hitColor;
}

void render (
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Object>> &lights
    )
{
    Mat44f cameraToWorld = Mat44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

    float imageAspectRatio = width / (float)height;
    float scale = tan(angleOfView / 2 * M_PI / 180);


    Vec3f minWorld(kInfinity), maxWorld(-kInfinity);

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

            buffer[j * width + i] = castRay(orig, dir, objects, lights, 0);
        }
    }


    // save to file
    std::ofstream ofs;
    ofs.open("./simpleshapes.ppm");
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
    vector<unique_ptr<Object>> objects;
    vector<unique_ptr<Object>> lights;

    // generate a scene made of random spheres
    uint32_t numSpheres = 35;
    gen.seed(0);
    for (uint32_t i = 0; i < numSpheres; ++i) {
        Vec3f randPos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, (0.5 + dis(gen) * 10));
        float randRadius = (0.5 + dis(gen) * 0.5);
        objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));
    }

    render(objects, lights);
}

