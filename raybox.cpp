/* Follow code in ray box. */
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

using namespace std;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);

class Ray {
public:
    Ray(Vec3f orig, Vec3f dir) : orig(orig), dir(dir)
    {
        invdir = 1 / dir;
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
    }
    Vec3f orig, dir;
    Vec3f invdir;
    int sign[3];
};

class AABBox {

public:

    AABBox(Vec3f b1, Vec3f b2) {
        bounds[0] = b1;
        bounds[1] = b2;
    }

    bool intersect(const Ray &r, float &t) const
    {
        float tmin, tmax, tymin, tymax;

        tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x;
        tmax = (bounds[1- r.sign[0]].x - r.orig.x) * r.invdir.x;
        tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y;
        tymax = (bounds[1- r.sign[1]].y - r.orig.y) * r.invdir.y;

        if ((tmin > tymax) || (tymin > tmax))
            return false;

        if (tymin > tmin)
            tmin = tymin;

        if (tymax < tmax)
            tmax = tymax;

        float tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z;
        float tzmax = (bounds[1- r.sign[2]].z - r.orig.z) * r.invdir.z;

        if ((tmin > tzmax) || (tzmin > tmax))
            return false;

        if (tzmin > tmin)
            tmin = tzmin;

        if (tzmax < tmax)
            tmax = tzmax;

        t = tmin;

        return true;

    }

    Vec3f bounds[2];

};

int main()
{
    AABBox box(Vec3f(-1), Vec3f(1));

    gen.seed(0);

    for (int i = 0; i < 32; i++)
    {
        Vec3f dir(2 * dis(gen) - 1, 2 * dis(gen) - 1, 2 * dis(gen) - 1);
        dir.normalize();

        Ray r(Vec3f(0), dir);

        float t;

        if (box.intersect(r, t))
        {
            Vec3f h = r.orig + dir * t;

            cerr << r.orig << " " << h << endl;
        }
    }
    return 0;
}
