/* Follow code in Camera rays. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <sstream> 
#include <cmath> 
#include <iomanip> 
#include <limits> 
#include <vector> 
#include <random> 
#include <chrono> 
#include "stdint.h"
#include "string.h"

#include "geometry.h"
#include "vertexdata.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif

#define CULLING 1
//#define SPHERE

using namespace std;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);

float kEpsilon = 1e-8;

const float kInfinity = std::numeric_limits<float>::max();
static const Vec3f kDefaultBackgroundColor = Vec3f(0.235294, 0.67451, 0.843137);


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
    virtual bool intersect(Vec3f &orig, Vec3f &dir, float &t, uint32_t &, Vec2f &) const = 0;
    virtual void getSurfaceData(const uint32_t &triIndex, Vec2f &uv, Vec3f &hitNormal, Vec2f &hitTextureCoordinates) const = 0;

    Vec3f color;
};

class Light
{
public:
    Light() {}
    virtual ~Light() {}
};

uint32_t width = 640;
uint32_t height = 480;
float angleOfView  = 50.0393;

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

class TriangleMash : public Object
{
public:
    TriangleMash(
            uint32_t nfaces,
            unique_ptr<uint32_t []> &faceIndex,
            unique_ptr<uint32_t []> &vertsIndex,
            unique_ptr<Vec3f []> &verts,
            unique_ptr<Vec3f []> &normals,
            unique_ptr<Vec2f []> &st
    ):numTris(0)
    {

        uint32_t k = 0, maxVerIndex = 0;

        for (int i = 0; i < nfaces; i++)
        {
            numTris += faceIndex[i] - 2;

            for (int j = 0; j < faceIndex[i]; j++)
            {
                maxVerIndex = max(vertsIndex[k + j], maxVerIndex);
            }

            k += faceIndex[i];
        }


        maxVerIndex += 1;

        //cerr << numTris << " " << maxVerIndex << endl;

        P = unique_ptr<Vec3f []> (new Vec3f[maxVerIndex]);

        for (int i = 0; i < maxVerIndex; i++)
        {
            P[i] = verts[i];
        }

        trisIndex = unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);

        uint32_t l = 0;

        N = unique_ptr<Vec3f []>(new Vec3f [numTris * 3]);
        texCoordinates = unique_ptr<Vec2f []>(new Vec2f [numTris * 3]);

        for (int i = 0, k = 0; i < nfaces; i++)
        {
            for (int j = 0; j < faceIndex[i] - 2; j++)
            {
                trisIndex[l] = vertsIndex[k];
                trisIndex[l + 1] = vertsIndex[k + j + 1];
                trisIndex[l + 2] = vertsIndex[k + j + 2];
#ifdef SPHERE
                N[l] = normals[vertsIndex[k]];
                N[l + 1] = normals[vertsIndex[k + j + 1]];
                N[l + 2] = normals[vertsIndex[k + j + 2]];
                texCoordinates[l] = st[vertsIndex[k]];
                texCoordinates[l + 1] = st[vertsIndex[k + j + 1]];
                texCoordinates[l + 2] = st[vertsIndex[k + j + 2]];
#else
               N[l] = normals[k];
                N[l + 1] = normals[k + j + 1];
                N[l + 2] = normals[k + j + 2];
                texCoordinates[l] = st[k];
                texCoordinates[l + 1] = st[k + j + 1];
                texCoordinates[l + 2] = st[k + j + 2];

#endif

                //std::cerr << l << " " <<k << " " << j <<  " " << texCoordinates[l] << texCoordinates[l + 1] << texCoordinates[l + 2] << std::endl;


                l += 3;
            }

            k += faceIndex[i];
        }

    };

    bool intersect(Vec3f &orig, Vec3f &dir, float &tNear, uint32_t &triIndex, Vec2f &uv) const
    {
        uint32_t j = 0;

        bool isect = false;

        for (int i = 0; i < numTris; i++)
        {
            Vec3f &v0 = P[trisIndex[j]];
            Vec3f &v1 = P[trisIndex[j + 1]];
            Vec3f &v2 = P[trisIndex[j + 2]];

            float t = kInfinity, u, v;

            if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v) && t < tNear)
            {
                tNear = t;
                uv.x = u;
                uv.y = v;
                triIndex = i;
                isect = true;
            }

            j += 3;
        }

        return isect;
    }

    void getSurfaceData(const uint32_t &triIndex, Vec2f &uv, Vec3f &hitNormal, Vec2f &hitTextureCoordinates) const
    {

        Vec3f v0 = P[trisIndex[triIndex * 3]];
        Vec3f v1 = P[trisIndex[triIndex * 3 + 1]];
        Vec3f v2 = P[trisIndex[triIndex * 3 + 2]];

        hitNormal = (v1 - v0).crossProduct(v2 - v0);
        hitNormal.normalize();

        Vec2f st0 = texCoordinates[triIndex * 3];
        Vec2f st1 = texCoordinates[triIndex * 3 + 1];
        Vec2f st2 = texCoordinates[triIndex * 3 + 2];

        //cerr << st0 << st1 << st2 << uv << endl;
        hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;

#ifndef SPHERE
        Vec3f &n0 = N[triIndex * 3];
        Vec3f &n1 = N[triIndex * 3 + 1];
        Vec3f &n2 = N[triIndex * 3 + 2];
        hitNormal = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2;
#endif

    }
    uint32_t numTris;

    unique_ptr<Vec3f []> P;
    unique_ptr<Vec3f []> N;
    unique_ptr<uint32_t []> trisIndex;
    unique_ptr<Vec2f []> texCoordinates;

};

TriangleMash *loadPolyMeshFromFile(const char * file)
{
    ifstream ifs;

    try
    {

        ifs.open(file);

        if (ifs.fail())
        {
            throw;
        }

        stringstream ss;
        ss << ifs.rdbuf();

        uint32_t numFaces;

        ss >> numFaces;

        unique_ptr<uint32_t []> faceIndex(new uint32_t[numFaces]);

        uint32_t vertsIndexArraySize = 0;

        for (int i = 0; i < numFaces; i++)
        {
            ss >> faceIndex[i];

            vertsIndexArraySize += faceIndex[i];
        }


    std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[vertsIndexArraySize]);
    uint32_t vertsArraySize = 0; // reading vertex index array 
    for (uint32_t i = 0; i < vertsIndexArraySize; ++i) { 
        ss >> vertsIndex[i]; 
        if (vertsIndex[i] > vertsArraySize) 
        vertsArraySize = vertsIndex[i]; 
        } 
        
        vertsArraySize += 1; // reading vertices 
        std::unique_ptr<Vec3f []> verts(new Vec3f[vertsArraySize]); 
        for (uint32_t i = 0; i < vertsArraySize; ++i) 
        { ss >> verts[i].x >> verts[i].y >> verts[i].z; } // reading normals 
        
        std::unique_ptr<Vec3f []> normals(new Vec3f[vertsIndexArraySize]); 
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) { 
        ss >> normals[i].x >> normals[i].y >> normals[i].z; } // reading st coordinates 
        
        std::unique_ptr<Vec2f []> st(new Vec2f[vertsIndexArraySize]); 
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) { ss >> st[i].x >> st[i].y; }

    return new TriangleMash(numFaces, faceIndex, vertsIndex, verts, normals, st);

    }
    catch(...)
    {
    
    }

    ifs.close();

    return nullptr;
}

TriangleMash * generatePolyMesh(float rad, uint32_t divs)
{
    uint32_t numVertices = (divs -1)*divs + 2;

    unique_ptr<Vec3f []> P(new Vec3f[numVertices]);
    unique_ptr<Vec3f []> N(new Vec3f[numVertices]);
    unique_ptr<Vec2f []> st(new Vec2f[numVertices]);

    float u = -M_PI / 2;
    float v = -M_PI;
    float du = M_PI /divs;
    float dv = 2 * M_PI /divs;

    P[0] = N[0] = Vec3f(0, -rad, 0);

    uint32_t k = 1;

    for (uint32_t i = 0; i < divs - 1; i++)
    {
        u += du;
        v = -M_PI;

        for (uint32_t j = 0; j < divs; j++)
        {
            P[k] = N[k] = Vec3f(rad * cos(u) * cos(v), rad * sin(u), rad*cos(u) * sin(v));

            st[k].x = u/M_PI + 0.5;
            st[k].y = v * 0.5 / M_PI + 0.5;

            v += dv;
            k++;
        }
    }

    P[k] = N[k] = Vec3f(0, rad, 0);

    uint32_t npolys = divs * divs;

    unique_ptr<uint32_t []> faceIndex(new uint32_t[npolys]);
    unique_ptr<uint32_t []> vertsIndex(new uint32_t[(6 + (divs - 1) * 4) * divs]);

    // create the connectivity lists
    uint32_t vid = 1, numV = 0, l = 0; 
    k = 0; 
    for (uint32_t i = 0; i < divs; i++) { 
        for (uint32_t j = 0; j < divs; j++) { 
            if (i == 0) { 
                faceIndex[k++] = 3; 
                vertsIndex[l] = 0; 
                vertsIndex[l + 1] = j + vid; 
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid : j + vid + 1; 
                l += 3; 
            } 
            else if (i == (divs - 1)) { 
                faceIndex[k++] = 3; 
                vertsIndex[l] = j + vid + 1 - divs; 
                vertsIndex[l + 1] = vid + 1; 
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs; 
                l += 3; 
            } 
            else { 
                faceIndex[k++] = 4; 
                vertsIndex[l] = j + vid + 1 - divs; 
                vertsIndex[l + 1] = j + vid + 1; 
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 : j + vid + 2; 
                vertsIndex[l + 3] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs; 
                l += 4; 
            } 
            numV++; 
        } 
        vid = numV; 
    } 

    return new TriangleMash(npolys, faceIndex, vertsIndex, P, N, st);
}

bool trace(
    Vec3f &orig,
    Vec3f &dir,
    vector<unique_ptr<Object>> &objects,
    float &tNear,
    uint32_t &index,
    Vec2f &uv,
    Object *&hitObject
    )
{
    vector<unique_ptr<Object>>::const_iterator iter = objects.begin();

    hitObject = nullptr;

    /* Check all objects. */
    for (; iter != objects.end(); iter++)
    {
        float t = kInfinity;
        uint32_t indexTriangle;
        Vec2f uvTriangle;

        if ((*iter)->intersect(orig, dir, t, indexTriangle, uvTriangle) && t < tNear)
        {
            hitObject = iter->get();
            tNear = t;
            index = indexTriangle;
            uv = uvTriangle;

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
    Vec3f hitColor = kDefaultBackgroundColor;
    Object *hitObject = nullptr;
    float t = kInfinity;
    uint32_t index = 0;
    Vec2f uv;

    if (trace(orig, dir, objects, t, index, uv, hitObject))
    {
        Vec3f Phit = orig + dir * t;
        Vec3f hitNormal;
        Vec2f tex;
        hitObject->getSurfaceData(index, uv, hitNormal, tex);
        Vec3f d = -dir;

        float NdotView = max(0.f, hitNormal.dot(d));


        int M = 10;

        float checker = (fmod(tex.x * M, 1.0) > 0.5) ^ (fmod(tex.y * M, 1.0) > 0.5); 
        //hitColor = max(0.f, Nhit.dot(dir1)) * mix(hitObject->color, hitObject->color * 0.8, pattern);
        float c = 0.3 *(1 - checker) + 0.7 * checker;
        //cerr << NdotView << " " << c << " "<< tex <<  endl;
        hitColor = c * NdotView;
    }

    return hitColor;
}

void render (
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Object>> &lights
    )
{
    Mat44f cameraToWorld;

    Mat44f tmp = Mat44f(0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295, 0.624695, 0, -1.63871, -5.747777, -40.400412, 1);

    cameraToWorld = tmp.inverse();


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
    ofs.open("./mesh.ppm");
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
#ifdef SPHERE
    TriangleMash *mesh = generatePolyMesh(10, 100);
#else
    TriangleMash *mesh = loadPolyMeshFromFile("./cow.geo");
#endif
    objects.push_back(unique_ptr<Object>(mesh));




    auto timeStart = chrono::high_resolution_clock::now();

    render(objects, lights);

    auto timeEnd = chrono::high_resolution_clock::now();
    auto passedTime = chrono::duration<double, std::milli>(timeEnd - timeStart).count();

    cerr << mesh->numTris << " " << passedTime << endl;
}

