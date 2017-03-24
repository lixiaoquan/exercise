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

//#define CULLING 1
//#define SPHERE
//#define DEBUG

using namespace std;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);

float kEpsilon = 1e-8;

bool p = false;

/* Options. */
uint32_t width        = 1024;
uint32_t height       = 747;
float    angleOfView  = 36.87;
uint32_t maxDepth     = 5;
float    bias         = 0.0001;
Mat44f cameraToWorld;


const float kInfinity = std::numeric_limits<float>::max();
static const Vec3f kDefaultBackgroundColor = Vec3f(0.235294, 0.67451, 0.843137);

enum MaterialType { kPhong };
enum RayType {kPrimaryRay, kShadowRay};

Mat44f kIdentity = Mat44f();

float clamp(const float &a, const float &b, const float &c)
{
    return max(a, min(b, c));
}

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{
    return a * (1 - mixValue) + b * mixValue;
}

bool solveQuadratic(float a, float b, float c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;

    if (discr < 0)
    {
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

float deg2rad( float deg)
{
    return deg * M_PI / 180;
}

inline float modulo(const float &f)
{
    return f - std::floor(f);
}



bool rayTriangleIntersect(
    const Vec3f &orig,
    const Vec3f &dir,
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


    float invDet = 1/ det;

    Vec3f tvec = orig - v0;

    u = tvec.dot(pvec) * invDet;

#ifdef DEBUG
    if (p)
    {
     //   cout << " u " << u << endl;
    }
#endif


    if (u < 0 || u > 1)
    {
        return false;
    }

    Vec3f qvec = tvec.crossProduct(v0v1);

    v = dir.dot(qvec) * invDet;

    if (v < 0 || u + v > 1) return false;

    t = v0v2.dot(qvec) * invDet;


#ifdef DEBUG
    if (p)
    {
      //  cout << " t " << t << endl;
    }
#endif

    return (t > 0) ? true : false;
}

class Object
{
public:
    Object(Mat44f &ObjectToWorld) : objectToWorld(ObjectToWorld), worldToObject(objectToWorld.inverse()) {}

    virtual ~Object() {}

    virtual bool intersect(
        const Vec3f &orig,
        const Vec3f &dir,
        float &t,
        uint32_t &,
        Vec2f &
    ) const = 0;

    virtual void getSurfaceData(
        const Vec3f &hitPoint,
        const Vec3f &viewDirection,
        const uint32_t &triIndex,
        Vec2f &uv,
        Vec3f &hitNormal,
        Vec2f &hitTextureCoordinates
    ) const = 0;

    Vec3f   color;
    Mat44f  objectToWorld, worldToObject;
    float   ior = 1;
    Vec3f   albedo = 0.18;
    bool    smoothShading = true;
    MaterialType type = kPhong;
    float   Kd = 0.8;
    float   Ks = 0.2;
    float   n = 10;
};

class Light
{
public:
    Light(const Mat44f &LightToWorld, const Vec3f &c = 1, const float &i = 1) : lightToWorld(LightToWorld), color(c), intensity(i)  {}
    virtual ~Light() {}

    virtual void illumiate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const = 0;

    Vec3f color;
    float intensity;
    Mat44f lightToWorld;
};

class DistantLight : public Light
{
    Vec3f dir;

public:
    DistantLight(const Mat44f &LightToWorld, const Vec3f &c = 1, const float &i = 1) : Light(LightToWorld, c, i) {
        lightToWorld.multDirMatrix(Vec3f(0, 0, -1), dir);
        dir.normalize();
    }

    void illumiate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const
    {
        lightDir = dir;
        lightIntensity = color * intensity;
        distance = kInfinity;
    }
};

class Sphere : public Object
{
public:
    Sphere(Mat44f &ObjectToWorld, const float &r): Object(ObjectToWorld), radius(r), radius2(r * r) {
        ObjectToWorld.multiple(Vec3f(0), center);
        cerr << center << endl;
        cerr << r << endl;
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &triIndex, Vec2f &uv) const
    {
        Vec3f L = orig - center;
        float a = dir.dot(dir);
        float b = 2 * dir.dot(L);
        float c = L.dot(L) - radius2;
        float t0, t1;

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

    void getSurfaceData(
        const Vec3f &hitPoint,
        const Vec3f &viewDirection,
        const uint32_t &triIndex,
        Vec2f &uv,
        Vec3f &hitNormal,
        Vec2f &tex
        ) const

    { 
        hitNormal = hitPoint - center; 
        hitNormal.normalize(); 
        // In this particular case, the normal is simular to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of Phit.
        // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
        // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
        tex.x = (1 + atan2(hitNormal.z, hitNormal.x) / M_PI) * 0.5; 
        tex.y = acosf(hitNormal.y) / M_PI; 
    } 

private:
    float radius;
    float radius2;
    Vec3f center;
};



class TriangleMash : public Object
{
public:
    TriangleMash(
            Mat44f &ObjectToWorld,
            uint32_t nfaces,
            unique_ptr<uint32_t []> &faceIndex,
            unique_ptr<uint32_t []> &vertsIndex,
            unique_ptr<Vec3f []> &verts,
            unique_ptr<Vec3f []> &normals,
            unique_ptr<Vec2f []> &st
    ):Object(ObjectToWorld), numTris(0)
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
            objectToWorld.multiple(verts[i], P[i]);
        }

        trisIndex = unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);

        uint32_t l = 0;

        N = unique_ptr<Vec3f []>(new Vec3f [numTris * 3]);
        texCoordinates = unique_ptr<Vec2f []>(new Vec2f [numTris * 3]);

        Mat44f transformNormals = worldToObject.transpose();

#ifdef DEBUG
        cout << transformNormals << endl;
#endif

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
    //            cout << normals[k] << normals[k + j + 1] << normals[k + j + 2] << endl;
                transformNormals.multDirMatrix(normals[k], N[l]);
                transformNormals.multDirMatrix(normals[k + j + 1], N[l + 1]);
                transformNormals.multDirMatrix(normals[k + j + 2], N[l + 2]);
              //  cout << N[l] << N[l + 1] << N[l + 2] << endl;
                N[l].normalize();
                N[l + 1].normalize();
                N[l + 2].normalize();
              //  cout << N[l] << N[l + 1] << N[l + 2] << endl;
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

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &triIndex, Vec2f &uv) const
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

    void getSurfaceData(
        const Vec3f &hitPoint,
        const Vec3f &viewDirection,
        const uint32_t &triIndex,
        Vec2f &uv,
        Vec3f &hitNormal,
        Vec2f &hitTextureCoordinates
        ) const
    {

        if (smoothShading)
        {
            Vec3f &n0 = N[triIndex * 3];
            Vec3f &n1 = N[triIndex * 3 + 1];
            Vec3f &n2 = N[triIndex * 3 + 2];
            hitNormal = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2;

#ifdef DEBUG
            if (p)
            {
                cout << "get normal" << endl;
                cout << n0 << n1 << n2 << endl;
                cout << "get normal end" << endl;
            }
#endif
        }
        else
        {
            Vec3f v0 = P[trisIndex[triIndex * 3]];
            Vec3f v1 = P[trisIndex[triIndex * 3 + 1]];
            Vec3f v2 = P[trisIndex[triIndex * 3 + 2]];

            hitNormal = (v1 - v0).crossProduct(v2 - v0);
        }

        hitNormal.normalize();

        Vec2f st0 = texCoordinates[triIndex * 3];
        Vec2f st1 = texCoordinates[triIndex * 3 + 1];
        Vec2f st2 = texCoordinates[triIndex * 3 + 2];

        //cerr << st0 << st1 << st2 << uv << endl;
        hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;

    }
    uint32_t numTris;

    unique_ptr<Vec3f []> P;
    unique_ptr<Vec3f []> N;
    unique_ptr<uint32_t []> trisIndex;
    unique_ptr<Vec2f []> texCoordinates;

};

TriangleMash *loadPolyMeshFromFile(const char * file, Mat44f & ObjectToWorld)
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

    return new TriangleMash(ObjectToWorld, numFaces, faceIndex, vertsIndex, verts, normals, st);

    }
    catch(...)
    {
        cerr << "Error happens when reading file" << endl;
    }

    ifs.close();

    return nullptr;
}

#if 0
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
#endif

bool trace(
    const Vec3f &orig,
    const Vec3f &dir,
    vector<unique_ptr<Object>> &objects,
    float &tNear,
    uint32_t &index,
    Vec2f &uv,
    Object *&hitObject,
    RayType Type = kPrimaryRay
    )
{
    vector<unique_ptr<Object>>::const_iterator iter = objects.begin();

    hitObject = nullptr;

#ifdef DEBUG
    if (Type == kShadowRay)
    {
        //cout << orig << " " << dir << " " << objects.size() << " " << tNear << " " << endl;

        //p = true;
    }
#endif

    /* Check all objects. */
    for (; iter != objects.end(); iter++)
    {
        float t = kInfinity;
        uint32_t indexTriangle;
        Vec2f uvTriangle;

        if ((*iter)->intersect(orig, dir, t, indexTriangle, uvTriangle) && t < tNear)
        {
#if 0
            if (Type == kShadowRay && (*iter)->type == kReflectionAndRefraction)
            {
                continue;
            }
#endif

            hitObject = iter->get();
            tNear = t;
            index = indexTriangle;
            uv = uvTriangle;
        }
    }

    p = false;
    return hitObject != nullptr;
}

Vec3f reflect(const Vec3f &I, Vec3f &N)
{
    return I - 2 * I.dot(N) * N;
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, I.dot(N));
    float etai = 1, etat = ior;
    Vec3f n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}

void fresnel(const Vec3f &I, const Vec3f &N, const float &ior, float &kr) 
{
    float cosi = clamp(-1, 1, I.dot(N));
    float etai = 1, etat = ior;
    if (cosi > 0) { std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

/* Compute the color. */
Vec3f castRay(
    const Vec3f &orig,
    Vec3f &dir,
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Light>> &lights,
    uint32_t depth
    )
{
    Vec3f hitColor = 0;
    Object *hitObject = nullptr;
    float t = kInfinity;
    uint32_t index = 0;
    Vec2f uv;

    if (depth > maxDepth)
    {
        return kDefaultBackgroundColor;
    }

    if (trace(orig, dir, objects, t, index, uv, hitObject))
    {
#ifdef DEBUG
       // cout << orig << " " << dir << " " << t << endl;
#endif
        Vec3f hitPoint = orig + dir * t;
        Vec3f hitNormal;
        Vec2f tex;

#ifdef DEBUG
        if (hitObject->type == kReflectionAndRefraction)
        {
            cout << "index " << index << " uv " << uv << endl;
            p = true;
        }
#endif

        hitObject->getSurfaceData(hitPoint, dir, index, uv, hitNormal, tex);
#ifdef DEBUG
        //cout << "hitPoint"<< hitPoint << " hitNormal " << hitNormal<< endl;
#endif

        Vec3f diffuse = 0, specular = 0;

        /************
         * Shading *
        *************/
        switch (hitObject->type)
        {
            case kPhong:
                for (uint32_t i = 0; i < lights.size(); i++)
                {
                    Vec3f lightDir, lightIntensity;
                    Object *hitShadow = nullptr;

                    lights[i]->illumiate(hitPoint, lightDir, lightIntensity, t);

                    //cout << "hitPoint " << hitPoint << endl;

                    bool vis = !trace(hitPoint + hitNormal * bias, -lightDir, objects, t, index, uv, hitShadow, kShadowRay);

                    diffuse += vis * hitObject->albedo * lightIntensity * max(0.f, hitNormal.dot(-lightDir));

                    //cout << "albedo " << hitObject->albedo << " lightIntensity " << lightIntensity << endl;

                    Vec3f R = reflect(lightDir, hitNormal);

                    specular += vis * lightIntensity * pow(max(0.f, R.dot(-dir)), hitObject->n);
                }

                hitColor = diffuse * hitObject->Kd + specular * hitObject->Ks;
                break;
            default:
                break;
            }
    }
    else
    {
        hitColor = kDefaultBackgroundColor;
    }

    return hitColor;
}

void render (
    vector<unique_ptr<Object>> &objects,
    vector<unique_ptr<Light>> &lights
    )
{
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
    ofs.open("./phong.ppm");
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
    vector<unique_ptr<Light>> lights;

    bias = 0.0001;
    cameraToWorld[3][2] = 12;
    cameraToWorld[3][1] = 1;
    
    Matrix44f xform;
    xform[0][0] = 1;
    xform[1][1] = 1;
    xform[2][2] = 1;
    TriangleMash *mesh = loadPolyMeshFromFile("geometry/plane.geo", xform);
    if (mesh != nullptr) {
        mesh->smoothShading = false;
        objects.push_back(std::unique_ptr<Object>(mesh));
    }

    /* Balls. */
    float w[5] = {0.04, 0.08, 0.1, 0.15, 0.2};
    for (int i = -4, n = 2, k = 0; i <= 4; i+= 2, n *= 5, k++) {
        Matrix44f xformSphere;
        xformSphere[3][0] = i;
        xformSphere[3][1] = 1;
        Sphere *sph = new Sphere(xformSphere, 0.9);
        sph->n = n;
        sph->Ks = w[k];
        objects.push_back(std::unique_ptr<Object>(sph));
    }


    Matrix44f l2w(11.146836, -5.781569, -0.0605886, 0, -1.902827, -3.543982, -11.895445, 0, 5.459804, 10.568624, -4.02205, 0, 0, 0, 0, 1);

    lights.push_back(unique_ptr<Light>(new DistantLight(l2w, 1, 5)));

    auto timeStart = chrono::high_resolution_clock::now();

    render(objects, lights);

    auto timeEnd = chrono::high_resolution_clock::now();
    auto passedTime = chrono::duration<double, std::milli>(timeEnd - timeStart).count();

    cerr << mesh->numTris << " " << passedTime << endl;
}

