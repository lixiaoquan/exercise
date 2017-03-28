/* Follow code in Monte Carlo Simulation. */
#include <iostream> 
#include <cstdlib> 
#include <cstdio> 
#include <fstream> 
#include <cmath> 
#include <algorithm> 
#include <iomanip> 
#include <limits> 
#include <vector> 
#include "stdlib.h"
#include <random> 
#include "stdint.h"
#include "string.h"

#include "geometry.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.1415926
#endif

double getCosTheta(const double &g) // sampling the H-G scattering phase function
{
    if (g == 0) return 2 * drand48() - 1;
    double mu = (1 - g * g) / (1 - g + 2 * g * drand48());
    return (1 + g * g - mu * mu) / (2.0 * g);
}

void spin(double &mu_x, double &mu_y, double &mu_z, const double &g)
{
    double costheta = getCosTheta(g);
    double phi = 2 * M_PI * drand48();
    double sintheta = sqrt(1.0 - costheta * costheta); // sin(theta)
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    if (mu_z == 1.0) {
        mu_x = sintheta * cosphi;
        mu_y = sintheta * sinphi;
        mu_z = costheta;
    }
    else if (mu_z == -1.0) {
        mu_x = sintheta * cosphi;
        mu_y = -sintheta * sinphi;
        mu_z = -costheta;
    }
    else {
        double denom = sqrt(1.0 - mu_z * mu_z);
        double muzcosphi = mu_z * cosphi;
        double ux = sintheta * (mu_x * muzcosphi - mu_y * sinphi) / denom + mu_x * costheta;
        double uy = sintheta * (mu_y * muzcosphi + mu_x * sinphi) / denom + mu_y * costheta;
        double uz = -denom * sintheta * cosphi + mu_z * costheta;
        mu_x = ux, mu_y = uy, mu_z = uz;
    }
}


void MCSimulation(double * records, double * reflects, double *absorption, uint32_t size, uint32_t pIndex)
{
    uint32_t nphotons = 100000; //sim that much of nphotons.

    /* Beer's Law Parameters. */
    double sigma_a = 1, sigma_s = 2;

    if (pIndex == 0)
    {
        sigma_a = 1;
        sigma_s = 2;
    }
    else if (pIndex == 1)
    {
        sigma_a = 1.2;
        sigma_s = 2.4;
    }
    else
    {
        sigma_a = 1.4;
        sigma_s = 2.8;
    }

    double sigma_t = sigma_a + sigma_s;

    double Rd = 0, Tt = 0;

    double d = 0.5, slabsize = 0.5, g = 0.75;

    const short m = 10;

    for (int n = 0; n < nphotons; n++)
    {
        double w = 1;
        double x = 0, y = 0, z = 0, mux = 0, muy = 0, muz = 1;

        while (w != 0)
        {
            /* Random walking a distance. */
            double s = -log(drand48()) / sigma_t;

            double distToBoundary = 0;

            if (muz > 0)
            {
                distToBoundary = (d - z) / muz;
            }
            else if (muz < 0)
            {
                distToBoundary = (- z) / muz;
            }

            if (s > distToBoundary)
            {
                if (muz > 0) Tt += w; else Rd += w;

                int xi = (int)((x + slabsize / 2) / slabsize * size);
                int yi = (int)((y + slabsize / 2) / slabsize * size);
                if (muz > 0 && xi >= 0 && x < size && yi >= 0 && yi < size) {
                    records[yi * size + xi] += w;
                }

                if (muz < 0 && xi >= 0 && x < size && yi >= 0 && yi < size) {
                    reflects[yi * size + xi] += w;
                }

                break;
            }

            x += s * mux;
            y += s * muy;
            z += s * muz;

            double dw = sigma_a / sigma_t;

            w -= dw;

            w = max(0.0, w);

            if (w < 0.001)
            {
                if (drand48() > 1.0 / m)
                {
                    int xi = (int)((x + slabsize / 2) / slabsize * size);
                    int yi = (int)((y + slabsize / 2) / slabsize * size);
                    if (xi >= 0 && x < size && yi >= 0 && yi < size) {
                        absorption[yi * size + xi] += (1-w);
                    }

                    break;
                }
                else w *= m;
            }

            spin(mux, muy, muz, g);
        }
    }

#if 0
    printf("%f %f\n", Rd/ nphotons, Tt / nphotons);
#endif

}

int main()
{
    const uint32_t size = 512;
    double *records[3];
    float *pixels[3];

    for (int i = 0; i < 3; i++)
    {
        records[i] = new double[size * size * 3];
        memset(records[i], 0x0, sizeof(double) * size * size * 3);
        pixels[i] = new float[size * size * 3];
    }


    for (int c = 0; c < 3; c++)
    {
        uint32_t npass = 1;

        while (npass < 64)
        {
            MCSimulation(records[0], records[1], records[2], size, c);

            npass++;
        }

        for (int i = 0; i < size * size ; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                pixels[j][i * 3 + c] = records[j][i] / (npass - 1);
            }
        }
    }

    for (int j = 0; j < 3; j++)
    {
        char name[100];

        sprintf(name, "./mcsim_%d.ppm", j);
        // save image to file
        std::ofstream ofs;
        ofs.open(name, std::ios::out | std::ios::binary);
        ofs << "P6\n" << size << " " << size << "\n255\n";
        for (uint32_t i = 0; i < size * size; ++i) {
            unsigned char r = (unsigned char)(255 * std::min(1.0f, pixels[j][i * 3]));
            unsigned char g = (unsigned char)(255 * std::min(1.0f, pixels[j][i * 3 + 1]));
            unsigned char b = (unsigned char)(255 * std::min(1.0f, pixels[j][i * 3 + 2]));
            ofs << r << g << b;
        }

        ofs.close();
    }


    for (int j = 0; j < 3; j++)
    {
        delete [] records[j];
        delete [] pixels[j];
    }

    return 0;
}
