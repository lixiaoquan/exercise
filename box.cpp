/* A sample from A Gentle Introduction to Computer Graphics*/

#include <iostream>
#include <iomanip>

using namespace std;

#define VERTICIES_NUM 8
#define WIDTH 512
#define HEIGHT 512

int main()
{
    float vertices[VERTICIES_NUM][3] = {
        {1, -1, -5},
        {1, -1, -3},
        {1,  1, -5},
        {1,  1, -3},
        {-1, -1, -5},
        {-1, -1, -3},
        {-1,  1, -5},
        {-1,  1, -3},
    };

    for (int i = 0; i < VERTICIES_NUM; i++)
    {
        /* prospective divide. */
        float iz = - 1 / vertices[i][2];
        float screen_x = vertices[i][0] * iz;
        float screen_y = vertices[i][1] * iz;

        /* NDC. */
        float ndc_x = (screen_x + 1) / 2;
#if 0
        float ndc_y = (1 - screen_y) / 2;
#else
        float ndc_y = (1 + screen_y) / 2;
#endif

        /* Raster coordinate. */
        float raster_x = ndc_x * WIDTH;
        float raster_y = ndc_y * HEIGHT;

        cout << fixed << setprecision(6) << raster_x << " " << raster_y << "\n";
    } 

    return 0;
}
