#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DFT_int 10
#define dt 0.01        //Time step (s)
#define t_s 500

int main()
{
    int i, j = 0, t;
    int DFT_size = t_s / DFT_int;
    double f = 1;
    double w = 2 * M_PI * f;
    double sin_sum, cos_sum, amplitude;
    double mean = 0;
    static double sample[t_s / DFT_int];
    FILE *stream5;

    stream5 = fopen("DFT_test.csv", "w");

    fprintf(stream5, "angular frequency, Im (sin), Re (cos), Magnitude, time, sample\n");

    for (i = 0; i < DFT_size; i++)
    {
        sample[i] = 0;
    }

    for (t = 0; t <= t_s; t++)
    {
        printf("Start adding sample %i\n", t);
        sample[j] += sin(w * t * dt);
        printf("%lf\n", sample[j]);
        printf("Done adding sample %i\n", t);

        if (t % DFT_int == 0 && t != 0)
        {
            sample[j] = sample[j] / DFT_int;
            j++;
        }
    }

    //This component is to remove the peak at frequency 0 (which represent the direct current / constant offset, and is meaningless in our case)
    for (i = 0; i < DFT_size; i++)
    {
        mean += sample[i];
    }
    mean = mean / DFT_size;
    for (i = 0; i < DFT_size; i++)
    {
        sample[i] -= mean;
    }

    //Performs DFT
    for (i = 0; i < DFT_size; i++)
    {
        sin_sum = 0;
        cos_sum = 0;
        amplitude = 0;
        for (t = 0; t < DFT_size; t++)
        {
            sin_sum += -sample[t] * sin(i * 2 * M_PI * t / DFT_size);
            cos_sum += sample[t] * cos(i * 2 * M_PI * t / DFT_size);
        }
        sin_sum /= DFT_size;
        cos_sum /= DFT_size;
        amplitude = sqrt(sin_sum * sin_sum + cos_sum * cos_sum);
        fprintf(stream5, "%lf, %lf, %lf, %lf, %lf, %lf\n", i * (1 / (DFT_int * dt)) / DFT_size, sin_sum, cos_sum, amplitude, i * DFT_int * dt, sample[i]);
    }

    fclose(stream5);
    return 0;
}
