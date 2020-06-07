/*
CHAMBER RANGE X = -L to 0
OSCILLATOR AT X = -L
SPECULAR WALL AT X = 0
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//Box definition
#define L 100.0       //Length of the square reservoir (angstrom)
#define n 10          //Number of particles in one row
#define N 1000        //Total number of particles (By rho * V_box)
#define rho 7.6458e-8 //Density of air (eV*ps^2/angstrom^5)
//Data recording
#define sum_int 100   //Interval for printing the sum
#define DFT_int 500  //Record DFT samples every 500 frames
//Constant definition
#define kB 8.617333e-5  //Boltzmann constant (eV/K)
#define T 300.0         //Initialization temperature (K)
#define u 9.3149e8      //Atomic mass (eV/c^2)
#define v_c 2997924.58  //Speed of light (angstrom/ps)
#define dt 0.01        //Time step (ps)
#define sigma 3.4       //Zero-potential point (angstrom)
#define sigma2 11.56    //Zero-potential point^2 (angstrom)
#define epsilon 1.03e-2 //Depth of the potential (eV)
#define rc 8.5          //Cut-off distance (angstrom)
//Time definition
#define t_s 100000 //Total simulation time
#define t1 20000 //First transient state (Right after initialization)
#define t2 30000 //Second transient state (After oscillator is on)

void initialize(double[], double[], double[], double[], double[], double[], double[], double[], double[], double[], double[], double[], double[], double[], double[]);
void DFT(double[], double[], double[]);

int main()
{
    printf("starting program\n");
    static double x[N];
    static double y[N];
    static double z[N];
    static double vx[N];
    static double vy[N];
    static double vz[N];
    //f0 = f(t), f1 = f(t+dt)
    static double fx0[N], fy0[N], fz0[N], fx1[N], fy1[N], fz1[N];
    double r2, rx, ry, rz;
    //Potential Energy, Kinetic Energy, average KE of one particle, Total Energy, Temperature during the simulation
    double PE = 0, KE = 0, KE_avg = 0, TE = 0, T_exp = 0;
    //Argon mass, rescaled such that unit is (eV * ps^2 /angstrom^2)
    double mass = 39.948 * u / (v_c * v_c);
    //scale factor for speed of argon molecule (angstrom/ps)
    double v_s = sqrt(kB * T / mass);
    //Angular frequency
    double w = 2 * M_PI * 20 / ((t_s - t1) * dt);
    //double w = M_PI / L;
    double x_temp = 0;
    //DFT
    static double f_wall[(int)((t_s - t2) / DFT_int)];
    static double vx_tsum[(int)((t_s - t2) / DFT_int)];
    static double temperature[(int)((t_s - t2) / DFT_int)];
    //Truncated LJ-potential
    double cutoff = 4 * epsilon * (pow(sigma / rc, 12) - pow(sigma / rc, 6));
    //Others
    int t, i, j, k;
    int h = 0;
    printf("done declaring variables and array\n");
    FILE *stream;
    FILE *stream2;
    FILE *stream4;
    stream = fopen("box.xyz", "w");
    stream2 = fopen("data.csv", "w");
    stream4 = fopen("temp_energy.csv", "w");
    printf("finish set up stream\n");
    if (stream == NULL || stream2 == NULL || stream4 == NULL)
    {
        printf("Error occurred when opening a file.");
    }
    else
    {
        initialize(x, y, z, vx, vy, vz, fx0, fy0, fz0, fx1, fy1, fz1, f_wall, vx_tsum, temperature);

        fprintf(stream, "%i\n", N);
        fprintf(stream, "Initialized as 300K, time = 0 ps\n", t * dt);
        fprintf(stream2, "ps, f_wall, total velocity\n");
        fprintf(stream4, "ps, Total PE, Total KE, Total Energy, KE_avg, T_exp\n");
        for (i = 0; i < N; i++)
        {
            fprintf(stream, "Ar %lf %lf %lf\n", x[i], y[i], z[i]);
        }

        //Find the forces f(t)  just after initialization (Only run once)
        //Updating each particle
        for (i = 0; i < (N - 1); i++)
        {
            //Calculation of f(t) by the truncated LJ-potential
            for (j = i + 1; j < N; j++)
            {
                //Particle doesn't interact with itself
                rx = x[i] - x[j];
                ry = y[i] - y[j];
                rz = z[i] - z[j];
                r2 = rx * rx + ry * ry + rz * rz;
                //Truncated distance check
                if (r2 < rc * rc)
                {
                    double force = 24 * epsilon * (2 * pow(sigma2 / r2, 6) - pow(sigma2 / r2, 3)) / r2;
                    fx0[i] += force * rx;
                    fy0[i] += force * ry;
                    fz0[i] += force * rz;
                    fx0[j] -= force * rx;
                    fy0[j] -= force * ry;
                    fz0[j] -= force * rz;
                    PE += 4 * epsilon * (pow(sigma2 / r2, 6) - pow(sigma2 / r2, 3)) - cutoff;
                }
            }
        }
        //Time evolution from 0.01 ~ 300 ps (1 tick = 0.01 ps)
        for (t = 0; t <= t_s; t++)
        {
            if (t % DFT_int == 0)
            {
                printf("now at step %d (%lf%%)\n", t, (double) 100 * t / t_s);
                printf("Time\t\t PE\t\t KE\t\t Total energy\t Average KE\t Temperature\n");
            }
            //Updating every f(t) by the previous f(t+dt) for t > 0
            //Reset f(t) and f(t+dt) for each tick

            for (i = 0; i < N; i++)
            {
                fx0[i] = fx1[i];
                fy0[i] = fy1[i];
                fz0[i] = fz1[i];
                fx1[i] = 0;
                fy1[i] = 0;
                fz1[i] = 0;
            }

            if (t % sum_int == 0)
            {
                //Averaging total sum within the interval
                PE = PE / sum_int;
                KE = KE / sum_int;
                //Find total energy, average KE and temperature
                TE = PE + KE;
                KE_avg = KE / N;
                T_exp = 2 * KE_avg / (3 * kB);
                //Start to accumulate temperature after transient of oscillation
                if (t >= t2)
                {
                    temperature[h] += T_exp;
                }
                //Printing results
                fprintf(stream4, "%lf, %lf, %lf, %lf, %lf, %lf\n", t * dt, PE, KE, TE, KE_avg, T_exp);
                printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", t * dt, PE, KE, TE, KE_avg, T_exp);
                //Reset
                PE = 0;
                KE = 0;
                TE = 0;
                KE_avg = 0;
            }

            if (t % sum_int == 0)
            {                                                                    //Capture one frame per 100 ticks (per ps)
                fprintf(stream, "%i\n", N);                                      //Number of atoms
                fprintf(stream, "Initialized as 300K, time = %lf ps\n", t * dt); //Description
            }
            //Updating each particle position
            for (i = 0; i < N; i++)
            {
                //Verlet method for updating position
                x_temp = x[i];
                x[i] = x[i] + vx[i] * dt + fx0[i] * dt * dt / (2 * mass);
                y[i] = y[i] + vy[i] * dt + fy0[i] * dt * dt / (2 * mass);
                z[i] = z[i] + vz[i] * dt + fz0[i] * dt * dt / (2 * mass);

                //Boundary conditions (Hard wall) (x within -L to 0)
                if (x[i] < -L)
                {
                    x[i] = -L - (x[i] + L);
                    vx[i] = -vx[i];
                    //Turn on oscillator after initial transient
                    if (t >= t1)
                        {
                            vx[i] += v_s * sin(w * t * dt);
                        }
                }

                else if (x[i] > 0)
                {
                    x[i] = -x[i];
                    vx[i] = -vx[i];
                    //Start to accumulate f_wall after transient of oscillation
                    if(t >= t2)
                        {
                            f_wall[h] += mass * fabs(2 * vx[i]) / dt;
                        }
                }

                //Hard wall in other directions
                if (y[i] < 0)
                {
                    y[i] = -y[i];
                    vy[i] = -vy[i];
                }
                else if (y[i] > L)
                {
                    y[i] = L - (y[i] - L);
                    vy[i] = -vy[i];
                }
                if (z[i] < 0)
                {
                    z[i] = -z[i];
                    vz[i] = -vz[i];
                }
                else if (z[i] > L)
                {
                    z[i] = L - (z[i] - L);
                    vz[i] = -vz[i];
                }

                //Start to accumulate total velocity after transient of oscillation
                if (t >= t2)
                {
                    vx_tsum[h] += vx[i];
                }

                if (t % sum_int == 0)
                {                                                          //Capture one frame per 100 ticks (per ps)
                    fprintf(stream, "Ar %lf %lf %lf\n", x[i], y[i], z[i]); //Coordinates
                }
            }
            //Updating each particle force
            for (i = 0; i < (N - 1); i++)
            {
                //Calculation of f(t) by the truncated LJ-potential
                for (j = i + 1; j < N; j++)
                {
                    //Particle doesn't interact with itself
                    rx = x[i] - x[j];
                    ry = y[i] - y[j];
                    rz = z[i] - z[j];
                    r2 = rx * rx + ry * ry + rz * rz;
                    //Truncated distance check
                    if (r2 < rc * rc)
                    {
                        double force = 24 * epsilon * (2 * pow(sigma2 / r2, 6) - pow(sigma2 / r2, 3)) / r2;
                        fx1[i] += force * rx;
                        fy1[i] += force * ry;
                        fz1[i] += force * rz;
                        fx1[j] -= force * rx;
                        fy1[j] -= force * ry;
                        fz1[j] -= force * rz;
                        PE += 4 * epsilon * (pow(sigma2 / r2, 6) - pow(sigma2 / r2, 3)) - cutoff;
                    }
                }
            }

            //Updating each particle velocity
            for (i = 0; i < N; i++)
            {
                //Verlet method for updating velocity
                vx[i] = vx[i] + (fx0[i] + fx1[i]) * dt / (2 * mass);
                vy[i] = vy[i] + (fy0[i] + fy1[i]) * dt / (2 * mass);
                vz[i] = vz[i] + (fz0[i] + fz1[i]) * dt / (2 * mass);
                KE += mass * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]) / 2;
            }
            if (t >= t2 && t % DFT_int == 0)
            {
                //printf("trying to print to stream2\n");
                f_wall[h] = f_wall[h] / DFT_int;
                vx_tsum[h] = vx_tsum[h] / DFT_int;
                temperature[h] = temperature[h] / DFT_int;
                fprintf(stream2, "%lf, %lf, %lf\n", t * dt, f_wall[h], vx_tsum[h]);
                if (t != t2)
                    h++;
            }
        }
        //Printing results

        DFT(f_wall, vx_tsum, temperature);

        printf("File generated.\n");
    }

    fclose(stream);
    fclose(stream2);
    fclose(stream4);

    return 0;
}

void initialize(double x[], double y[], double z[], double vx[], double vy[], double vz[], double fx0[], double fy0[], double fz0[], double fx1[], double fy1[], double fz1[], double f_wall[], double vx_tsum[], double temperature[])
{

    int i, j, k, m = 0;
    //n = 8, L = 100
    double d = L / n;
    double y1, y2, y3, y4, y5, y6;
    //Mass of one Argon molecule, rescaled such that unit is (eV*ps/angstrom)
    double mass = 39.948 * u / (v_c * v_c);
    //Sum of every particles velocity for each component, for center of mass correction
    double vx_s = 0, vy_s = 0, vz_s = 0;

    srand(time(NULL));

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                //Arranging the particles in lattice form, the position is in (angstrom)
                x[m] = -i * d - d / 2;
                y[m] = j * d + d / 2;
                z[m] = k * d + d / 2;
                //Hard coded to prevent generating 0 such that (0 < y <= 1)
                y1 = (double)(rand() + 1) / (RAND_MAX);
                y2 = (double)(rand() + 1) / (RAND_MAX);
                y3 = (double)(rand() + 1) / (RAND_MAX);
                y4 = (double)(rand() + 1) / (RAND_MAX);
                y5 = (double)(rand() + 1) / (RAND_MAX);
                y6 = (double)(rand() + 1) / (RAND_MAX);
                //Applying Muller-Box transform with rescaling
                vx[m] = sqrt(kB * T / mass) * sqrt(-2 * log(y1)) * cos(2 * M_PI * y2);
                vy[m] = sqrt(kB * T / mass) * sqrt(-2 * log(y3)) * cos(2 * M_PI * y4);
                vz[m] = sqrt(kB * T / mass) * sqrt(-2 * log(y5)) * cos(2 * M_PI * y6);
                vx_s += vx[m];
                vy_s += vy[m];
                vz_s += vz[m];
                //Index from 0 ~ 999
                m++;
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        //Center of mass correction
        vx[i] = vx[i] - vx_s / N;
        vy[i] = vy[i] - vy_s / N;
        vz[i] = vz[i] - vz_s / N;
        //Initialization of f(t) and f(t+dt)
        fx0[i] = 0;
        fy0[i] = 0;
        fz0[i] = 0;
        fx1[i] = 0;
        fy1[i] = 0;
        fz1[i] = 0;
    }

    printf("set f_walls to 0\n");
    for (i = 0; i < (t_s - t2) / DFT_int; i++)
    {
        f_wall[i] = 0;
        vx_tsum[i] = 0;
        temperature[i] = 0;
    }
    printf("done initialize\n");
    return;
}

void DFT(double f_wall[], double vx_tsum[], double temperature[])
{
    int i, t;
    //DFT start after transient state
    int DFT_size = (t_s - t2) / DFT_int;
    double mean_temp = 0, mean_vx = 0, mean_f = 0;
    double sin_sum_temp, cos_sum_temp, amplitude_temp;
    double sin_sum_f, cos_sum_f, amplitude_f;
    double sin_sum_vx, cos_sum_vx, amplitude_vx;
    FILE *stream5;

    stream5 = fopen("DFT.csv", "w");

    fprintf(stream5, "Angular frequency, Im (temp), Re (temp), Magnitude (temp), Im (f), Re (f), Magnitude (f), Im (vx), Re (vx), Magnitude (vx)\n");

    //This component is to remove the peak at frequency 0 (which represent the direct current / constant offset, and is meaningless in our case)
    for (i = 0; i < DFT_size; i++)
    {
        mean_temp += temperature[i];
        mean_f += f_wall[i];
        mean_vx += vx_tsum[i];
    }
    mean_temp = mean_temp / DFT_size;
    mean_f = mean_f / DFT_size;
    mean_vx = mean_vx / DFT_size;
    for (i = 0; i < DFT_size; i++)
    {
        temperature[i] -= mean_temp;
        f_wall[i] -= mean_f;
        vx_tsum[i] -= mean_vx;
    }

    for (i = 0; i < DFT_size; i++)
    {
        sin_sum_temp = 0;
        cos_sum_temp = 0;
        amplitude_temp = 0;
        sin_sum_f = 0;
        cos_sum_f = 0;
        amplitude_f = 0;
        sin_sum_vx = 0;
        cos_sum_vx = 0;
        amplitude_vx = 0;
        for (t = 0; t < DFT_size; t++)
        {
            sin_sum_temp += -temperature[t] * sin(i * 2 * M_PI * t / DFT_size);
            cos_sum_temp += temperature[t] * cos(i * 2 * M_PI * t / DFT_size);
            sin_sum_f += -f_wall[t] * sin(i * 2 * M_PI * t / DFT_size);
            cos_sum_f += f_wall[t] * cos(i * 2 * M_PI * t / DFT_size);
            sin_sum_vx += -vx_tsum[t] * sin(i * 2 * M_PI * t / DFT_size);
            cos_sum_vx += vx_tsum[t] * cos(i * 2 * M_PI * t / DFT_size);
        }
        sin_sum_temp /= DFT_size;
        cos_sum_temp /= DFT_size;
        amplitude_temp = sqrt(sin_sum_temp * sin_sum_temp + cos_sum_temp * cos_sum_temp);
        sin_sum_f /= DFT_size;
        cos_sum_f /= DFT_size;
        amplitude_f = sqrt(sin_sum_f * sin_sum_f + cos_sum_f * cos_sum_f);
        sin_sum_vx /= DFT_size;
        cos_sum_vx /= DFT_size;
        amplitude_vx = sqrt(sin_sum_vx * sin_sum_vx + cos_sum_vx * cos_sum_vx);
        fprintf(stream5, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", i * (1 / (DFT_int * dt)) / DFT_size, sin_sum_temp, cos_sum_temp, amplitude_temp, sin_sum_f, cos_sum_f, amplitude_f, sin_sum_vx, cos_sum_vx, amplitude_vx);
    }

    fclose(stream5);
    return;
}
