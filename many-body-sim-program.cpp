#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
using namespace std;

struct Body{
    Body(){
        vel[0] = 0;
        vel[1] = 0;
        vel[2] = 0;
    }
    double pos[3];
    double vel[3];
};
double dot(double A[],double B[]) {
    double C = 0.0;
    for (int i = 0; i < 3; i++)
    {
        C += A[i]*B[i];
    }
    return C;
}

double magnitude(double s[]) {
    return sqrt(dot(s,s));
}

void collision(Body &a, Body &b) {
    double dij[3] = {
        a.pos[0] - b.pos[0],
        a.pos[1] - b.pos[1],
        a.pos[2] - b.pos[2]
    };
    double dji[3] = {
        b.pos[0] - a.pos[0],
        b.pos[1] - a.pos[1],
        b.pos[2] - a.pos[2]};
    double a_p[3], b_p[3], a_n[3], b_n[3];
    for(int i = 0 ; i < 3; ++i){
        a_p[i] = dot(a.vel, dji) * dji[i] / magnitude(dji);
        b_p[i] = dot(b.vel, dij) * dij[i] / magnitude(dij);
        a_n[i] = a.vel[i] - a_p[i];
        b_n[i] = b.vel[i] - b_p[i];
    }
    for (int i = 0; i < 3; i++)
    {
        a.vel[i] = a_n[i] + b_p[i];
        b.vel[i] = b_n[i] + a_p[i];
    }
}

void calculate_force(int n, Body particle[], double dt, int pi, int num_threads, double force[])
{
    double t[num_threads][3];
#pragma omp parallel
    {
        for (int i = 0; i < 3; i++)
        {
            t[omp_get_thread_num()][i] = 0.0;
        }
#pragma omp for schedule(static)
        for (int pj = 0; pj < n; pj++)
        {
            if (pj != pi)
            {
                const double dij[3] = {
                    particle[pj].pos[0] - particle[pi].pos[0],
                    particle[pj].pos[1] - particle[pi].pos[1],
                    particle[pj].pos[2] - particle[pi].pos[2]};

                double dist2 = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
                if (dist2 <= 1)
                {
                    collision(particle[pi], particle[pj]);
                    dist2 = 1.0;
                }
                double idist2 = 1.0 / dist2;
                idist2 *= (1.0 / sqrt(dist2));
                // F = C * m * m / ||x2 - x1||^2 * (x2 - x1) / ||x2 - x1||
                t[omp_get_thread_num()][0] += idist2 * dij[0];
                t[omp_get_thread_num()][1] += idist2 * dij[1];
                t[omp_get_thread_num()][2] += idist2 * dij[2];
            }
        }
#pragma omp for
        for (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < num_threads; i++)
            {
                force[j] += t[i][j];
            }
        }
    }
}

void simulate(int n, Body particle[], double dt, int num_threads){
    double boundary[3] = {200.0, 100.0, 400.0};

    for (int pi = 0; pi < n; pi++)
    {

        double force[3] = {0.0, 0.0, 0.0};

        // Calculate total force
        calculate_force(n, particle, dt, pi, num_threads, force);
        for (int i = 0; i < 3; i++)
        {
            particle[pi].vel[i] += force[i] * 0.5 * dt;
            particle[pi].pos[i] += particle[pi].vel[i] * dt;
            particle[pi].vel[i] += force[i] * 0.5 * dt;
        }
        for (int i = 0; i < 3; i++)
        {
            if (particle[pi].pos[i] > boundary[i]){
                particle[pi].pos[i] = 2 * boundary[i] - particle[pi].pos[i];
                particle[pi].vel[i] *= -1;
            }
            else if (particle[pi].pos[i] < 0){
                particle[pi].pos[i] *= -1;
                particle[pi].vel[i] *= -1;
            }
        }
    }
}
int main()
{
    ifstream infile;
    infile.open("inp.txt");
    double l, w, d, r, m, dt;
    int n;
    if (!infile.is_open()){
        printf("Failed to open input file\n");
        exit(1);
    }
    infile >> l;
    infile >> w;
    infile >> d;
    infile >> n;
    infile >> r;
    infile >> m;
    infile >> dt;
    Body bodies[n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            infile >> bodies[i].pos[j];
        }
    }
    infile.close();
    int steps;
    cout << "Number of steps : ";
    cin >> steps;
    int num_threads;
    cout << "Number of threads : ";
    cin >> num_threads;
    omp_set_num_threads(num_threads);
    ofstream outfile;
    outfile.open("test_run.txt");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            outfile << bodies[i].pos[j] << " ";
        }
        outfile << "\n";
    }
    double total_time = 0.0;
    for (int step = 0; step < steps; step++)
    {
        double start = omp_get_wtime();
        simulate(n, bodies, dt, num_threads);
        double end = omp_get_wtime();
        if ((step + 1) % 100 == 0)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    outfile << bodies[i].pos[j] << " ";
                }
                outfile << "\n";
            }
            outfile << "\n";
        }
        total_time += (end - start);
        if(step%1000 == 0){
            cout << "step " << step << " completed\n";
        }
    }
    cout << "Time taken (parallel part) for " << num_threads << " threads : " << total_time << endl;
    cout << "Average Time taken for each step" << num_threads << " threads : " << (total_time/steps) << endl;
    return 0;
}
