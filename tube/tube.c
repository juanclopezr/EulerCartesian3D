#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tube.h"

int main(int argc, char **argv)
{
    int i;
    FLOAT t=0, dt;
    density = malloc(N*sizeof(FLOAT));
    speed = malloc(N*sizeof(FLOAT));
    pressure = malloc(N*sizeof(FLOAT));
    energy = malloc(N*sizeof(FLOAT));
    
    FLOAT *t_density, *t_speed, *t_pressure, *t_energy;
    t_density = malloc(N*sizeof(FLOAT));
    t_speed = malloc(N*sizeof(FLOAT));
    t_pressure = malloc(N*sizeof(FLOAT));
    t_energy = malloc(N*sizeof(FLOAT));
    x = linspace(0, 1.0, N);
    
    init_tube(energy, density, speed, pressure, x);
    
    FLOAT dx, max_speed, abs_speed, check = 0, sound_speed;
    dx = x[1] - x[0];
    int check_at = 0.9/dx, low, up;
    
    dt = 1e-6;
    max_speed = 0;    
    low = check_at - 1;
    up = check_at + 1;
    
    //if animation is wanted
    //char name[100];
    //int j = 0;
    
    while(check < 0.1)
    {
        //sprintf(name, "%d_data.dat", j); // if animation is wanted
        for(i = 1; i<N-1; i++)
        {
            lax(i, dt, dx, t_density, t_speed, t_energy, t_pressure);
        }
        max_speed = 0;
        for(i = 1; i<N-1; i++)
        {
            density[i] = t_density[i];
            speed[i] = t_speed[i];
            energy[i] = t_energy[i];
            pressure[i] = t_pressure[i];
            sound_speed = sqrt(GAMMA*pressure[i]/density[i]);
            if(sound_speed + speed[i] > max_speed)
            {
                max_speed = sound_speed + speed[i];
            }
        }
        dt = dx/max_speed;
        t += dt;
        //print_tube(name, t); // if animation is wanted
        check = pressure[low] - pressure[up];
        
        //j += 1;
    }
    print_tube("tube.dat", t);
    free(x); 
    free(energy); free(density); free(speed); free(pressure);
    free(t_energy); free(t_density); free(t_speed); free(t_pressure);
	return 0;
}

FLOAT calculate_pressure(FLOAT energy, FLOAT density, FLOAT speed)
{
    return (GAMMA - 1) * (energy - 0.5*density*speed*speed);
}

FLOAT calculate_energy(FLOAT pressure, FLOAT density, FLOAT speed)
{
    return pressure/(GAMMA - 1) + 0.5*density*speed*speed;
}

FLOAT *linspace(FLOAT init, FLOAT last, FLOAT steps)
{
    int i;
    FLOAT dx = (last-init)/(steps-1);
    FLOAT *array = malloc(steps*sizeof(FLOAT));
    
    array[0] = init;
    for(i=1; i<steps; i++)
    {
        array[i] = array[i-1]+dx;
    }
    
    return array;
}

void init_tube(FLOAT *energy, FLOAT *density, FLOAT *speed, FLOAT *pressure, FLOAT *space)
{
    int i;
    FLOAT x;
    for(i=0; i<N; i++)
    {
        x = space[i];
        if(x <= X0)
        {
            density[i] = RHO1;
            pressure[i] = P1;
        }
        else
        {
            density[i] = RHO5;
            pressure[i] = P5;
        }
        speed[i] = U1;
        energy[i] = calculate_energy(pressure[i], density[i], speed[i]);
    }
}

FLOAT *calc_F(FLOAT *u)
{
    FLOAT rho, v, E, p;  
    rho = u[0];
    v = u[1]/rho;
    E = u[2];
    p = calculate_pressure(E, rho, v);
   
    FLOAT *vector = malloc(3*sizeof(FLOAT));
    vector[0] = rho*v;
    vector[1] = vector[0]*v + p;
    vector[2] = v*(E + p);
    return vector;
}

FLOAT *calc_U(FLOAT rho, FLOAT v, FLOAT E)
{
    FLOAT *vector = malloc(3*sizeof(FLOAT));
    vector[0] = rho;
    vector[1] = rho*v;
    vector[2] = E;
    return vector;
}

void from_vector_to_relevant(FLOAT *vector)
{
    vector[1] *= 1.0/vector[0];
}

void lax(int pos, FLOAT dt, FLOAT dx, FLOAT *t_density, FLOAT *t_speed, FLOAT *t_energy, FLOAT *t_pressure)
{
    int i;
    FLOAT ratio = dt/dx;
    FLOAT *U, *U_before, *half1, *half2, *U_next;
    FLOAT *F, *F_before, *F_half1, *F_half2, *F_next;
    
    half1 = malloc(3*sizeof(FLOAT));
    half2 = malloc(3*sizeof(FLOAT));

    U = calc_U(density[pos], speed[pos], energy[pos]);
    U_next = calc_U(density[pos+1], speed[pos+1], energy[pos+1]);
    U_before = calc_U(density[pos-1], speed[pos-1], energy[pos-1]);

    F = calc_F(U);
    F_next = calc_F(U_next);
    F_before = calc_F(U_before);

    for(i=0; i<3; i++)
    {
        half1[i] = 0.5*(U[i] + U_next[i]) - 0.5*ratio*(F_next[i] - F[i]);
        half2[i] = 0.5*(U_before[i] + U[i]) - 0.5*ratio*(F[i] - F_before[i]);
    }
    
    F_half1 = calc_F(half1);
    F_half2 = calc_F(half2);
    for(i=0; i<3; i++)
    {
        U[i] = U[i] - ratio*(F_half1[i] - F_half2[i]);
    }
    
    from_vector_to_relevant(U);
    
    t_density[pos] = U[0];
    t_speed[pos] = U[1];
    t_energy[pos] = U[2];
    t_pressure[pos] = calculate_pressure(U[2], U[0], U[1]);
    
    free(U); free(U_next); free(U_before);
    free(F); free(F_next); free(F_before);
    free(half1); free(half2); free(F_half1);
    free(F_half2);
}

void print_tube(const char *name, double t)
{
    int i;
    FILE *output = fopen(name, "w");
    fprintf(output, "gamma = %f\n", GAMMA);
    fprintf(output, "rho1 = %f\n", RHO1);
    fprintf(output, "rho5 = %f\n", RHO5);
    fprintf(output, "u1 = %f\n", U1);
    fprintf(output, "u5 = %f\n", U1);
    fprintf(output, "p1 = %f\n", P1);
    fprintf(output, "p5 = %f\n", P5);
    fprintf(output, "t = %f\n", t);
    fprintf(output, "x0 = %f\n", X0);
    for(i=0; i<N; i++)
    {
        fprintf(output, "%f %f %f %f %f\n", x[i], energy[i], density[i], speed[i], pressure[i]);
    }
    fclose(output);
}
