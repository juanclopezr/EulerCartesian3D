#include "struct.h"
#include "solver.h"
#include "init.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Gives the initial condition for the problem
void init_cond(physics_grid *P, U_grid *U, F_grid *F)
{
    int i, j, k, index, middle;
    double big_E;
    physics_cell P_;
    U_vector U_;
    F_vector F_x, F_y, F_z;
    big_E = calculateThermalEnergy(P);
    int halfx, halfy, halfz;

    halfx = halfy = halfz = 0.5*P->N_x;
    middle = transform3d(halfx, halfy, halfz);
    for(i = 0; i<P->N_x; i++)
    {
        for(j = 0; j<P->N_y; j++)
        {
            for(k = 0; k<P->N_z; k++)
            {
                index = transform3d(i, j, k);
                P_ = P->P[index];
                U_ = U->U[index];
                F_x = F->F_x[index];
                F_y = F->F_y[index];
                F_z = F->F_z[index];
                P_.rho = RHO;
                P_.p = PRESSURE;
                P_.u = VX;
                P_.v = VY;
                P_.w = VZ;
                P_.E = big_E;
                if(index == middle)
                {
                    P_.E = EXP_ENERGY;
                    P_.p = calcp(P_.rho, P_.E);
                }
                // if ((i == halfx) || (i == halfx-1) || (i == halfx+1))
                // {
                //     if ((j == halfx) || (j == halfx-1) || (j == halfx+1))
                //     {
                //         if ((k == halfx) || (k == halfx-1) || (k == halfx+1))
                //         {
                //             P_.E = EXP_ENERGY;
                //             P_.p = calcp(P_.rho, P_.E);
                //         }
                //     }
                // }
                fromP_calcU(P_, &U_);
                fromU_calcF(U_, &F_x, &F_y, &F_z);
                P->P[index] = P_;
                U->U[index] = U_;
                F->F_x[index] = F_x;
                F->F_y[index] = F_y;
                F->F_z[index] = F_z;
            }
        }
    }
}

void solve_SEDOV(physics_grid *P, U_grid *U, U_grid *U_temp, F_grid *F)
{
    FLOAT T, dt, speed;
    FLOAT *density;
    T = 0;
    dt = 1e-7;
    int i = 0;
    int pos_10, pos_60, pos_120;
    int ver_10 = 1, ver_60 = 1, ver_120 = 1;
    FLOAT last_10, last_60, last_120;
    pos_10 = transform3d(10.0/P->delta_x + 0.5*P->N_x, 10.0/P->delta_y + 0.5*P->N_y, 10.0/P->delta_z + 0.5*P->N_z);
    pos_60 = transform3d(60.0/P->delta_x + 0.5*P->N_x, 60.0/P->delta_y + 0.5*P->N_y, 60.0/P->delta_z + 0.5*P->N_z);
    pos_120 = transform3d(120.0/P->delta_x + 0.5*P->N_x, 120.0/P->delta_y + 0.5*P->N_y, 120.0/P->delta_z + 0.5*P->N_z);

    last_10 = U->U[pos_10].components[0];
    last_60 = U->U[pos_60].components[0];
    last_120 = U->U[pos_120].components[0];
    while (1)
    {
        speed = calculateNextU(P, U, *U_temp, dt);
        if((speed > 0) && (speed < 3e8))
        // if(speed > 0)
        {
            dt = 0.5*P->delta_x/speed;
        }
        T += dt;
        if((last_10 != U->U[pos_10].components[0]) && ver_10)
        {
            updatePhysics(P, *U);
            print_rho(0, P);
            ver_10 = 0;
        }
        if((last_60 != U->U[pos_60].components[0]) && ver_60)
        {
            updatePhysics(P, *U);
            print_rho(1, P);
            ver_60 = 0;
        }
        if((last_120 != U->U[pos_120].components[0]) && ver_120)
        {
            updatePhysics(P, *U);
            print_rho(2, P);
            ver_120 = 0;
            break;
        }
        last_10 = U->U[pos_10].components[0];
        last_60 = U->U[pos_60].components[0];
        last_120 = U->U[pos_120].components[0];

        printf("%f\n", T);
        i += 1;
    }
}
