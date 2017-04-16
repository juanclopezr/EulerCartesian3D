#include "struct.h"
#include "solver.h"
#include "init.h"
#include <stdio.h>
#include <stdlib.h>

//Gives the initial condition for the problem
void init_cond(physics_grid *P, U_grid *U, F_grid *F)
{
    int i, j, k, index, middle;
    double big_E;
    physics_cell P_;
    U_vector U_;
    F_vector F_x, F_y, F_z;
    big_E = calculateThermalEnergy(P);

    middle = transform3d(0.5*P->N_x, 0.5*P->N_y, 0.5*P->N_z);
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
                P_.u = P_.v = P_.w = 0.0;
                P_.E = big_E;
                if (middle == index)
                {
                    P_.E = EXP_ENERGY;
                }
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
    T = 0;
    dt = 1e-6;
    int i = 0;
    while (i < 10)
    {
        speed = calculateNextU(P, U_temp, U, dt);
        dt = P->delta_x/speed;
        T += dt;
        printf("%f %f\n", dt, speed);
        i += 1;
    }
    updatePhysics(P, U);
}
