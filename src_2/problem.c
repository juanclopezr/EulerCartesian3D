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
                P_.u = VX;
                P_.v = VY;
                P_.w = VZ;
                P_.E = big_E;
                if (middle == index)
                {
                    P_.E = EXP_ENERGY;
                    P_.p = calcp(P_.rho, P_.E);
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
    dt = 1e-7;
    int i = 0;
    printf("%.12f\n",VX);
    while (i < 200)
    {
        speed = calculateNextU(P, U, *U_temp, dt);
        if(speed > 0)
        {
            dt = 0.5*P->delta_x/speed;
        }
        T += dt;

        printf("%f\n", speed);

        updatePhysics(P, *U);
        print_E(i, P);
        i += 1;
    }

}
