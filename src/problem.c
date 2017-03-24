#include "struct.h"
#include "solver.h"
#include <stdio.h>

//Gives the initial condition for the problem
void init_cond(U_grid *U)
{
    int i, j, k;
    for(i=0; i<U->N_x; i++)
    {
        for(j=0; j<U->N_y; j++)
        {
            for(k=0; k<U->N_z; k++)
            {
                U->U[transform_U(U, i, j, k, 0)] = RHO;
            }
        }
    }
	U->U[transform_U(U,(U->N_x-1)/2,(U->N_y-1)/2,(U->N_z-1)/2, 4)] = RHO*1e6;
}

//Run the whole simulation
void prob_solve(U_grid *U, F_grid *F, physics_grid *P, FLOAT T)
{
	FLOAT t = 0;
	FLOAT dt;
	int i,j,k;
	while(t<T)
	{
		dt = calcdt(P,calcsps(U));
		for(i=0;i<F->N_x;i++)
		{
			for(j=0;j<F->N_y;j++)
			{
				for(k=0;k<F->N_z;k++)
				{
					calcF(F,U,i,j,k);
				}
			}
		}
		
		for(i=0;i<U->N_x;i++)
		{
			for(j=0;j<U->N_y;j++)
			{
				for(k=0;k<U->N_z;k++)
				{
					newU(U,F,P,i,j,k,dt);
				}
			}
		}
		printf("%f %f\n", t, dt);
		t += dt;
	}
}
