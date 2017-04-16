#include "struct.h"
#include <stdio.h>
#include "solver.h"

void print_L(physics_grid *G){
  fprintf(stdout, "L_x = %f\n", G->L_x);
  fprintf(stdout, "L_y = %f\n", G->L_y);
  fprintf(stdout, "L_z = %f\n", G->L_z);
}

void print_E(physics_grid *P)
{
    FILE *energy = fopen("energy.dat", "w");
    int i, j, k, index;
    for(i = 0; i<P->N_x; i++)
    {
        for(j = 0; j<P->N_y; j++)
        {
            for(k = 0; k<P->N_z; k++)
            {
                index = transform3d(i, j, k);
                fprintf(energy, "%f\n", P->P[index].E);
            }
        }
    }
    fclose(energy);
}
