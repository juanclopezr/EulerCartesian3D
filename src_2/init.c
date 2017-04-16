#include "struct.h"
#include <stdlib.h>
#include <stdio.h>

void init_to_zero(FLOAT *p, int n_points){
    int i;
    for(i=0;i<n_points;i++){
        p[i] = 0.0;
  }
}


physics_grid * create_physics_grid(void){
  physics_grid *G;
  if(!(G = malloc(sizeof(physics_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }

  G->L_x=0.0;
  G->L_y=0.0;
  G->L_z=0.0;
  G->delta_x=0.0;
  G->delta_y=0.0;
  G->delta_z=0.0;
  G->N_x=0;
  G->N_y=0;
  G->N_z=0;
  G->N_cells=0;
  G->P=NULL;
  return G;
}


U_grid * create_U_grid(void){
  U_grid *G;
  if(!(G = malloc(sizeof(U_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->U=NULL;
  return G;
}

F_grid * create_F_grid(void){
  F_grid *G;
  if(!(G = malloc(sizeof(F_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->F_x=NULL;
  G->F_y=NULL;
  G->F_z=NULL;
  return G;
}

U_vector * create_U_vector(void){
    U_vector *U;
    if(!(U = malloc(sizeof(U_vector)))){
      fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
      exit(0);
    }
    U->components = malloc((NDIM + 2)*sizeof(FLOAT));
    init_to_zero(U->components, (NDIM + 2));
    return U;
}

F_vector * create_F_vector(void){
    F_vector *F;
    if(!(F = malloc(sizeof(F_vector)))){
      fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
      exit(0);
    }
    F->components = malloc((NDIM + 2)*sizeof(FLOAT));
    init_to_zero(F->components, (NDIM + 2));
    return F;
}

physics_cell * create_P_vector(void){
    physics_cell *P;
    if(!(P = malloc(sizeof(physics_cell)))){
      fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
      exit(0);
    }
    return P;
}

void init_problem(physics_grid *P, U_grid *U1, U_grid *U2, F_grid *F, int problem){

    P->L_x = 250.0;
    P->L_y = 250.0;
    P->L_z = 250.0;
    P->delta_x = 2.0;
    P->delta_y = 2.0;
    P->delta_z = 2.0;
    P->N_x = (int)(P->L_x/P->delta_x);
    P->N_y = (int)(P->L_y/P->delta_y);
    P->N_z = (int)(P->L_z/P->delta_z);
    P->N_cells = P->N_x * P->N_y * P->N_z;

    U1->N_x = P->N_x;
    U1->N_y = P->N_y;
    U1->N_z = P->N_z;
    U1->N_cells = P->N_cells;

    U2->N_x = P->N_x;
    U2->N_y = P->N_y;
    U2->N_z = P->N_z;
    U2->N_cells = P->N_cells;

    F->N_x = P->N_x;
    F->N_y = P->N_y;
    F->N_z = P->N_z;
    F->N_cells = P->N_cells;

    if(!(P->P = malloc(P->N_cells * sizeof(physics_cell)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
    }

    if(!(U1->U = malloc(U1->N_cells * sizeof(U_vector)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
    }

    if(!(U2->U = malloc(U2->N_cells * sizeof(U_vector)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
    }

    if(!(F->F_x = malloc(F->N_cells * sizeof(F_vector)))){
    fprintf(stderr, "Problem with Fx allocation");
    exit(1);
    }

    if(!(F->F_y = malloc(F->N_cells * sizeof(F_vector)))){
    fprintf(stderr, "Problem with Fy allocation");
    exit(1);
    }

    if(!(F->F_z = malloc(F->N_cells * sizeof(F_vector)))){
    fprintf(stderr, "Problem with Fz allocation");
    exit(1);
    }

    int i;
    for(i = 0; i<P->N_cells; i++)
    {
    //     // U_vector *U_ = malloc(sizeof(U_vector));
    //     // U_vector *U__ = malloc(sizeof(U_vector));
    //     // F_vector *Fx = malloc(sizeof(F_vector)), *Fy = malloc(sizeof(F_vector)), *Fz = malloc(sizeof(F_vector));
    //     // physics_cell *P_ = malloc(sizeof(physics_cell));
    //     // U_vector *U_ = create_U_vector();
    //     // U_vector *U__ = create_U_vector();
    //     // physics_cell *P_ = create_P_vector();
    //     // F_vector *Fx = create_F_vector(), *Fy = create_F_vector(), *Fz = create_F_vector();
        U1->U[i].components = malloc((NDIM+2)*sizeof(FLOAT));
        U2->U[i].components = malloc((NDIM+2)*sizeof(FLOAT));
        F->F_x[i].components = malloc((NDIM+2)*sizeof(FLOAT));
        F->F_y[i].components = malloc((NDIM+2)*sizeof(FLOAT));
        F->F_z[i].components = malloc((NDIM+2)*sizeof(FLOAT));
    //     // U1->U[i] = *U_;
    //     // U2->U[i] = *U__;
    //     // P->P[i] = *P_;
    //     // F->F_x[i] = *Fx;
    //     // F->F_y[i] = *Fy;
    //     // F->F_z[i] = *Fz;
    }

}

FLOAT calculateThermalEnergy(physics_grid *P)
{
    FLOAT volume = (P->delta_x)*(P->delta_y)*(P->delta_z);
    return PRESSURE*volume*3.0/2.0;
}

void destruct_U_vector(U_vector *U)
{
    free(U->components);
    free(U);
}

void destruct_F_vector(F_vector *F)
{
    free(F->components);
    free(F);
}

void destruct_grids(physics_grid *P, U_grid *U1, U_grid *U2, F_grid *F)
{
    int i;
    for(i=0; i<P->N_cells; i++)
    {
        if(i!=P->N_cells)
        {
            free(U1->U[i].components);
            free(U2->U[i].components);
            free(F->F_x[i].components);
            free(F->F_y[i].components);
            free(F->F_z[i].components);
        }

    }
    free(U1->U);
    free(U2->U);
    free(F->F_x);
    free(F->F_y);
    free(F->F_z);
    free(U1);
    free(U2);
    free(F);
    free(P->P);
    free(P);
}
