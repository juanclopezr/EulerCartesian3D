#ifndef STRUCT_H
#define STRUCT_H

#define GAMMA 1.4

#define SEDOV 1
#define NDIM 3

#define EXP_ENERGY 1e13 // J/kg
#define PRESSURE 101325.0 // Pa
#define TEMPERATURE 300 // Kelvin
#define RHO 1.0 // kg/m3
#define VX 0.0 // m/s
#define VY 0.0 // m/s
#define VZ 0.0 // m/s

#define FLOAT double

int Nx, Ny, Nz;

typedef struct physics_grid_str
{
    FLOAT L_x;
    FLOAT L_y;
    FLOAT L_z;
    FLOAT delta_x;
    FLOAT delta_y;
    FLOAT delta_z;
    int N_x;
    int N_y;
    int N_z;
    int N_cells;
    struct physics_str *P;
} physics_grid;

typedef struct U_grid_str{
    int N_x;
    int N_y;
    int N_z;
    int N_cells;
    struct U_vector_str *U;
} U_grid;

typedef struct F_grid_str{
    int N_x;
    int N_y;
    int N_z;
    int N_cells;
    struct F_vector_str *F_x;
    struct F_vector_str *F_y;
    struct F_vector_str *F_z;
} F_grid;

typedef struct physics_str{
    FLOAT rho;
    FLOAT p;
    FLOAT u;
    FLOAT w;
    FLOAT v;
    FLOAT E;
} physics_cell;

typedef struct U_vector_str{
    FLOAT *components;
    // FLOAT *components[NDIM+2];
} U_vector;

typedef struct F_vector_str{
    FLOAT *components;
    // FLOAT *components[NDIM+2];
} F_vector;
#endif
