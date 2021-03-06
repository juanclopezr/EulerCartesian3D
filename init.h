void init_problem(physics_grid *P, U_grid *U1, U_grid *U2, F_grid *F, int problem);
physics_grid * create_physics_grid(void);
U_grid * create_U_grid(void);
F_grid * create_F_grid(void);
U_vector * create_U_vector(void);
F_vector * create_F_vector(void);
physics_cell * create_P_vector(void);
FLOAT calculateThermalEnergy(physics_grid *P);
void destruct_U_vector(U_vector *U);
void destruct_F_vector(F_vector *F);
void destruct_grids(physics_grid *P, U_grid *U1, U_grid *U2, F_grid *F);
