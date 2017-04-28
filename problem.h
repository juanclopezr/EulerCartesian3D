// void prob_solve(U_grid *U, F_grid *F, physics_grid *P, FLOAT T);
void init_cond(physics_grid *P, U_grid *U, F_grid *F);
void solve_SEDOV(physics_grid *P, U_grid *U, U_grid *U_temp, F_grid *F);
