int transform3d(int pos_x, int pos_y, int pos_z);

FLOAT calce(FLOAT E, FLOAT u, FLOAT v, FLOAT w);
FLOAT calcp(FLOAT rho, FLOAT e);
FLOAT calch(FLOAT E, FLOAT p, FLOAT rho);
FLOAT calcs(FLOAT h);

FLOAT absolute(FLOAT value);
FLOAT max_speed(U_vector *U);

void fromU_calcP(U_vector U, physics_cell *P);
void fromP_calcU(physics_cell P, U_vector *U);
void fromU_calcFx(U_vector U, F_vector *F);
void fromU_calcFy(U_vector U, F_vector *F);
void fromU_calcFz(U_vector U, F_vector *F);
void fromU_calcF(U_vector U, F_vector *Fx, F_vector *Fy, F_vector *Fz);

F_vector *calculateF_half(U_grid *U, int i, int j, int k, int i_, int j_, int k_);
U_vector *calculateU_half(U_vector *U1, U_vector *U2);
F_vector *calculateF_diff(F_vector *F_plus, F_vector *F_minus);

FLOAT calculateNextU(physics_grid *P, U_grid *U_new, U_grid *U_past, double dt);
void updatePhysics(physics_grid *P, U_grid *U);
