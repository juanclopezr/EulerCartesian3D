int transform_U(U_grid *U, int pos_x, int pos_y, int pos_z, int prop);
//int index(U_grid *U, int pos_x, int pos_y, int pos_z)
FLOAT extract_rho(U_grid *U, int pos_x, int pos_y, int pos_z);
FLOAT extract_u(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho);
FLOAT extract_v(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho);
FLOAT extract_w(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho);
FLOAT extract_E(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho);
FLOAT calce(FLOAT E, FLOAT u, FLOAT v, FLOAT w);
FLOAT calcp(FLOAT rho, FLOAT e);
FLOAT calch(FLOAT E, FLOAT p, FLOAT rho);
FLOAT calcs(FLOAT h);
FLOAT calcsps(U_grid *U);
FLOAT calcdt(physics_grid *P, FLOAT sps);
int transform_F(F_grid *F, int pos_x, int pos_y, int pos_z, int pos_g, int prop);
void calcF(F_grid *F, U_grid *U, int pos_x, int pos_y, int pos_z);
void newU(U_grid *U, F_grid *F, physics_grid *P, int pos_x, int pos_y, int pos_z, FLOAT dt);
