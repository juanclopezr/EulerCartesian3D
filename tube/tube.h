#define FLOAT double
#define N 500
#define GAMMA 1.4
#define RHO1 1.0
#define RHO5 0.125
#define P1 1.0
#define P5 0.1
#define U1 0.0
#define X0 0.5

FLOAT *linspace(FLOAT init, FLOAT last, FLOAT steps);
FLOAT calculate_pressure(FLOAT energy, FLOAT density, FLOAT speed);
FLOAT calculate_energy(FLOAT pressure, FLOAT density, FLOAT speed);
void init_tube(FLOAT *energy, FLOAT *density, FLOAT *speed, FLOAT *pressure, FLOAT *space);
void lax(int pos, FLOAT dt, FLOAT dx, FLOAT *t_density, FLOAT *t_speed, FLOAT *t_energy, FLOAT *t_pressure);
void print_tube(const char *name, double t);

FLOAT *calc_F(FLOAT *U);
FLOAT *calc_U(FLOAT rho, FLOAT v, FLOAT e);
FLOAT *density, *speed, *pressure, *energy, *x;
