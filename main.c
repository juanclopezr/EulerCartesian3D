#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "init.h"
#include "io.h"
#include "solver.h"
#include "problem.h"
#include "tube.h"

int main(int argc, char **argv)
{
    if(argc)
    {
        if(atoi(argv[1]) == SOD)
        {
            solve_SOD();
        }
        else if(atoi(argv[1]) == SEDOV)
        {
            physics_grid * P_state;
            U_grid * U_state, * U_temp;
            F_grid  * F_state;

            P_state = create_physics_grid();
            U_state = create_U_grid();
            U_temp = create_U_grid();
            F_state = create_F_grid();

            init_problem(P_state, U_state, U_temp, F_state, SEDOV);

            // constants used by transform3d function
            Nx = P_state->N_x;
            Ny = P_state->N_y;
            Nz = P_state->N_z;

            init_cond(P_state, U_state, F_state);
            solve_SEDOV(P_state, U_state, U_temp, F_state);
            destruct_grids(P_state, U_state, U_temp, F_state);
        }
    }
    else
    {
        printf("args missing\n");
    }
    return 0;
}
