#include <stdio.h>
#include "struct.h"
#include "init.h"
#include "io.h"
#include "solver.h"
#include "problem.h"

int main(int argc, char **argv){
    physics_grid * P_state;
    U_grid * U_state;
    F_grid  * F_state;

    P_state = create_physics_grid();
    U_state = create_U_grid();
    F_state = create_F_grid();

    init_problem(P_state, U_state, F_state, SEDOV);
    init_cond(U_state);
    prob_solve(U_state, F_state, P_state, 0.1);
    //print_L(P_state);

    return 0;
}
