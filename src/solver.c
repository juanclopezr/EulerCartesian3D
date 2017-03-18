#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Extracts the value of rho for a given cell
FLOAT extract_rho(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z];
}

FLOAT eaxtract_u(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 1]/extract_rho(U, pos_x, pos_y, pos_z);
}

FLOAT eaxtract_v(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 2]/extract_rho(U, pos_x, pos_y, pos_z);
}

FLOAT eaxtract_w(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 3]/extract_rho(U, pos_x, pos_y, pos_z);
}

FLOAT eaxtract_E(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 4]/extract_rho(U, pos_x, pos_y, pos_z);
}
