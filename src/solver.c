#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Extracts the value of rho for a given cell
FLOAT extract_rho(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z];
}

//Extracts the value of u for a given cell
FLOAT eaxtract_u(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 1]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of v for a given cell
FLOAT eaxtract_v(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 2]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of w for a given cell
FLOAT eaxtract_w(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 3]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of E for a given cell
FLOAT eaxtract_E(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[(pos_x*U->M_x + pos_y)*U->N_y + pos_z*U->N_z + 4]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Calculates the value of e
FLOAT calce(FLOAT E, FLOAT u, FLOAT v, FLOAT w)
{
	return E - 0.5*(u*u+v*v+w*w);
}

//Calculates the value of p
FLOAT calcp(FLOAT rho, FLOAT e)
{
	return rho*e*(GAMMA-1);
}

//Calculates h
FLOAT calh(FLOAT E, FLOAT p, FLOAT rho)
{
	return E + p/rho;
}

//Calculates cs
FLOAT calcs(FLOAT h)
{
	return sqrt(h*(GAMMA+1))
}

/*TODO
 * Find the general max value of (u+cs,v+cs,w+cs) of all cells
 * This is necessary for calculating the maximum sound speed*/

//Calculates dt
FLOAT dt(FLOAT dx, FLOAT sps)
{
	return 0.5*dx/sps;
}

/*TODO
 * Write the function to calculate the F-vectors.
 * Use the fourth order interpolation scheme to calculate 'F-half' with neighboring cell values.
 * Define the iteration function which will be used in the main.*/
 
