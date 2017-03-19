#include <math.h>
#include <stdio.h>
#include <stdlib.h>


((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z + i
//Extracts the value of rho for a given cell
FLOAT extract_rho(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z];
}

//Extracts the value of u for a given cell
FLOAT extract_u(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z + 1]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of v for a given cell
FLOAT extract_v(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z + 2]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of w for a given cell
FLOAT extract_w(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z + 3]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of E for a given cell
FLOAT extract_E(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z + 4]/extract_rho(U, pos_x, pos_y, pos_z);
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
FLOAT calch(FLOAT E, FLOAT p, FLOAT rho)
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

//Finds the maximum sounds speed of the whole array
FLOAT calcsps(U_grid *U, FLOAT cs)
{
	return 0;
}

//Calculates dt
FLOAT calcdt(FLOAT dx, FLOAT sps)
{
	return 0.5*dx/sps;
}

/*TODO
 * Write the function to calculate the F-vectors.
 * Use the fourth order interpolation scheme to calculate 'F-half' with neighboring cell values.
 * Define the iteration function which will be used in the main.*/
 
//Calculate the F-vectors for the corresponding U-vectors

void calcF(F_grid *F, U_grid *U)
{
	int i;
	int j;
	int k;
	FLOAT rho,u,v,w,E,e,p;
	//FLOAT h,cs,sps,dt; To be used in iterations
	for(i=0;i<U->N_x;i++)
	{
		for(j=0;j<U->N_y;j++)
		{
			for(k=0;k<U->N_y;k++)
			{
				rho = extract_rho(U,i,j,k);
				u = extract_u(U,i,j,k);
				v = extract_v(U,i,j,k);
				w = extract_w(U,i,j,k);
				E = extract_E(U,i,j,k);
				e = calce(E,u,v,w);
				p = calcp(rho,e);
				//h = calch(E,p,rho);
				//cs = calcs(h);
				//sps = calcsps(U,cs);
				//dt = calcdt(dx,sps); To be used in iterations
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z)*NDIM] = rho*u;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 1)*NDIM = rho*v;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 2)*NDIM = rho*w;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z)*NDIM + 1 = rho*u*u+p;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 1)*NDIM + 1 = rho*u*v;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 2)*NDIM + 1 = rho*u*w;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z)*NDIM + 2 = rho*u*v;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 1)*NDIM + 2 = rho*v*v+p;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 2)*NDIM + 2 = rho*v*w;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z)*NDIM + 3 = rho*u*w;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 1)*NDIM + 3 = rho*v*w;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 2)*NDIM + 3 = rho*w*w+p;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z)*NDIM + 4 = rho*(E+p/rho)*u;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 1)*NDIM + 4 = rho*(E+p/rho)*v;
				F[(((i*U->N_x + j)*U->N_y + k)*U->N_z + 2)*NDIM + 4 = rho*(E+p/rho)*w;
			}
		}
	}
}
