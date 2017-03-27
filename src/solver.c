#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"


//Transform function for accessing the U-indices
int transform_U(U_grid *U, int pos_x, int pos_y, int pos_z, int prop)
{
	return ((prop*U->N_z+pos_z)*U->N_y+pos_y)*U->N_x+pos_x;
}

//Extracts the value of rho for a given cell
FLOAT extract_rho(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U->U[transform_U(U,pos_x,pos_y,pos_z, 0)];
}

//Extracts the value of u for a given cell
FLOAT extract_u(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho)
{
	return U->U[transform_U(U,pos_x,pos_y,pos_z, 1)]/rho;
}

//Extracts the value of v for a given cell
FLOAT extract_v(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho)
{
	return U->U[transform_U(U,pos_x,pos_y,pos_z, 2)]/rho;
}

//Extracts the value of w for a given cell
FLOAT extract_w(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho)
{
	return U->U[transform_U(U,pos_x,pos_y,pos_z, 3)]/rho;
}

//Extracts the value of E for a given cell
FLOAT extract_E(U_grid *U, int pos_x, int pos_y, int pos_z, FLOAT rho)
{
	return U->U[transform_U(U,pos_x,pos_y,pos_z, 4)]/rho;
}

//Calculates the value of e
FLOAT calce(FLOAT E, FLOAT u, FLOAT v, FLOAT w)
{
	return E - 0.5*(pow(u, 2) + pow(v, 2) + pow(w, 2));
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
	return sqrt(h*(GAMMA+1));
}

//Finds the maximum sound speed of the whole array
FLOAT calcsps(U_grid *U)
{
	FLOAT sps_max = 0;
	FLOAT temp, E, u, v, w, e, p, rho;
	int i, j, k;
	for(i=0; i<U->N_x; i++)
	{
		for(j=0; j<U->N_y; j++)
		{
			for(k=0; k<U->N_z; k++)
			{
				rho = extract_rho(U,i,j,k);
				E = extract_E(U,i,j,k,rho);
				u = extract_u(U,i,j,k,rho);
				v = extract_v(U,i,j,k,rho);
				w = extract_w(U,i,j,k,rho);
				e = calce(E,u,v,w);
				p = calcp(rho,e);
                
                //printf("%f %f %f %f %f %f %f\n", rho, E, u, v, w, e, p);
				temp = 0;
				if(u>temp)
				{
					temp = u;
				}
				if(v>temp)
				{
					temp = v;
				}
				if(w>temp)
				{
					temp = w;
				}
				temp += calcs(calch(E,p,rho));
				if(temp>sps_max)
				{
					sps_max=temp;
				}
			}
		}
	}
    printf("%f\n", sps_max);
	return sps_max;
}

//Calculates dt
FLOAT calcdt(physics_grid *P, FLOAT sps)
{
	return 0.5*P->delta_x/sps;
}

//Transform function for accessing the F-indices
int transform_F(F_grid *F, int pos_x, int pos_y, int pos_z, int pos_g, int prop)
{
	return (((prop*NDIM+pos_g)*F->N_z+pos_z)*F->N_y+pos_y)*F->N_x+pos_x;
}

//Calculate F-vector
void calcF(F_grid *F, U_grid *U, int pos_x, int pos_y, int pos_z)
{
	FLOAT rho,u,v,w,E,e,p;
	rho = extract_rho(U,pos_x,pos_y,pos_z);
	u = extract_u(U,pos_x,pos_y,pos_z,rho);
	v = extract_v(U,pos_x,pos_y,pos_z,rho);
	w = extract_w(U,pos_x,pos_y,pos_z,rho);
	E = extract_E(U,pos_x,pos_y,pos_z,rho);
	e = calce(E,u,v,w);
	p = calcp(rho,e);
	F->F[transform_F(F,pos_x,pos_y,pos_z,0,0)] = rho*u;
	F->F[transform_F(F,pos_x,pos_y,pos_z,1,0)] = rho*v;
	F->F[transform_F(F,pos_x,pos_y,pos_z,2,0)] = rho*w;
	F->F[transform_F(F,pos_x,pos_y,pos_z,0,1)] = rho*pow(u, 2) + p;
	F->F[transform_F(F,pos_x,pos_y,pos_z,1,1)] = rho*u*v;
	F->F[transform_F(F,pos_x,pos_y,pos_z,2,1)] = rho*u*w;
	F->F[transform_F(F,pos_x,pos_y,pos_z,0,2)] = rho*u*v;
	F->F[transform_F(F,pos_x,pos_y,pos_z,1,2)] = rho*pow(v, 2) + p;
	F->F[transform_F(F,pos_x,pos_y,pos_z,2,2)] = rho*v*w;
	F->F[transform_F(F,pos_x,pos_y,pos_z,0,3)] = rho*u*w;
	F->F[transform_F(F,pos_x,pos_y,pos_z,1,3)] = rho*v*w;
	F->F[transform_F(F,pos_x,pos_y,pos_z,2,3)] = rho*pow(w, 2) + p;
	F->F[transform_F(F,pos_x,pos_y,pos_z,0,4)] = rho*(E+p/rho)*u;
	F->F[transform_F(F,pos_x,pos_y,pos_z,1,4)] = rho*(E+p/rho)*v;
	F->F[transform_F(F,pos_x,pos_y,pos_z,2,4)] = rho*(E+p/rho)*w;
}

//One iteration for U-update
void newU(U_grid *U, F_grid *F, physics_grid *P, int pos_x, int pos_y, int pos_z, FLOAT dt)
{
	FLOAT *F_nhjk,*F_phjk,*F_inhk,*F_iphk,*F_ijnh,*F_ijph;
	if(!(F_nhjk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT))))
	{
		fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		exit(0);
	}
	if(!(F_phjk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT))))
	{
		fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		exit(0);
	}
	if(!(F_inhk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT))))
	{
		fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		exit(0);
	}
	if(!(F_iphk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT))))
	{
		fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		exit(0);
	}
	if(!(F_ijnh = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT))))
	{
		fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		exit(0);
	}
	if(!(F_ijph = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT))))
	{
		fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		exit(0);
	}
	int i;
	
	/*
	 * Error might be here.
	 * Posible reasons include:
	 * -> Interpolation scheme
	 * -> Bad if cases
	 */
	if(pos_x > 0 && pos_x < F->N_x) //middle case
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_nhjk[i] = (-F->F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]);
			F_phjk[i] = (-F->F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]);
		}
	}
	else if(pos_x == 0) //lower corner case
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_nhjk[i] = (-F->F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
			F_phjk[i] = (6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]);
		}
	}
	else if(pos_x == F->N_x) //higher corner case
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_nhjk[i] = (6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]);
			F_phjk[i] = (-F->F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
		}
	}
	
	if(pos_y > 0 && pos_y < F->N_y)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_inhk[i] = (-F->F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]);
			F_iphk[i] = (-F->F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]);
		}
	}
	else if(pos_y == 0)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_inhk[i] = (-F->F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
			F_iphk[i] = (6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]);
		}
	}
	else if(pos_y == F->N_y)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_inhk[i] = (6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]);
			F_iphk[i] = (-F->F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
		}
	}
	
	if(pos_z > 0 && pos_z < F->N_z)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_ijnh[i] = (-F->F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]);
			F_ijph[i] = (-F->F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]);
		}
	}
	else if(pos_z == 0)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_ijnh[i] = (-F->F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
			F_ijph[i] = (6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]);
		}
	}
	else if(pos_z == F->N_z)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_ijnh[i] = (6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F->F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]);
			F_ijph[i] = (-F->F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]+6*F->F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
		}
	}
	
	for(i=0;i<NDIM+2;i++)
	{
		U->U[transform_U(U,pos_x,pos_y,pos_z,i)] += dt*(F_nhjk[i]+F_inhk[i]+F_ijnh[i]-F_phjk[i]-F_iphk[i]-F_ijph[i])/P->delta_x;
	}
	
	free(F_nhjk);
	free(F_phjk);
	free(F_inhk);
	free(F_iphk);
	free(F_ijnh);
	free(F_ijph);
}
