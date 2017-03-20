#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//Transform function for accessing the U-indices
int transform_U(U_grid *U, int pos_x, int pos_y, int pos_z, int prop)
{
	return ((pos_x*U->N_x + pos_y)*U->N_y + pos_z)*U->N_z + prop
}

//Extracts the value of rho for a given cell
FLOAT extract_rho(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[transform_U(U,pos_x,pos_y,pos_z,0)];
}

//Extracts the value of u for a given cell
FLOAT extract_u(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[transform_U(U,pos_x,pos_y,pos_z,1)]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of v for a given cell
FLOAT extract_v(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[transform_U(U,pos_x,pos_y,pos_z,2)]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of w for a given cell
FLOAT extract_w(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[transform_U(U,pos_x,pos_y,pos_z,3)]/extract_rho(U, pos_x, pos_y, pos_z);
}

//Extracts the value of E for a given cell
FLOAT extract_E(U_grid *U, int pos_x, int pos_y, int pos_z)
{
	return U[transform_U(U,pos_x,pos_y,pos_z,4)]/extract_rho(U, pos_x, pos_y, pos_z);
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

//Transform function for accessing the F-indices
int transform_F(F_grid *F, int pos_x, int pos_y, int pos_z, int pos_g, int prop)
{
	return (((pos_x*F->N_x + pos_y)*F->N_y + pos_z)*F->N_z + pos_g)*NDIM + prop
}

/*TODO
 * Use the fourth order interpolation scheme to calculate 'F-half' with neighboring cell values.
 * Define the iteration function which will be used in the main.*/
 
//For serial calculation of F-vectors.
/*void calcF(F_grid *F, U_grid *U)
{
	int i;
	int j;
	int k;
	FLOAT rho,u,v,w,E,e,p;
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
				F[transform_F(F,i,j,k,0,0)] = rho*u;
				F[transform_F(F,i,j,k,1,0)] = rho*v;
				F[transform_F(F,i,j,k,2,0)] = rho*w;
				F[transform_F(F,i,j,k,0,1)] = rho*u*u+p;
				F[transform_F(F,i,j,k,1,1)] = rho*u*v;
				F[transform_F(F,i,j,k,2,1)] = rho*u*w;
				F[transform_F(F,i,j,k,0,2)] = rho*u*v;
				F[transform_F(F,i,j,k,1,2)] = rho*v*v+p;
				F[transform_F(F,i,j,k,2,2)] = rho*v*w;
				F[transform_F(F,i,j,k,0,3)] = rho*u*w;
				F[transform_F(F,i,j,k,1,3)] = rho*v*w;
				F[transform_F(F,i,j,k,2,3)] = rho*w*w+p;
				F[transform_F(F,i,j,k,0,4)] = rho*(E+p/rho)*u;
				F[transform_F(F,i,j,k,1,4)] = rho*(E+p/rho)*v;
				F[transform_F(F,i,j,k,2,4)] = rho*(E+p/rho)*w;
			}
		}
	}
}*/

//Calculate F-vector
void calcF(F_grid *F, U_grid *U, int pos_x, int pos_y, int pos_z)
{
	FLOAT rho,u,v,w,E,e,p;
	rho = extract_rho(U,pos_x,pos_y,pos_z);
	u = extract_u(U,pos_x,pos_y,pos_z);
	v = extract_v(U,pos_x,pos_y,pos_z);
	w = extract_w(U,pos_x,pos_y,pos_z);
	E = extract_E(U,pos_x,pos_y,pos_z);
	e = calce(E,u,v,w);
	p = calcp(rho,e);
	F[transform_F(F,pos_x,pos_y,pos_z,0,0)] = rho*u;
	F[transform_F(F,pos_x,pos_y,pos_z,1,0)] = rho*v;
	F[transform_F(F,pos_x,pos_y,pos_z,2,0)] = rho*w;
	F[transform_F(F,pos_x,pos_y,pos_z,0,1)] = rho*u*u+p;
	F[transform_F(F,pos_x,pos_y,pos_z,1,1)] = rho*u*v;
	F[transform_F(F,pos_x,pos_y,pos_z,2,1)] = rho*u*w;
	F[transform_F(F,pos_x,pos_y,pos_z,0,2)] = rho*u*v;
	F[transform_F(F,pos_x,pos_y,pos_z,1,2)] = rho*v*v+p;
	F[transform_F(F,pos_x,pos_y,pos_z,2,2)] = rho*v*w;
	F[transform_F(F,pos_x,pos_y,pos_z,0,3)] = rho*u*w;
	F[transform_F(F,pos_x,pos_y,pos_z,1,3)] = rho*v*w;
	F[transform_F(F,pos_x,pos_y,pos_z,2,3)] = rho*w*w+p;
	F[transform_F(F,pos_x,pos_y,pos_z,0,4)] = rho*(E+p/rho)*u;
	F[transform_F(F,pos_x,pos_y,pos_z,1,4)] = rho*(E+p/rho)*v;
	F[transform_F(F,pos_x,pos_y,pos_z,2,4)] = rho*(E+p/rho)*w;
}

	//FLOAT h,cs,sps,dt; To be used in iterations
	//h = calch(E,p,rho);
	//cs = calcs(h);
	//sps = calcsps(U,cs);

//One iteration for U-update
void newU(U_grid *U, F_grid, int pos_x, int pos_y, int pos_z, FLOAT dt, FLOAT dx)
{
	FLOAT *F_nhjk,*F_phjk,*F_inhk,*F_iphk,*F_ijnh,*Fijph;
	F_nhjk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT));
	F_phjk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT));
	F_inhk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT));
	F_iphk = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT));
	F_ijnh = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT));
	F_ijph = (FLOAT *)malloc((NDIM+2)*sizeof(FLOAT));
	int i;
	if(pos_x > 0 && pos_x < F->N_x)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_nhjk[i] = (-F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]);
			F_phjk[i] = (-F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]);
		}
	}
	else if(pos_x == 0)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_nhjk[i] = (-F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
			F_phjk[i] = (6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x+1,pos_y,pos_z,0,i)]);
		}
	}
	else if(pos_x == F->N_x)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_nhjk[i] = (6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]);
			F_phjk[i] = (-F[transform_F(F,pos_x-1,pos_y,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
		}
	}
	
	if(pos_y > 0 && pos_y < F->N_y)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_inhk[i] = (-F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]);
			F_iphk[i] = (-F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]);
		}
	}
	else if(pos_y == 0)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_inhk[i] = (-F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
			F_iphk[i] = (6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y+1,pos_z,0,i)]);
		}
	}
	else if(pos_y == F->N_y)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_inhk[i] = (6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]);
			F_iphk[i] = (-F[transform_F(F,pos_x,pos_y-1,pos_z,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
		}
	}
	
	if(pos_z > 0 && pos_z < F->N_z)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_ijnh[i] = (-F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]);
			F_ijph[i] = (-F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]);
		}
	}
	else if(pos_z == 0)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_ijnh[i] = (-F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
			F_ijph[i] = (6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y,pos_z+1,0,i)]);
		}
	}
	else if(pos_z == F->N_z)
	{
		for(i=0;i<NDIM+2;i++)
		{
			F_ijnh[i] = (6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]+3*F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]);
			F_ijph[i] = (-F[transform_F(F,pos_x,pos_y,pos_z-1,0,i)]+6*F[transform_F(F,pos_x,pos_y,pos_z,0,i)]);
		}
	}
	
	for(i=0;i<NDIM+2;i++)
	{
		U[transform_U(U,pos_x,pos_y,pos_z,i)] += dt*(F_nhjk[i]+F_inhk[i]+F_ijnh[i]-F_phjk[i]-F_iphk[i]-F_ijph[i])/dx;
	}
	
	free(F_nhjk);
	free(F_phjk);
	free(F_inhk);
	free(F_iphk);
	free(F_ijnh);
	free(F_ijph);
}
