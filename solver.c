#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "init.h"

//Calculates the value of e
FLOAT calce(FLOAT E, FLOAT u, FLOAT v, FLOAT w)
{
	// return E_rho - 0.5*(pow(u, 2) + pow(v, 2) + pow(w, 2));
	FLOAT speeds = 0.5*(pow(u, 2) + pow(v, 2) + pow(w, 2));
	FLOAT ans = E - speeds;
	return ans;
}

//Calculates the pressure value
FLOAT calcp(FLOAT rho, FLOAT e)
{
	return rho*e*(GAMMA-1);
}

//Calculates h
FLOAT calch(FLOAT E, FLOAT p, FLOAT rho)
{
	// return p/rho;
	return E + p/rho;
}

FLOAT absolute(FLOAT value)
{
	if(value >= 0)
	{
		return value;
	}
	return -value;
}

//Calculates cs
FLOAT calcs(FLOAT h)
{
	// return sqrt(h*GAMMA);
	return sqrt(h*(GAMMA+1));
}

// Transform_3d
int transform3d(int pos_x, int pos_y, int pos_z)
{
	return Ny*Nx*pos_x + Nx*pos_y + pos_z;
}

void fromU_calcP(U_vector U, physics_cell *P)
{
	FLOAT e;
	P->rho = U.components[0];
	P->E = U.components[4]/P->rho;
	P->u = U.components[1]/P->rho;
	P->v = U.components[2]/P->rho;
	P->w = U.components[3]/P->rho;
	e = calce(P->E, P->u, P->v, P->w);
	P->p = calcp(P->rho, e);
}

void fromP_calcU(physics_cell P, U_vector *U)
{
	U->components[0] = P.rho;
	U->components[1] = (P.rho)*(P.u);
	U->components[2] = (P.rho)*(P.v);
	U->components[3] = (P.rho)*(P.w);
	U->components[4] = P.E*P.rho;
}

void fromU_calcFx(U_vector U, F_vector *F)
{
	FLOAT u, v, w, p, e;
	u = U.components[1]/U.components[0];
	v = U.components[2]/U.components[0];
	w = U.components[3]/U.components[0];
	e = calce(U.components[4]/U.components[0], u, v, w);
	p = calcp(U.components[0], e);
	F->components[0] = U.components[1];
	F->components[1] = U.components[1]*u + p;
	F->components[2] = U.components[2]*u;
	F->components[3] = U.components[3]*u;
	F->components[4] = (U.components[4] + p)*u;
}

void fromU_calcFy(U_vector U, F_vector *F)
{
	FLOAT u, v, w, p, e;
	u = U.components[1]/U.components[0];
	v = U.components[2]/U.components[0];
	w = U.components[3]/U.components[0];
	e = calce(U.components[4]/U.components[0], u, v, w);
	p = calcp(U.components[0], e);
	F->components[0] = U.components[2];
	F->components[1] = U.components[1]*v;
	F->components[2] = U.components[2]*v + p;
	F->components[3] = U.components[3]*v;
	F->components[4] = (U.components[4] + p)*v;
}

void fromU_calcFz(U_vector U, F_vector *F)
{
	FLOAT u, v, w, p, e;
	u = U.components[1]/U.components[0];
	v = U.components[2]/U.components[0];
	w = U.components[3]/U.components[0];
	e = calce(U.components[4]/U.components[0], u, v, w);
	p = calcp(U.components[0], e);
	F->components[0] = U.components[3];
	F->components[1] = U.components[1]*w;
	F->components[2] = U.components[2]*w;
	F->components[3] = U.components[3]*w + p;
	F->components[4] = (U.components[4] + p)*w;
}
//
FLOAT max_speed(U_vector *U)
{
	FLOAT rho, E, u, v, w, p, e, cs, max;
	rho = U->components[0];
	E = U->components[4]/rho;
	u = U->components[1]/rho;
	v = U->components[2]/rho;
	w = U->components[3]/rho;
	e = calce(E, u, v, w);
	p = calcp(rho, e);
	cs = calcs(calch(E, p, rho));
	u = absolute(u);
	v = absolute(v);
	w = absolute(w);
	max = u;
	if (v > max)
	{
		max = v;
	}
	if (w > max)
	{
		max = w;
	}
	return max + cs;
}

void fromU_calcF(U_vector U, F_vector *Fx, F_vector *Fy, F_vector *Fz)
{
	fromU_calcFx(U, Fx);
	fromU_calcFy(U, Fy);
	fromU_calcFz(U, Fz);
}

U_vector *calculateU_half(int ver, U_vector U1, U_vector U2)
{
	int i;
	U_vector *U = create_U_vector();


	if (ver==1)
	{
		for(i = 0; i<(NDIM+2); i++)
		{
			U->components[i] = 0.5*(U1.components[i] + U2.components[i]);
		}
	}
	else{
		for(i = 0; i<(NDIM+2); i++)
		{
			U->components[i] = U1.components[i];
		}
	}

	return U;
}

F_vector *calculateF_half(U_grid *U, int i, int j, int k, int i_, int j_, int k_)
{
	F_vector *F = create_F_vector();
	U_vector U1;
	U_vector U2;
	U_vector *U_;

	int ver = 0;
	U1 = U->U[transform3d(i,j,k)];

	if((i+i_>=0) && (i+i_< U->N_x))
	{
		if((j+j_>=0) && (j+j_< U->N_y))
		{
			if((k+k_>=0) && (k+k_< U->N_z))
			{
				ver = 1;
				U2 = U->U[transform3d(i + i_, j + j_, k + k_)];
			}
		}
	}
	U_ = calculateU_half(ver, U1, U2);
	if(i_ != 0)
	{
		fromU_calcFx(*U_, F);
	}
	else if(j_ != 0)
	{
		fromU_calcFy(*U_, F);
	}
	else if(k_ != 0)
	{
		fromU_calcFz(*U_, F);
	}

	destruct_U_vector(U_);
	return F;
}

F_vector *calculateF_diff(F_vector F_plus, F_vector F_minus)
{
	int i;
	F_vector *F = create_F_vector();
	for(i = 0; i<(NDIM+2); i++)
	{
		F->components[i] = F_plus.components[i] - F_minus.components[i];
	}
	return F;
}

void updateU_grid(U_grid *U, U_grid U_temp)
{
	int i, j, k, l, index;
	for(i=0; i<U->N_x; i++)
	{
		for(j=0; j<U->N_y; j++)
		{
			for(k=0; k<U->N_z; k++)
			{
				index = transform3d(i,j,k);
				for(l=0; l<(NDIM+2); l++)
				{
					U->U[index].components[l] = U_temp.U[index].components[l];
				}
			}
		}
	}
}

FLOAT calculateNextU(physics_grid *P, U_grid *U, U_grid U_temp, double dt)
{
	int i, j, k, l, index;
	FLOAT dtdx, dtdy, dtdz;
	dtdx = dt/P->delta_x;
	dtdy = dt/P->delta_y;
	dtdz = dt/P->delta_z;
	F_vector *Fx, *Fxl, *Fxr;
	F_vector *Fy, *Fyl, *Fyr;
	F_vector *Fz, *Fzl, *Fzr;
	// FLOAT values;
	FLOAT max, speed;
	U_vector UV_temp;
	max = 0;
	for(i=0; i<U->N_x; i++)
	{
		for(j=0; j<U->N_y; j++)
		{
			for(k=0; k<U->N_z; k++)
			{
				index = transform3d(i, j, k);
				Fxr = calculateF_half(U, i, j, k, 1, 0, 0);
				Fxl = calculateF_half(U, i, j, k, -1, 0, 0);
				Fyr = calculateF_half(U, i, j, k, 0, 1, 0);
				Fyl = calculateF_half(U, i, j, k, 0, -1, 0);
				Fzr = calculateF_half(U, i, j, k, 0, 0, 1);
				Fzl = calculateF_half(U, i, j, k, 0, 0, -1);

				Fx = calculateF_diff(*Fxr, *Fxl);
				Fy = calculateF_diff(*Fyr, *Fyl);
				Fz = calculateF_diff(*Fzr, *Fzl);
				for(l=0; l<(NDIM+2); l++)
				{
					// values = U->U[transform3d(i+1, j, k)].components[l] +
					// 		U->U[transform3d(i-1, j, k)].components[l] +
					// 		U->U[transform3d(i, j+1, k)].components[l] +
					// 		U->U[transform3d(i, j-1, k)].components[l] +
					// 		U->U[transform3d(i, j, k+1)].components[l] +
					// 		U->U[transform3d(i, j, k-1)].components[l] -
					// 		6*U->U[index].components[l];
					// U_temp.U[index].components[l] = U->U[index].components[l] + dtdx*values;
					U_temp.U[index].components[l] = U->U[index].components[l]
					+ ((dtdx)*Fx->components[l] + (dtdy)*Fy->components[l]
					+ (dtdz)*Fz->components[l]);
				}

				UV_temp = U_temp.U[index];
				speed = max_speed(&UV_temp);
				if(speed > max)
				{
					max = speed;
				}

				destruct_F_vector(Fx);
				destruct_F_vector(Fxl);
				destruct_F_vector(Fxr);
				destruct_F_vector(Fy);
				destruct_F_vector(Fyl);
				destruct_F_vector(Fyr);
				destruct_F_vector(Fz);
				destruct_F_vector(Fzl);
				destruct_F_vector(Fzr);
			}
		}
	}
	updateU_grid(U, U_temp);
	return max;
}

void updatePhysics(physics_grid *P, U_grid U)
{
	int i, j, k, index;
	physics_cell P_;
	U_vector U_;
	for(i=0; i<U.N_x; i++)
	{
		for(j=0; j<U.N_y; j++)
		{
			for(k=0; k<U.N_z; k++)
			{
				index = transform3d(i,j,k);
				P_ = P->P[index];
				U_ = U.U[index];
				fromU_calcP(U_, &P_);
				P->P[index] = P_;
			}
		}
	}
}

FLOAT *calc_density(physics_grid P)
{
	FLOAT *density = malloc(0.5*P.N_x*sizeof(FLOAT));
	int i, j, k, l, i_c, j_c, k_c, index, contador;
	FLOAT temp_density;
	i_c = 0.5*P.N_x;
	j_c = 0.5*P.N_x;
	k_c = 0.5*P.N_x;
	for(l=0; l<i_c; l++)
	{
		temp_density = 0;
		contador = 0;
		for(i=-l; i<l; i++)
		{
			for(j=-l; j<l; j++)
			{
				for(k=-l; k<l; k++)
				{
					if(abs(i)+abs(j)+abs(k) >= l)
					{
						index = transform3d(i+i_c, j+j_c, k+k_c);
						temp_density += P.P[index].rho;
						contador += 1;
					}
				}
			}
		}
		density[l] = temp_density/contador;
	}
	return density;
}
