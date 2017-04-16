#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "init.h"

//Calculates the value of e
FLOAT calce(FLOAT E, FLOAT u, FLOAT v, FLOAT w)
{
	return E - 0.5*(pow(u, 2) + pow(v, 2) + pow(w, 2));
}

//Calculates the pressure value
FLOAT calcp(FLOAT rho, FLOAT e)
{
	return rho*e*(GAMMA-1);
}

//Calculates h
FLOAT calch(FLOAT E, FLOAT p, FLOAT rho)
{
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
	return sqrt(h*(GAMMA+1));
}

// Transform_3d
int transform3d(int pos_x, int pos_y, int pos_z)
{
	return Ny*Nx*pos_x + Nx*pos_y + pos_z;
}

void fromU_calcP(U_vector *U, physics_cell *P)
{
	FLOAT e;
	P->rho = U->components[0];
	P->E = U->components[4]/P->rho;
	P->u = U->components[1]/P->rho;
	P->v = U->components[2]/P->rho;
	P->w = U->components[3]/P->rho;
	e = calce(P->E, P->u, P->v, P->w);
	P->p = calcp(P->rho, e);
}

void fromP_calcU(physics_cell *P, U_vector *U)
{
	U->components[0] = P->rho;
	U->components[1] = (P->rho)*(P->u);
	U->components[2] = (P->rho)*(P->v);
	U->components[3] = (P->rho)*(P->w);
	U->components[4] = (P->rho)*(P->E);
}

void fromU_calcFx(U_vector *U, F_vector *F)
{
	FLOAT u, v, w, p, e;
	u = U->components[1]/U->components[0];
	v = U->components[2]/U->components[0];
	w = U->components[3]/U->components[0];
	e = calce(U->components[4]/U->components[0], u, v, w);
	p = calcp(U->components[0], e);
	F->components[0] = U->components[1];
	F->components[1] = U->components[1]*u + p;
	F->components[2] = U->components[2]*u;
	F->components[3] = U->components[3]*u;
	F->components[4] = (U->components[4] + p)*u;
}

void fromU_calcFy(U_vector *U, F_vector *F)
{
	FLOAT u, v, w, p, e;
	u = U->components[1]/U->components[0];
	v = U->components[2]/U->components[0];
	w = U->components[3]/U->components[0];
	e = calce(U->components[4]/U->components[0], u, v, w);
	p = calcp(U->components[0], e);
	F->components[0] = U->components[2];
	F->components[1] = U->components[1]*v;
	F->components[2] = U->components[2]*u + p;
	F->components[3] = U->components[3]*v;
	F->components[4] = (U->components[4] + p)*v;
}

void fromU_calcFz(U_vector *U, F_vector *F)
{
	FLOAT u, v, w, p, e;
	u = U->components[1]/U->components[0];
	v = U->components[2]/U->components[0];
	w = U->components[3]/U->components[0];
	e = calce(U->components[4]/U->components[0], u, v, w);
	p = calcp(U->components[0], e);
	F->components[0] = U->components[3];
	F->components[1] = U->components[1]*w;
	F->components[2] = U->components[2]*w;
	F->components[3] = U->components[3]*w + p;
	F->components[4] = (U->components[4] + p)*w;
}

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

void fromU_calcF(U_vector *U, F_vector *Fx, F_vector *Fy, F_vector *Fz)
{
	fromU_calcFx(U, Fx);
	fromU_calcFy(U, Fy);
	fromU_calcFz(U, Fz);
}

U_vector *calculateU_half(U_vector *U1, U_vector *U2)
{
	int i;
	U_vector *U = create_U_vector();
	for(i = 0; i<(NDIM+2); i++)
	{
		U->components[i] = 0.5*(U1->components[i] + U2->components[i]);
	}
	// destruct_U_vector(U1);
	// destruct_U_vector(U2);
	return U;
}

F_vector *calculateF_half(U_grid *U, int i, int j, int k, int i_, int j_, int k_)
{
	F_vector *F = create_F_vector();
	U_vector *U1 = create_U_vector();
	U_vector *U2 = create_U_vector();
	U_vector *U_;
	int l;

	for(l=0; l<(NDIM+2); l++)
	{
		U1->components[l] = U->U[transform3d(i,j,k)].components[l];
	}

	if((i+i_>=0) && (i+i_< U->N_x))
	{
		if((j+j_>=0) && (j+j_< U->N_y))
		{
			if((k+k_>=0) && (k+k_< U->N_z))
			{
				for(l=0; l<(NDIM+2); l++)
				{
					U2->components[l] = U->U[transform3d(i + i_, j + j_, k + k_)].components[l];
				}
			}
		}
	}

	U_ = calculateU_half(U1, U2);

	if(i_ != 0)
	{
		fromU_calcFx(U_, F);
	}
	else if(j_ != 0)
	{
		fromU_calcFy(U_, F);
	}
	else if(k_ != 0)
	{
		fromU_calcFz(U_, F);
	}

	// destruct_U_vector(U_);
	return F;
}

F_vector *calculateF_diff(F_vector *F_plus, F_vector *F_minus)
{
	int i;
	F_vector *F = create_F_vector();
	for(i = 0; i<(NDIM+2); i++)
	{
		F->components[i] = F_minus->components[i] - F_plus->components[i];
	}
	// destruct_F_vector(F_plus);
	// destruct_F_vector(F_minus);
	return F;
}

FLOAT calculateNextU(physics_grid *P, U_grid *U_new, U_grid *U_past, double dt)
{
	int i, j, k, l, index;
	double dx, dy, dz;
	dx = P->delta_x;
	dy = P->delta_y;
	dz = P->delta_z;
	F_vector *Fx;
	F_vector *Fy;
	F_vector *Fz;
	for(i=0; i<U_past->N_x; i++)
	{
		for(j=0; j<U_past->N_y; j++)
		{
			for(k=0; k<U_past->N_z; k++)
			{
				index = transform3d(i, j, k);
				Fx = calculateF_diff(calculateF_half(U_past, i, j, k, 1, 0, 0), \
					calculateF_half(U_past, i, j, k, -1, 0, 0));
				Fy = calculateF_diff(calculateF_half(U_past, i, j, k, 0, 1, 0), \
					calculateF_half(U_past, i, j, k, 0, -1, 0));
				Fz = calculateF_diff(calculateF_half(U_past, i, j, k, 0, 0, 1), \
					calculateF_half(U_past, i, j, k, 0, 0, -1));
				for(l=0; l<(NDIM+2); l++)
				{
					U_new->U[index].components[l] = U_past->U[index].components[l] \
					+ (dt/dx)*Fx->components[l] + (dt/dy)*Fy->components[l] \
					+ (dt/dz)*Fz->components[l];
				}
				destruct_F_vector(Fx);
				destruct_F_vector(Fy);
				destruct_F_vector(Fz);
			}
		}
	}


	FLOAT max, speed;
	U_vector *U_temp = create_U_vector();
	max = 0;
	// Update values
	for(i=0; i<U_past->N_x; i++)
	{
		for(j=0; j<U_past->N_y; j++)
		{
			for(k=0; k<U_past->N_z; k++)
			{
				index = transform3d(i,j,k);
				for(l=0; l<(NDIM+2); l++)
				{
					U_past->U[index].components[l] = U_new->U[index].components[l];
				}
				*U_temp = U_past->U[index];
				speed = max_speed(U_temp);
				if(speed > max)
				{
					max = speed;
				}
			}
		}
	}
	// destruct_U_vector(U_temp);
	return max;
}

void updatePhysics(physics_grid *P, U_grid *U)
{
	int i, j, k, index;
	physics_cell *P_ = create_P_vector();
	U_vector *U_ = create_U_vector();
	for(i=0; i<U->N_x; i++)
	{
		for(j=0; j<U->N_y; j++)
		{
			for(k=0; k<U->N_z; k++)
			{
				index = transform3d(i,j,k);
				*P_ = P->P[index];
				*U_ = U->U[index];
				fromU_calcP(U_, P_);
				P->P[index] = *P_;
				U->U[index] = *U_;
			}
		}
	}
	// free(P_);
	// destruct_U_vector(U_);
}
