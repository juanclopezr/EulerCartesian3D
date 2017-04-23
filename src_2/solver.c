#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "init.h"

//Calculates the value of e
// FLOAT calce(FLOAT E, FLOAT rho, FLOAT u, FLOAT v, FLOAT w)
FLOAT calce(FLOAT E_rho, FLOAT u, FLOAT v, FLOAT w)
{
	// return E_rho - 0.5*(pow(u, 2) + pow(v, 2) + pow(w, 2));
	FLOAT speeds = 0.5*(pow(u, 2) + pow(v, 2) + pow(w, 2));
	FLOAT ans = E_rho - speeds;
	// if (ans > 0)
	// {
	// 	return ans;
	// }
	// printf("%f %f %f\n", E_rho, speeds, ans);
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
	return p/rho;
	// return E + p/rho;
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
	return sqrt(h*GAMMA);
	// return sqrt(h*(GAMMA+1));
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
	P->E = U.components[4];
	P->u = U.components[1]/P->rho;
	P->v = U.components[2]/P->rho;
	P->w = U.components[3]/P->rho;
	e = calce(P->E/P->rho, P->u, P->v, P->w);
	P->p = calcp(P->rho, e);
}

void fromP_calcU(physics_cell P, U_vector *U)
{
	U->components[0] = P.rho;
	U->components[1] = (P.rho)*(P.u);
	U->components[2] = (P.rho)*(P.v);
	U->components[3] = (P.rho)*(P.w);
	U->components[4] = P.E;
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
	E = U->components[4];
	u = U->components[1]/rho;
	v = U->components[2]/rho;
	w = U->components[3]/rho;
	e = calce(E/rho, u, v, w);
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

/*F_vector *calculateF_half(U_grid *U, int i, int j, int k, int i_, int j_, int k_)
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
}*/
F_vector *calculateF_half(U_grid *U, int i, int j, int k, int i_, int j_, int k_)
{
	F_vector *F = create_F_vector();
	F_vector *Fp = create_F_vector();
	U_vector U1;
	int m;
	FLOAT temp;
	
	U1 = U->U[transform3d(i, j, k)];
	
	if(i_ != 0)
	{
		fromU_calcFx(U1, F);
	}
	else if(j_ != 0)
	{
		fromU_calcFy(U1, F);
	}
	else if(k_ != 0)
	{
		fromU_calcFz(U1, F);
	}
	
	if((i+i_>=0) && (i+i_< U->N_x))
	{
		if((j+j_>=0) && (j+j_< U->N_y))
		{
			if((k+k_>=0) && (k+k_< U->N_z))
			{
				U1 = U->U[transform3d(i + i_, j + j_, k + k_)];
			}
		}
	}
	
	if(i_ != 0)
	{
		fromU_calcFx(U1, Fp);
	}
	else if(j_ != 0)
	{
		fromU_calcFy(U1, Fp);
	}
	else if(k_ != 0)
	{
		fromU_calcFz(U1, Fp);
	}
	
	temp = 0.5*(F->components[0]+Fp->components[0]);
	F->components[0] = temp;
	temp = 0.5*(F->components[1]+Fp->components[1]);
	F->components[1] = temp;
	temp = 0.5*(F->components[2]+Fp->components[2]);
	F->components[2] = temp;
	temp = 0.5*(F->components[3]+Fp->components[3]);
	F->components[3] = temp;
	temp = 0.5*(F->components[4]+Fp->components[4]);
	F->components[4] = temp;
	
	destruct_F_vector(Fp);
	
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
	double dx, dy, dz;
	dx = P->delta_x;
	dy = P->delta_y;
	dz = P->delta_z;
	F_vector *Fx, *Fxl, *Fxr;
	F_vector *Fy, *Fyl, *Fyr;
	F_vector *Fz, *Fzl, *Fzr;

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
					U_temp.U[index].components[l] = U->U[index].components[l]
					- (dt/dx)*Fx->components[l] + (dt/dy)*Fy->components[l]
					- (dt/dz)*Fz->components[l];
				}

				UV_temp = U_temp.U[index];
				speed = max_speed(&UV_temp);
				if(speed > max)
				{
					max = speed;
				}
				//max = speed;

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
