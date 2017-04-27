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

U_vector *calcualteU_fourth(U_grid *U, int i, int j, int k, int comp, int compp, int disp, int term, FLOAT dt)
{
	U_vector *Ucal = create_U_vector();
	U_vector U1;
	F_vector *Fcal = create_F_vector();
	
	if(comp==1)
	{
		if(term==1)
		{
			U1 = U->U[transform3d(i+disp,j,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j+1,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j,k+1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j+1,k+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i+disp,j+1,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k+1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
		else if(term==2)
		{
			U1 = U->U[transform3d(i+disp,j,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j-1,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j,k+1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j-1,k+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i+disp,j-1,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k+1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
		else if(term==3)
		{
			U1 = U->U[transform3d(i+disp,j,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j+1,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j,k-1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j+1,k-1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i+disp,j+1,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k-1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
		else
		{
			U1 = U->U[transform3d(i+disp,j,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j-1,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j,k-1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+disp,j-1,k-1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i+disp,j-1,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,1,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k-1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+disp,j,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
	}
	else if(comp==2)
	{
		if(term==1)
		{
			U1 = U->U[transform3d(i,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j+disp,k+1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j+disp,k+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i+1,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+1,j+disp,k+1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+1,j+disp,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
		else if(term==2)
		{
			U1 = U->U[transform3d(i,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j+disp,k+1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j+disp,k+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i-1,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i-1,j+disp,k+1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i-1,j+disp,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
		else if(term==3)
		{
			U1 = U->U[transform3d(i,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j+disp,k-1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j+disp,k-1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i+1,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+1,j+disp,k-1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+1,j+disp,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
		else
		{
			U1 = U->U[transform3d(i,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j+disp,k)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j+disp,k-1)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j+disp,k-1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i-1,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j+disp,k,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i-1,j+disp,k-1,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i-1,j+disp,k,2,comp,compp,term,dt);
			FromU_calcFz(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_z);
			subsUF(U1,Fcal);
		}
	}
	//fix
	else
	{
		if(term==1)
		{
			U1 = U->U[transform3d(i,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j+1,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j+1,k+disp+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+1,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j+1,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
		}
		else if(term==2)
		{
			U1 = U->U[transform3d(i,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j+1,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j+1,k+disp+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i-1,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j+1,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
		}
		else if(term==3)
		{
			U1 = U->U[transform3d(i,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j-1,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i+1,j-1,k+disp+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i+1,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j-1,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
		}
		else
		{
			U1 = U->U[transform3d(i,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i,j-1,k+disp)];
			addU(Ucal,U1);
			U1 = U->U[transform3d(i-1,j-1,k+disp+1)];
			addU(Ucal,U1);
			U1 = calculateU_sixth(U,i,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i-1,j,k+disp,1,comp,compp,term,dt);
			FromU_calcFx(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_x);
			subsUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			addUF(U1,Fcal);
			U1 = calculateU_sixth(U,i,j-1,k+disp,2,comp,compp,term,dt);
			FromU_calcFy(U1,Fcal);
			factor(Fcal,0.25*dt/U->delta_y);
			subsUF(U1,Fcal);
		}
	}
	destruct_F_vector(Fcal);
	factorU(Ucal,0.25);
	return Ucal;
}
	

U_vector *calculateU_half(U_grid *U, int i, int j, int k, int edgex, int edgey, int edgez, int comp, int term, int dir, FLOAT dt)
{
	U_vector *Ucal = create_U_vector();
	U_vector U1;
	F_vector *Fcal = create_F_vector();
	
	U1 = U->U[transform3d(i,j,k)];
	addU(Ucal,U1);
	
	if(comp==1)
	{
		if(term==1)
		{
			if(edgex*dir<=0)
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
					}
					else
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
			}
		}
		else if(term==2)
		{
			if(edgex*dir<=0)
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
					}
					else
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
				}
			}
		}
		else if(term==3)
		{
			if(edgex*dir<=0)
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j+1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j+1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else
					{
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
				}
			}
		}
		else
		{
			if(edgex*dir<=0)
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j-1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,1,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,dir,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+dir,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+dir,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+dir,j,k)];
						addU(Ucal,U1);
					}
				}
			}
			else
			{
				if(edgey==0)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
				}
				else if(edgey==1)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else
					{
					}
				}
			}
		}
	}
	else if(comp==2)
	{
		if(term==1)
		{
			if(edgey*dir<=0)
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
					}
					else
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
			}
		}
		else if(term==2)
		{
			if(edgey*dir<=0)
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
					}
					else
					{
						U1 = U->U[transform3d(i,j,k+1)];
						addU(Ucal,U1);
					}
				}
			}
		}
		else if(term==3)
		{
			if(edgey*dir<=0)
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{	
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						
					}
					else
					{
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
				}
			}
		}
		else
		{
			if(edgey*dir<=0)
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+dir,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+dir,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = u->u[transform3d(i,j+dir,k)];
						addU(Ucal,U1);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
				}
				else if(edgex==1)
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k-1)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgez==0)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else if(edgez==1)
					{
						U1 = U->U[transform3d(i,j,k-1)];
						addU(Ucal,U1);
					}
					else
					{
					}
				}
			}
		}
	}
	else
	{
		if(term==1)
		{
			if(edgez*dir<=0)
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else if(edgey==1)
					{
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
			}
		}
		else if(term==2)
		{
			if(edgez*dir<=0)
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j+1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j+1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
					else if(edgey==1)
					{
					}
					else
					{
						U1 = u->u[transform3d(i,j+1,k)];
						addU(Ucal,U1);
					}
				}
			}
		}
		else if(term==3)
		{
			if(edgez*dir<=0)
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{	
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						
					}
					else if(edgey==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i+1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
					else if(edgey==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
					else
					{
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i+1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i+1,j,k)];
						addU(Ucal,U1);
					}
				}
			}
		}
		else
		{
			if(edgez*dir<=0)
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{	
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
						
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,1,1,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,1,dir,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,1,1,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						addUF(U1,Fcal);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i-1,j,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,2,2,0,term,dt);
						fromU_calcFy(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_y);
						subsUF(U1,Fcal);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
						U1 = U->U[transform3d(i,j-1,k+dir)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,1,2,0,term,dt);
						fromU_calcFx(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_x);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i,j,k+dir)];
						addU(Ucal,U1);
					}
				}
			}
			else
			{
				if(edgex==0)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
				}
				else if(edgex==1)
				{
					if(edgey==0)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else if(edgey==1)
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
						U1 = u->u[transform3d(i-1,j-1,k)];
						addU(Ucal,U1);
						U1 = calculateUfourth(U,i,j,k,3,2,0,term,dt);
						fromU_calcFz(U1,Fcal);
						factor(Fcal,0.5*dt/U->delta_z);
						subsUF(U1,Fcal);
					}
					else
					{
						U1 = U->U[transform3d(i-1,j,k)];
						addU(Ucal,U1);
					}
				}
				else
				{
					if(edgey==0)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
					else if(edgey==1)
					{
						U1 = u->u[transform3d(i,j-1,k)];
						addU(Ucal,U1);
					}
					else
					{
					}
				}
			}
		}
	}
	destruct_F_vector(Fcal);
	FactorU(Ucal,1./8.);
	return Ucal;
}
	

F_vector *calculateF_half(U_grid *U, int i, int j, int k, int i_, int j_, int k_, FLOAT dt)
{
	F_vector *F = create_F_vector();
	U_vector *U_;
	F_vector *Fcal = create_F_vector();
	
	int edgex = 0;
	int edgey = 0;
	int edgez = 0;
	
	if(i==0)
	{
		edgex = 1;
	}
	else if(i==U->N_x-1)
	{
		edgex = -1;
	}
	
	if(j==0)
	{
		edgey = 1;
	}
	else if(j==U->N_y-1)
	{
		edgey = -1;
	}
	
	if(k==0)
	{
		edgez = 1;
	}
	else if(k==U->N_z-1)
	{
		edgez = -1;
	}
	
	if(i_!=0)
	{
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,1,1,i_,dt);
		FromU_calcFx(*U_, Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,1,2,i_,dt);
		FromU_calcFx(*U_, Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,1,3,i_,dt);
		FromU_calcFx(*U_, Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,1,4,i_,dt);
		FromU_calcFx(*U_, Fcal);
		addF(F,Fcal);
		factor(F,0.25);
	}
	if(j_!=0)
	{
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,2,1,j_,dt);
		FromU_calcFy(*U_,Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,2,2,j_,dt);
		FromU_calcFy(*U_,Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,2,3,j_,dt);
		FromU_calcFy(*U_,Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,2,4,j_,dt);
		FromU_calcFy(*U_,Fcal);
		addF(F,Fcal);
		factor(F,0.25);
	}
	if(k_!=0)
	{
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,3,1,k_,dt);
		FromU_calcFz(*U_,Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,3,2,k_,dt);
		FromU_calcFz(*U_,Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,3,3,k_,dt);
		FromU_calcFz(*U_,Fcal);
		addF(F,Fcal);
		U_ = calculateU_half(U,i,j,k,edgex,edgey,edgez,3,4,k_,dt);
		FromU_calcFz(*U_,Fcal);
		addF(F,Fcal);
		factor(F,0.25);
	}
	destruct_F_vector(Fcal);
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
	double dtdx, dtdy, dtdz;
	dtdx = dt/P->delta_x;
	dtdy = dt/P->delta_y;
	dtdz = dt/P->delta_z;
	F_vector *Fx, *Fxl, *Fxr;
	F_vector *Fy, *Fyl, *Fyr;
	F_vector *Fz, *Fzl, *Fzr;

	FLOAT max, speed;
	U_vector UV_temp;
	max = 0;
	int middle = transform3d(0.5*U->N_x, 0.5*U->N_y, 0.5*U->N_z);
	for(i=0; i<U->N_x; i++)
	{
		for(j=0; j<U->N_y; j++)
		{
			for(k=0; k<U->N_z; k++)
			{
				index = transform3d(i, j, k);
				Fxr = calculateF_half(U, dtdx, i, j, k, 1, 0, 0);
				Fxl = calculateF_half(U, dtdx, i, j, k, -1, 0, 0);
				Fyr = calculateF_half(U, dtdy, i, j, k, 0, 1, 0);
				Fyl = calculateF_half(U, dtdy, i, j, k, 0, -1, 0);
				Fzr = calculateF_half(U, dtdz, i, j, k, 0, 0, 1);
				Fzl = calculateF_half(U, dtdz, i, j, k, 0, 0, -1);

				Fx = calculateF_diff(*Fxr, *Fxl);
				Fy = calculateF_diff(*Fyr, *Fyl);
				Fz = calculateF_diff(*Fzr, *Fzl);
				if(1)
				{
					for(l=0; l<(NDIM+2); l++)
					{
						U_temp.U[index].components[l] = U->U[index].components[l]
						+ ((dtdx)*Fx->components[l] + (dtdy)*Fy->components[l]
						+ (dtdz)*Fz->components[l]);
					}
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

void fromUgrid_calcFgrid(U_grid Ugrid, F_grid Fgrid)
{
	U_vector U;
	F_vector Fx, Fy, Fz;

	int i, j, k, index;
	for(i=0; U->N_x; i++)
	{
		for(j=0; U->N_y; j++)
		{
			for(k=0; U->N_z; k++)
			{
				index = transform3d(i, j, k);
				U = Ugrid.U[index];
				Fx = Fgrid.F_x[index];
				Fy = Fgrid.F_y[index];
				Fz = Fgrid.F_z[index];

				fromU_calcF(*U, *Fx, *Fy, *Fz);
			}
		}
	}
}
