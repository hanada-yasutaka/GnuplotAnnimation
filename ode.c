#include <stdio.h>

void Integrator(void (*gen)(double, double, double *, double *),
	 double t, double dt,double *x, double *xout);
void Euler(double t, double dt, double *x, double *xout);
void RK2(double t, double dt, double *x, double *xout);
void RK4(double t, double dt, double *x, double *xout);
void Symplectic(double t, double dt, double *x, double *xout);
void Symplectic2(double t, double dt, double *x, double *xout);

/*---- ユーザー定義関数 (ここから)：user defined functions ---- */

double Hamiltonian(double t, double q, double p)
{
	/* 
	definition of Hamiltonian
	t : time
	q : position
	p : momentum
	*/
	double omega=1.0;
	return 0.5*p*p + 0.5*omega*q*q;
}
double Diff_q(double t,double q, double p)
{
	/* 
	difinition of dq/dt
	t : time
	q : position
	p : momentum
	*/
	return p;
}
double Diff_p(double t, double q, double p)
{
	/* 
	difinition of dp/dt
	t : time
	q : position
	p : momentum
	*/
	double omega=1.0;
	return -omega*q;
}

/*---- ユーザー定義関数 (ここまで)：end user defined functions ---- */

int main(void)
{
	FILE *f;

	int i;
	int iteration;

	double x[2];
	double xx[2];	
	double t, dt;

	dt = 0.05;
	iteration=1000;

	x[0] = 1.0; // q
	x[1] = 0.0; // p

	f = fopen("trajectory_symp2.dat","w");

		for (i=0; i<iteration;i++){
		t = i*dt;
		fprintf(f,"%lf %lf %lf %le\n", t, x[0], x[1], Hamiltonian(t,x[0],x[1]));
		Integrator(Symplectic2, t,dt,x,xx);

		x[0] = xx[0];
		x[1] = xx[1];
	}
	fclose(f);
	return 0;
}

void Integrator(void (*gen)(double, double, double *, double *),
	 double t, double dt,double *x, double *xout)
{
	/*
	積分器の切り替えのWapper
	用意している積分器 genは
	Euler, RK2, RK4, symplectic, symplectic2
	です．
	*/
	gen(t, dt, x, xout);
}

void Euler(double t, double dt, double *x, double *xout)
{
	xout[0] = x[0] + dt * Diff_q(t,x[0],x[1]);
	xout[1] = x[1] + dt * Diff_p(t,x[0],x[1]);
}

void RK2(double t, double dt, double *x, double *xout)
{
	int i;
	double k1[2],k2[2];
	k1[0] = dt * Diff_q(t,x[0],x[1]);
	k1[1] = dt * Diff_p(t,x[0],x[1]);
	k2[0] = dt * Diff_q(t + dt/2.0, x[0] + k1[0]/2.0, x[1] + k1[1]/2.0);
	k2[1] = dt * Diff_p(t + dt/2.0, x[0] + k1[0]/2.0, x[1] + k1[1]/2.0);
	for(i=0;i<2;i++){
		xout[i] = x[i] + k2[i];
	}
}
void RK4(double t, double dt, double *x, double *xout)
{
	int i;
	double k1[2],k2[2],k3[2],k4[2];
	k1[0] = dt * Diff_q(t,x[0],x[1]);
	k1[1] = dt * Diff_p(t,x[0],x[1]);
	
	k2[0] = dt * Diff_q(t + dt/2.0, x[0] + k1[0]/2.0, x[1] + k1[1]/2.0);
	k2[1] = dt * Diff_p(t + dt/2.0, x[0] + k1[0]/2.0, x[1] + k1[1]/2.0);
	
	k3[0] = dt * Diff_q(t + dt/2.0, x[0] + k2[0]/2.0, x[1] + k2[1]/2.0);
	k3[1] = dt * Diff_p(t + dt/2.0, x[0] + k2[0]/2.0, x[1] + k2[1]/2.0);
	
	k4[0] = dt * Diff_q(t + dt, x[0] + k3[0], x[1] + k3[1]);
	k4[1] = dt * Diff_p(t + dt, x[0] + k3[0], x[1] + k3[1]);
	
	for(i=0;i<2;i++){
		xout[i] = x[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
	}
}

void Symplectic(double t, double dt, double *x, double *xout)
{
	xout[0] = x[0] + dt * Diff_q(t,   x[0], x[1]);
	xout[1] = x[1] + dt * Diff_p(t,xout[0], x[1]);
}

void Symplectic2(double t, double dt, double *x, double *xout)
{
	double dummy;
	dummy   = x[0]   + dt/2.0 * Diff_q(t,   x[0],    x[1]);
	xout[1] = x[1]   + dt     * Diff_p(t,  dummy,    x[1]);
	xout[0] = dummy  + dt/2.0 * Diff_q(t,  dummy, xout[1]);
}

