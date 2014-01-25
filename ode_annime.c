#include <stdio.h>
#include <unistd.h>

void Integrator(void (*gen)(double, double, double *, double *),
	 double t, double dt,double *x, double *xout);
void Euler(double t, double dt, double *x, double *xout);
void RK2(double t, double dt, double *x, double *xout);
void RK4(double t, double dt, double *x, double *xout);
void Symplectic(double t, double dt, double *x, double *xout);
void Symplectic2(double t, double dt, double *x, double *xout);

void SetInit(int num, double *x);
void GnuplotInit(FILE *gp, int num, double xmin, double xmax,double ymin, double ymax);
	
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
	return 0.5*p*p + 0.5*q*q;
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
//	return -4*q*(q*q-10.0);
}

/*---- ユーザー定義関数 (ここまで)：end user defined functions ---- */

int main(void)
{
	FILE *f,*ff;
	FILE *gp;
	gp = popen("gnuplot","w");
	char fname[100];

	int i,j;
	int sample,iteration;
	iteration=2000;	
	sample=10;

	double x[2], xx[2];
	double points[2*sample];
	double t, dt;

	dt = 0.01;


	SetInit(sample,points);
	
	f = fopen("trajectories.dat","w");
	ff= fopen("energies.dat","w");
	
	for (i=0;i<iteration;i++){
		t = i*dt;
		fprintf(ff,"%lf ", t);
		
		GnuplotInit(gp, sample, -10,10,-10,10);

		for (j=0;j<sample;j++){
			x[0] = points[2*j];
			x[1] = points[2*j+1];
						
			fprintf(gp,"%lf %lf\n", x[0],x[1]);
			fprintf(gp,"e\n");
			fprintf(f,"%.16e %.16e ", x[0], x[1]);
			fprintf(ff,"%lf ", Hamiltonian(t, x[0],x[1]) );			
			
			Integrator(RK4, t,dt,x,xx);

			points[2*j] = xx[0];
			points[2*j+1] = xx[1];

		}
		
		fprintf(f,"\n");
		fprintf(ff,"\n");		
		
		fflush(f);				
		fflush(ff);						
		fflush(gp);		
		
		if (i==0){
		    printf("Start : Push enter key! \n");
		    printf("Stop  : Push Control key + c \n");			
		    while(getchar() != '\n'); 
		}		

	}
	fclose(gp);
	fclose(f);
	fclose(ff);
	return 0;
}

void SetInit(int num, double *x)
{
	int i;
	for(i=0;i<2*num;i++){
		if (i%2==0)
			x[i] = 0.5*i;
		else
			x[i] = 0.0;
	}
}



void GnuplotInit(FILE *gp, int num, double xmin, double xmax,double ymin, double ymax)
{
	int i;
	fprintf(gp, "set grid\n");
	fprintf(gp, "set xlabel 'q' font 'Helvetica,20'\n");
	fprintf(gp, "set ylabel 'p' font 'Helvetica,20'\n");
	fprintf(gp, "set size square\n");
	fprintf(gp, "set xrange [%lf: %lf]\n",xmin,xmax);
	fprintf(gp, "set yrange [%lf: %lf]\n",ymin,ymax);
	fprintf(gp, "unset key\n");
	if (num ==1){
		fprintf(gp,"plot '-' with points pt 7 ps 4 lc -1 \n");
	}
	else{
		fprintf(gp,"plot '-' with points pt 7 ps 4 lc -1 ,");
		i=1;
		for(i=1;i<num-1;i++){
			fprintf(gp," '-' with points pt 7 ps 4 lc %d ,",i);
		}
		fprintf(gp," '-' with points pt 7 ps 4 lc %d \n", i);	
	}
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

