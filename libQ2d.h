#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define PRINT

#define D0 3./2.
#define D1 -2.
#define D2 1./2.
#define I0 5./12.
#define I1 8./12.
#define I2 -1./12.

struct pmct{
	int Nt,Ntexp,Nt2,Nc,rpt,den,itr;	//Nt: #(grid elements), Ntexp: log2(Nt), rpt, Nt2: Nt/2, Nsh: #(short-grid elements), den: Nt/Nsh , itr: #(cycles)
	double t,t0,eps,alpha;
	double p,s_eps,s,T,beta;			//p: (p+s)-spin, s_eps: weight of s-interaction, s: (p+s)-spin, T: temperature, beta: inverse temperature of the initial equilibrium
};

struct parr{							//see appendix C
	double dt,dmu;
	double *mu,*E;
	double **C,**Q,**f1,**f2;			//X can stend for R,F or Q
	double **dCh,**dQh,**df1h,**df2h;
	double **dCv,**dQv,**df1v,**df2v;
};

struct psys{
	char file[100],dir[100];			//file and directory name
	double *sh,*sh2;
};


void parameters_initialization(struct pmct *z,struct parr *x, struct psys *w, int argc, char *argv[]);

double power(double x,int p);
double f(double x, struct pmct *z);
double fd1(double x, struct pmct *z);
double fd2(double x, struct pmct *z);

void array_initialization(struct pmct *z,struct parr *x);
void contract(struct pmct *z,struct parr *x,double *dt,double *dmu);
void extrapolate(struct pmct *z,struct parr *x,int i);
void extrapolate_step(struct pmct *z,struct parr *x,int i, int j);

void initialarray(struct pmct *z,struct parr *x);
void write_C_0(struct pmct *z, struct parr *x, struct psys *w, int ini, int ifi);
void write_C(struct pmct *z, struct parr *x, struct psys *w, int j);
void final_write_C(struct pmct *z, struct parr *x, struct psys *w);

void snap_config(double **C,double time,struct pmct *z,struct psys *w);
void snap_vector(double **C, double **dC, char *name, double time, struct pmct *z, struct psys *w);

/**************************PaRaMeTeRs****************************/

void parameters_initialization(struct pmct *z,struct parr *x, struct psys *w, int argc, char *argv[]) {
	
	#ifdef N_PRINT
		printf("\nWithout printing!\nTo change see #define N_PRINT\n");
	#endif

	if(argc!=5) { printf("Input 4 parameters! $./FP <eps_s> <T> <Te> <log2(Nt)>\n"); exit(0); }

	/* system parameters */
	z->T = atof(argv[2]);		//  T
	z->beta = 1./atof(argv[3]);		//  beta'
	z->p     = 3.0;			//  f1(x) = x^p/2
	z->s_eps = atof(argv[1]);		//	weight of s-spin
	z->s     = 4.0;			//  f2(x) = x^s/2

	/* mct parameters */
	z->Ntexp = atoi(argv[4]);         	//  Nt=2^{Ntexp}  [5..10]
	z->Nt = 1 << z->Ntexp;
	z->Nt2= (int)z->Nt/2;

	z->den  = 1 << 5;			// Nt/Nc
	z->Nc = 32; //<< (z->Ntexp-5);
	z->itr = 34;				// number of cycle.

	z->t0 = 1.0E-5/z->T;			// Initial time window
	z->eps= 1.0E-10; 			// Accepted error distance between the solution and the discrete equation ->
	z->rpt= 100;				// 		in the iteration procedure with z->rpt the maximum number of iterations

	z->alpha = 1.;				// Coefficient of self-consistence iteration

	/* array parameters */
	x->dt = z->t0/z->Nt2;		// Grid time set to Initial Grid time
	x->dmu = 0.;				//IMPORTANT: local upgrade of mu


	/* shared memory */
	/*#ifndef N_PRINT
		int KEY = z->Ntexp*12;							//  for the printing procedure
		w->sh = start_double_shared_memory(KEY,z->Nt*z->Nt);			//PRINT_INITIALIZE
		w->sh2 = start_double_shared_memory(KEY+1,z->Nt*z->Nt);		//PRINT_INITIALIZE
	#endif*/

	/* output files */
	sprintf(w->file,"0.dat");
	sprintf(w->dir,"Q_eps%.2f_T%.5f_Tp%.5fNt%d",z->s_eps,z->T,1./z->beta,z->Nt);
    mkdir(w->dir, 0700);

    printf("\nDIRECTORY_OUTPUT: %s\n", w->dir);

	printf("\n-----------------------------------------------------START-------------------------------------------------------\n");
	printf("------SYSTEM: (%d+eps*%d)-spin glass, eps = %2.3f------------------------------------------------------------------\n",(int)z->p,(int)z->s,z->s_eps);
	printf("------PARAMETERS: T = %1.3f, T' = %1.3f, grid dimension = %d (Nsh = %d), initial time window = %.2e----------\n",z->T,1./z->beta,z->Nt,z->Nc,z->t0/z->T);
	printf("-----------------------------------------------------------------------------------------------------------------\n\n");

}

void write_parameters(struct pmct *z,struct parr *x, struct psys *w){
  FILE *fout;
  char fn[100];
  sprintf(fn,"%s/par.txt",w->dir);
  if((fout=fopen(fn, "w"))==NULL){
    fprintf(stderr,"write_parameters: Cannot open a outfile\n");
    exit(1);
  }

  fprintf(fout,"####[Output]##############################");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# Program name:   '%s'\n",thisfile);
  fprintf(fout,"# This directory name: '%s'\n",w->dir);
  fprintf(fout,"####[System-related parameters]############");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# z->T      =%.5f\n", z->T);
  fprintf(fout,"# z->Tp     =%.5f\n", 1./z->beta);
  fprintf(fout,"# z->Nt     =2^%d=%d\n",   z->Ntexp,z->Nt);
  fprintf(fout,"# z->Nc     =%d   \n",   z->Nc);
  fprintf(fout,"# z->p      =%.0f  \n",   z->p);
  fprintf(fout,"# z->s_eps  =%.2e  \n",   z->s_eps);
  fprintf(fout,"# z->s      =%.0f  \n",   z->s);
  fprintf(fout,"# z->itr    =%d    \n",   z->itr);
  fprintf(fout,"# z->t0     =%.2e \n",   z->t0);
  fprintf(fout,"# z->rpt    =%d   \n",   z->rpt);
  fprintf(fout,"# z->eps    =%.1e \n",   z->eps);
  fprintf(fout,"####[Ouput format]##########################");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# 0.dat:\n"); 
  fprintf(fout,"# 1\t2\t3\t4\t5\t6\t7\t8\t9");
  fprintf(fout,"# t\tC(t,0)\tdCv(t,0)\tdCh(t,0)\tQ(t,0)\tdQv(t,0)\tdQh(t,0)\tmu(t)\tE(t)");
  fprintf(fout,"# x.dat:\n"); 
  fprintf(fout,"# 1\t2\t3\t4\t5\t6\t7");
  fprintf(fout,"# t-t'\tC(t,t')\tdCv(t,t')\tdCh(t,t')\tQ(t,t')\tdQv(t,t')\tdQh(t,t')");
  fprintf(fout,"############################################");
  fprintf(fout,"##########################################\n");
  fclose(fout);
}

double power(double x,int p){
	double ff=1.; int i;
	for(i=0;i<p;i++) { ff *= x; }
	return ff;
}

double f(double x, struct pmct *z){
	return 0.5*power(x,z->p) + z->s_eps*0.5*power(x,z->s);
}

double fd1(double x, struct pmct *z){
	return 0.5*(z->p)*power(x,z->p-1.0) + z->s_eps*0.5*(z->s)*power(x,z->s-1.0);
}

double fd2(double x, struct pmct *z){
	return 0.5*(z->p)*(z->p-1.0)*power(x,z->p-2.0) + z->s_eps*0.5*(z->s)*(z->s-1.0)*power(x,z->s-2.0);
}

void array_initialization(struct pmct *z,struct parr *x) {
	/*------------------------------array initialization------------------------------*/
	int i;

	x->mu = (double *) calloc ((z->Nt+1),sizeof(double));
	x->E = (double *) calloc ((z->Nt+1),sizeof(double));
	x->C  = (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->Q  = (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->f1  = (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->f2  = (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->dCh= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->dQh= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->df1h= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->df2h= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->dCv= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->dQv= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->df1v= (double **) malloc ((z->Nt+1) * sizeof(double*));
	x->df2v= (double **) malloc ((z->Nt+1) * sizeof(double*));

	for(i=0;i<(z->Nt+1);i++) {
		x->C[i]  = (double *) calloc ((z->Nt+1),sizeof(double));
		x->Q[i]  = (double *) calloc ((z->Nt+1),sizeof(double));
		x->f1[i]  = (double *) calloc ((z->Nt+1),sizeof(double));
		x->f2[i]  = (double *) calloc ((z->Nt+1),sizeof(double));
		x->dCh[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->dQh[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->df1h[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->df2h[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->dCv[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->dQv[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->df1v[i]= (double *) calloc ((z->Nt+1),sizeof(double));
		x->df2v[i]= (double *) calloc ((z->Nt+1),sizeof(double));

	}
/*------------------------------end array initialization------------------------------*/
}

void initialarray(struct pmct *z,struct parr *x){
	int i,j;

	for(i=0;i<=z->Nt2;i++){
		x->mu[i] = 0.0;
		for(j=0;j<=i;j++){
			x->C[i][j]= 1.0 - (double)(i-j)*x->dt*z->T;	// very short time expansion --> time homogeneous initial C
			/*-(1-z->beta*z->T)*fd1(1,z)			ADDED LATER	*/							
			x->Q[i][j]= 0.0;														// very short time expansion --> FDT initially respected (Q=0)
			x->f1[i][j]= fd1(x->C[i][j],z);
			x->f2[i][j]= fd2(x->C[i][j],z);
		}
	}
	for(i=1;i<=z->Nt2;i++){
		for(j=0;j<i;j++){
			x->dCh[i][j]= 0.5*(x->C[i-1][j]+x->C[i][j]);
			x->dQh[i][j]= 0.5*(x->Q[i-1][j]+x->Q[i][j]);
			x->df1h[i][j]= 0.5*(x->f1[i-1][j]+x->f1[i][j]);
			x->df2h[i][j]= 0.5*(x->f2[i-1][j]+x->f2[i][j]);
			x->dCv[i][j]= 0.5*(x->C[i][j+1]+x->C[i][j]);
			x->dQv[i][j]= 0.5*(x->Q[i][j+1]+x->Q[i][j]);
			x->df1v[i][j]= 0.5*(x->f1[i][j+1]+x->f1[i][j]);
			x->df2v[i][j]= 0.5*(x->f2[i][j+1]+x->f2[i][j]);
		}
	}

	printf("\t (i)\t(err):(err_lmt)\t(rpt)\n");
}

void contract(struct pmct *z,struct parr *x,double *dt,double *dmu){
  int i,j;
  double Dl;
  i=z->Nt;
  for(j=z->Nt-z->Nc*2+1;j<=z->Nt-z->Nc;j++){
    Dl=(x->Q[i][j]-x->Q[i][j-1])*(I2*(x->f1[i][j+1]+x->f2[i][j+1]*x->C[i][j+1])
				+I1*(x->f1[i][j  ]+x->f2[i][j  ]*x->C[i][j  ])
				+I0*(x->f1[i][j-1]+x->f2[i][j-1]*x->C[i][j-1]) );
    (*dmu) += Dl;
  }
  for(i=1;i<=z->Nt2;i++){
    for(j=0;j<=i-1;j++){
      x->dCh[i][j]= 0.5*(x->dCh[2*i][2*j]+x->dCh[2*i-1][2*j]);
      x->dQh[i][j]= 0.5*(x->dQh[2*i][2*j]+x->dQh[2*i-1][2*j]);
      x->df1h[i][j]= 0.5*(x->df1h[2*i][2*j]+x->df1h[2*i-1][2*j]);
      x->df2h[i][j]= 0.5*(x->df2h[2*i][2*j]+x->df2h[2*i-1][2*j]);
      x->dCv[i][j]= 0.5*(x->dCv[2*i][2*j+1]+x->dCv[2*i][2*j]);
      x->dQv[i][j]= 0.5*(x->dQv[2*i][2*j+1]+x->dQv[2*i][2*j]);
      x->df1v[i][j]= 0.5*(x->df1v[2*i][2*j+1]+x->df1v[2*i][2*j]);
      x->df2v[i][j]= 0.5*(x->df2v[2*i][2*j+1]+x->df2v[2*i][2*j]);
    }
  }
  for(i=0;i<=z->Nt2;i++){
    for(j=0;j<=i;j++){
      x->C[i][j]= x->C[2*i][2*j];
      x->Q[i][j]= x->Q[2*i][2*j];
      x->f1[i][j]= x->f1[2*i][2*j];
      x->f2[i][j]= x->f2[2*i][2*j];
    }
  }
  (*dt) *= 2.0;
}

void extrapolate(struct pmct *z,struct parr *x,int i) {

	int j;
	// (1) copy the value for the top Nsh from the previous column
	for(j=i-z->Nc;j<=i;j++){
		x->C[i][j]= x->C[i-1][j-1];
		x->Q[i][j]= x->Q[i-1][j-1];
		x->f1[i][j]= x->f1[i-1][j-1];
		x->f2[i][j]= x->f2[i-1][j-1];
	}
	for(j=i-z->Nc;j<i;j++){
		x->dCh[i][j]= x->dCh[i-1][j-1];
		x->dQh[i][j]= x->dQh[i-1][j-1];
		x->df1h[i][j]= x->df1h[i-1][j-1];
		x->df2h[i][j]= x->df2h[i-1][j-1];
		x->dCv[i][j]= x->dCv[i-1][j-1];
		x->dQv[i][j]= x->dQv[i-1][j-1];
		x->df1v[i][j]= x->df1v[i-1][j-1];
		x->df2v[i][j]= x->df2v[i-1][j-1];
	}

  // (2) Prepare test values
	for(j=0;j<=i-z->Nc-1;j++){
		x->C[i][j]= x->C[i-1][j];
		x->Q[i][j]= x->Q[i-1][j];
		x->f1[i][j]= fd1(x->C[i][j],z);
		x->f2[i][j]= fd2(x->C[i][j],z);
	}
	for(j=0;j<=i-z->Nc-1;j++){
		x->dCh[i][j] = I2*x->C[i-2][j]+I1*x->C[i-1][j]+I0*x->C[i][j];
		x->dQh[i][j] = I2*x->Q[i-2][j]+I1*x->Q[i-1][j]+I0*x->Q[i][j];
		x->df1h[i][j]= I2*x->f1[i-2][j]+I1*x->f1[i-1][j]+I0*x->f1[i][j];
		x->df2h[i][j]= I2*x->f2[i-2][j]+I1*x->f2[i-1][j]+I0*x->f2[i][j];
		x->dCv[i][j] = I2*x->C[i][j+2]+I1*x->C[i][j+1]+I0*x->C[i][j];
		x->dQv[i][j] = I2*x->Q[i][j+2]+I1*x->Q[i][j+1]+I0*x->Q[i][j];
		x->df1v[i][j]= I2*x->f1[i][j+2]+I1*x->f1[i][j+1]+I0*x->f1[i][j];
		x->df2v[i][j]= I2*x->f2[i][j+2]+I1*x->f2[i][j+1]+I0*x->f2[i][j];
	}
}

void extrapolate_step(struct pmct *z,struct parr *x,int i, int j) {

	x->dCh[i][j]= I2*x->C[i-2][j]+I1*x->C[i-1][j]+I0*x->C[i][j];
	x->dQh[i][j]= I2*x->Q[i-2][j]+I1*x->Q[i-1][j]+I0*x->Q[i][j];
	x->df1h[i][j]= I2*x->f1[i-2][j]+I1*x->f1[i-1][j]+I0*x->f1[i][j];
	x->df2h[i][j]= I2*x->f2[i-2][j]+I1*x->f2[i-1][j]+I0*x->f2[i][j];
	x->dCv[i][j]= I2*x->C[i][j+2]+I1*x->C[i][j+1]+I0*x->C[i][j];
	x->dQv[i][j]= I2*x->Q[i][j+2]+I1*x->Q[i][j+1]+I0*x->Q[i][j];
	x->df1v[i][j]= I2*x->f1[i][j+2]+I1*x->f1[i][j+1]+I0*x->f1[i][j];
	x->df2v[i][j]= I2*x->f2[i][j+2]+I1*x->f2[i][j+1]+I0*x->f2[i][j];
}

void write_C_0(struct pmct *z, struct parr *x, struct psys *w, int ini, int ifi){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%s",w->dir,w->file);
	if((fout=fopen(fn, "at"))==NULL){
		fprintf(stderr,"write_C_0: Cannot open a outfile\n");
		exit(1);
	}

	for(i=ini;i<=ifi;i++) {
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.10e\t%1.10e\n",
			(double)i*x->dt*z->T,x->C[i][0],x->dCv[i][0],x->dCh[i][0],x->Q[i][0],x->dQv[i][0],x->dQh[i][0],x->mu[i],x->E[i]);
  	}
	fclose(fout);
}

void write_C(struct pmct *z, struct parr *x, struct psys *w, int j){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%.2e.dat",w->dir,j*x->dt*z->T);
	if((fout=fopen(fn, "w"))==NULL){
		fprintf(stderr," write_C: Cannot open a outfile\n");
		exit(1);
	}

	for(i=j+1;i<=z->Nt;i++) {
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",(double)(i-j)*x->dt*z->T,x->C[i][j],x->dCv[i][j],x->dCh[i][j],x->Q[i][j],x->dQv[i][j],x->dQh[i][j]);
  	}
	fclose(fout);
}

void final_write_C(struct pmct *z, struct parr *x, struct psys *w){
	int i,j=2;
	for(i=1;i<z->Ntexp;i++) {
		write_C(z,x,w,j); j = power(2,i); 
	} 
}

//#include <algorithm>

/**************************PrInT****************************/

#ifndef PRINT
#include "lib/share_memory.h"

void print_vector(double *sh, double **R, int l) {

	int j,k;
	double temp_sh;

	temp_sh = *(sh-2);
	*(sh-2) = 0.;
	for(k=0;k<l;k++) {
		for(j=0;j<l;j++) { *(sh+k*l+j) = R[k][j]; }
	}
	*(sh-2) = temp_sh+1.;
}
#endif

/**********************************************************/

void snap_config(double **C,double time,struct pmct *z,struct psys *w) {
	int i,j;
	FILE *f;
	char fn[100];
	sprintf(fn,"%s/config_time%.2e.txt",w->dir,time);
	if((f=fopen(fn, "w"))==NULL){
		fprintf(stderr,"save_config: Cannot open a outfile\n");
		exit(1);
	}

	int di = z->Nt/256;

	for (i=0; i<z->Nt+1; i+=di) {
		for (j=0; j<z->Nt; j+=di) {
			fprintf(f,"%.5f,",C[i][j]);
		}
		fprintf(f,"%.5f\n",C[i][j]);
	}
}

void snap_vector(double **C, double **dC, char *name, double time, struct pmct *z, struct psys *w) {
	int k;
	FILE *f;
	char fn[100];
	sprintf(fn,"%s/%s%.2e.txt",w->dir,name,time);
	printf("PRINTTT");
	if((f=fopen(fn, "w"))==NULL){
		fprintf(stderr,"save_config: Cannot open a outfile\n");
		exit(1);
	}

	for (k=0; k<z->Nt+1; k++) {
		fprintf(f,"%.5f %.5f\n",C[z->Nt][k],dC[z->Nt][k]);
	}
	fclose(f);
}

