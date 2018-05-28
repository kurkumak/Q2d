#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "segments.h"

#define PRINT

#ifndef I1

#define D1 3./2.
#define D2 -2.
#define D3 1./2.
#define I1 5./12.
#define I2 8./12.
#define I3 -1./12.

#endif

struct pmct{
	int Nt,Ntexp,Nt2,Nc,rpt,den,itr;	//Nt: #(grid elements), Ntexp: log2(Nt), rpt, Nt2: Nt/2, Nsh: #(short-grid elements), den: Nt/Nsh , itr: #(cycles)
	double t,t0,eps,alpha;
	double p,s_eps,s,T,beta;			//p: (p+s)-spin, s_eps: weight of s-interaction, s: (p+s)-spin, T: temperature, beta: inverse temperature of the initial equilibrium
	int *scmaxx;
	double f_err;	
};

struct parr{							//see appendix C
	double dt,dmu;
	struct segment *mu,*E;
	struct corner **C,**Q,**f1,**f2;			//X can stend for R,F or Q
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
void write_C_0(struct parr *x, struct psys *w, int ini, int ifi);
void write_C(struct pmct *z, struct parr *x, struct psys *w, int j);

void snap_config(double **C,double time,struct pmct *z,struct psys *w);
void snap_vector(double **C, double **dC, char *name, double time, struct pmct *z, struct psys *w);

void print_C_screen(int i,struct parr *x);

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
	z->Nc = 3; //<< (z->Ntexp-5);
	z->itr = 34;				// number of cycle.

	z->t0 = 1.0E-4;			// Initial time window
	z->eps= 1.0E-10; 			// Accepted error distance between the solution and the discrete equation ->
	z->rpt= 1000;				// 		in the iteration procedure with z->rpt the maximum number of iterations

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
	sprintf(w->dir,"QAd_eps%.2f_T%.5f_Tp%.5fNt%d",z->s_eps,z->T,1./z->beta,z->Nt);
    mkdir(w->dir, 0700);

    printf("\nDIRECTORY_OUTPUT: %s\n", w->dir);

	printf("\n-----------------------------------------------------START-------------------------------------------------------\n");
	printf("------SYSTEM: (%d+eps*%d)-spin glass, eps = %2.3f------------------------------------------------------------------\n",(int)z->p,(int)z->s,z->s_eps);
	printf("------PARAMETERS: T = %1.3f, T' = %1.3f, grid dimension = %d (Nsh = %d), initial time window = %.2e----------\n",z->T,1./z->beta,z->Nt,z->Nc,z->t0);
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
	return 0.5*power(x,z->p)/z->T + z->s_eps*0.5*power(x,z->s)/z->T;
}

double fd1(double x, struct pmct *z){
	return 0.5*(z->p)*power(x,z->p-1.0)/z->T + z->s_eps*0.5*(z->s)*power(x,z->s-1.0)/z->T;
}

double fd2(double x, struct pmct *z){
	return 0.5*(z->p)*(z->p-1.0)*power(x,z->p-2.0)/z->T + z->s_eps*0.5*(z->s)*(z->s-1.0)*power(x,z->s-2.0)/z->T;
}

void array_initialization(struct pmct *z,struct parr *x) {
	/*------------------------------array initialization------------------------------*/
	int i,j,itr;

	x->mu = (struct segment *) calloc ((z->Nt+1),sizeof(struct segment));
	x->E = (struct segment *) calloc ((z->Nt+1),sizeof(struct segment));
	x->C  = (struct corner **) malloc ((z->Nt+1) * sizeof(struct corner *));
	x->Q  = (struct corner **) malloc ((z->Nt+1) * sizeof(struct corner *));
	x->f1  = (struct corner **) malloc ((z->Nt+1) * sizeof(struct corner *));
	x->f2  = (struct corner **) malloc ((z->Nt+1) * sizeof(struct corner *));
	
	for(i=0;i<=z->Nt;i++) {
		x->C[i]  = (struct corner *) calloc (i,sizeof(struct corner));
		x->Q[i]  = (struct corner *) calloc (i,sizeof(struct corner));
		x->f1[i]  = (struct corner *) calloc (i,sizeof(struct corner));
		x->f2[i]  = (struct corner *) calloc (i,sizeof(struct corner));
	}

	for(i=0;i<=z->Nt;i++) {
		for(j=0;j<i;j++) {
			x->C[i][j].h  = (struct segment *) calloc(1,sizeof(struct segment));
			x->Q[i][j].h  = (struct segment *) calloc(1,sizeof(struct segment));
			x->f1[i][j].h  = (struct segment *) calloc(1,sizeof(struct segment));
			x->f2[i][j].h  = (struct segment *) calloc(1,sizeof(struct segment));
			x->C[i][j].v  = (struct segment *) calloc(1,sizeof(struct segment));
			x->Q[i][j].v  = (struct segment *) calloc(1,sizeof(struct segment));
			x->f1[i][j].v  = (struct segment *) calloc(1,sizeof(struct segment));
			x->f2[i][j].v  = (struct segment *) calloc(1,sizeof(struct segment));
		}
	}

	z->scmaxx = (int *) malloc ((z->Nt+1) * sizeof(int));
/*------------------------------end array initialization------------------------------*/
}

void initialarray(struct pmct *z,struct parr *x){
	int i,j;
	double dT=x->dt*z->T;

	for(i=1;i<=z->Nt2;i++){
		fill_segment(&(x->mu[i]),0.,0.,0.);
		/*x->C[i][i-1].h->x1 = 1.;
		x->C[i][i-1].v->x2 = 1.;
		x->Q[i][i-1].h->x1 = 0.;
		x->Q[i][i-1].v->x2 = 0.;				// very short time expansion --> FDT initially respected (Q=0)
		x->f1[i][i-1].h->x1 = fd1(1.,z);
		x->f1[i][i-1].v->x2 = fd1(1.,z);
		x->f2[i][i-1].h->x1 = fd2(1.,z);
		x->f2[i][i-1].v->x2 = fd2(1.,z);*/
	
		for(j=0;j<i;j++){
			fill_corner(&(x->C[i][j]),1.0-(i-j)*dT
									,1.0-(i-1-j)*dT
									,1.0-(i-j-1)*dT
									,1.0-(i-0.5-j)*dT
									,1.0-(i-j-0.5)*dT);
			fill_corner(&(x->Q[i][j]),0.,0.,0.,0.,0.);

			fill_corner(&(x->f1[i][j])
									,fd1(1.0-(i-j)*dT,z)
									,fd1(1.0-(i-1-j)*dT,z)
									,fd1(1.0-(i-j-1)*dT,z)
									,0.5*(fd1(1.0-(i-j)*dT,z)+fd1(1.0-(i-1-j)*dT,z))
									,0.5*(fd1(1.0-(i-j)*dT,z)+fd1(1.0-(i-1-j)*dT,z)));

			fill_corner(&(x->f2[i][j])
									,fd2(1.0-(i-j)*dT,z)
									,fd2(1.0-(i-1-j)*dT,z)
									,fd2(1.0-(i-j-1)*dT,z)
									,0.5*(fd2(1.0-(i-j)*dT,z)+fd2(1.0-(i-1-j)*dT,z))
									,0.5*(fd2(1.0-(i-j)*dT,z)+fd2(1.0-(i-1-j)*dT,z)));
			/*x->C[i][j].h->x2 = 1.0 - (double)(i-j)*x->dt*z->T;	// very short time expansion --> time homogeneous initial C
			x->C[i][j].h->x_ = 
			x->C[i][j].v->x1 = x->C[i][j].h->x2;
			x->C[i][j].v->x_ = 
			/*-(1-z->beta*z->T)*fd1(1,z)			ADDED LATER	*/							
			/*x->Q[i][j].h->x2 = 0.0;
			x->Q[i][j].v->x1 = x->Q[i][j].h->x2;				// very short time expansion --> FDT initially respected (Q=0)
			x->Q[i][j].h->x_ = 0.0;

			x->f1[i][j].h->x2 = fd1(x->C[i][j],z);
			x->f1[i][j].v->x1 = x->f1[i][j].h->x2;
			x->f2[i][j].h->x2 = fd2(x->C[i][j],z);
			x->f2[i][j].v->x1 = x->f2[i][j].h->x2;*/
		}
				//print_C_screen(i,x);

	}

	printf("\t (i)\t(err):(err_lmt)\t(rpt)\n");
}

void contract(struct pmct *z,struct parr *x,double *dt,double *dmu){
	int i,j;
	double Dl;
	i=z->Nt;
	for(j=z->Nt-z->Nc*2+1;j<=z->Nt-z->Nc;j++){

		Dl =(x->Q[i][j].v->x1-x->Q[i][j-1].v->x1)*
			(I3*(x->f1[i][j+1].v->x1+x->f2[i][j+1].v->x1*x->C[i][j+1].v->x1)
			+I2*(x->f1[i][j  ].v->x1+x->f2[i][j  ].v->x1*x->C[i][j  ].v->x1)
			+I1*(x->f1[i][j-1].v->x1+x->f2[i][j-1].v->x1*x->C[i][j-1].v->x1));
		//Dl = product_d(x->Q[i][j-1].v,x->f1[i][j-1].v)+product3_d(x->Q[i][j-1].v,x->f2[i][j-1].v,x->C[i][j-1].v);
		    	/*Dl=(x->Q[i][j]-x->Q[i][j-1])*(I3*(x->f1[i][j+1]+x->f2[i][j+1]*x->C[i][j+1])
									+I2*(x->f1[i][j  ]+x->f2[i][j  ]*x->C[i][j  ])
									+I1*(x->f1[i][j-1]+x->f2[i][j-1]*x->C[i][j-1]) );*/
    	(*dmu) += Dl;
  	}

  	for(i=1;i<=z->Nt2;i++) {

  		//printf("\n\nPRIMA//\n");
  		//printf("C[%d][%d]=%f",i,0,x->C[i][0].v->x1);

  		//print_C_screen(i,x);

		for(j=0;j<i;j++) {
			contract_corners(&(x->C[2*i-1][2*j]), &(x->C[2*i][2*j+1]), &(x->C[2*i][2*j]), &(x->C[i][j]));
			contract_corners(&(x->Q[2*i-1][2*j]), &(x->Q[2*i][2*j+1]), &(x->Q[2*i][2*j]), &(x->Q[i][j]));
			contract_corners(&(x->f1[2*i-1][2*j]), &(x->f1[2*i][2*j+1]), &(x->f1[2*i][2*j]), &(x->f1[i][j]));
			contract_corners(&(x->f2[2*i-1][2*j]), &(x->f2[2*i][2*j+1]), &(x->f2[2*i][2*j]), &(x->f2[i][j]));
			//printf("%i %i",i,j); //getchar();
		}

	  	//printf("\n\nDOPO//\n");
  		//printf("C[%d][%d]=%f",i,0,x->C[i][0].v->x1); getchar();

		//print_C_screen(i,x);
	}

	for(i=z->Nt2+1;i<=z->Nt;i++) {
		for(j=0;j<i;j++) {
			empty_corner(&(x->C[i][j]));
			empty_corner(&(x->Q[i][j]));
			empty_corner(&(x->f1[i][j]));
			empty_corner(&(x->f2[i][j]));
		}
		//print_C_screen(i,x); fflush(stdout);

	}


  	/*for(i=1;i<=z->Nt2;i++){
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
  }*/
  (*dt) *= 2.0;
}

void extrapolate(struct pmct *z,struct parr *x,int i) {

	int j;
	// (1) copy the value for the top Nsh from the previous column
	for(j=i-z->Nc;j<i;j++){
		copy_corner(&(x->C[i-1][j-1]),&(x->C[i][j]));
		copy_corner(&(x->Q[i-1][j-1]),&(x->Q[i][j]));
		copy_corner(&(x->f1[i-1][j-1]),&(x->f1[i][j]));
		copy_corner(&(x->f2[i-1][j-1]),&(x->f2[i][j]));
	}
	/*
	for(j=i-z->Nc;j<i;j++){
		x->dCh[i][j]= x->dCh[i-1][j-1];
		x->dQh[i][j]= x->dQh[i-1][j-1];
		x->df1h[i][j]= x->df1h[i-1][j-1];
		x->df2h[i][j]= x->df2h[i-1][j-1];
		x->dCv[i][j]= x->dCv[i-1][j-1];
		x->dQv[i][j]= x->dQv[i-1][j-1];
		x->df1v[i][j]= x->df1v[i-1][j-1];
		x->df2v[i][j]= x->df2v[i-1][j-1];
	}*/

  // (2) Prepare test values
	for(j=i-z->Nc-1;j>=0;j--){
		extrapolate_corner(&(x->C[i-1][j]), &(x->C[i][j+1]), x->C[i-1][j].h->x2, &(x->C[i][j]));
		extrapolate_corner(&(x->Q[i-1][j]), &(x->Q[i][j+1]), x->Q[i-1][j].h->x2, &(x->Q[i][j]));
		extrapolate_corner(&(x->f1[i-1][j]), &(x->f1[i][j+1]), fd1(x->C[i][j].h->x1,z), &(x->f1[i][j]));
		extrapolate_corner(&(x->f2[i-1][j]), &(x->f2[i][j+1]), fd2(x->C[i][j].h->x1,z), &(x->f2[i][j]));
		/*x->C[i][j]= x->C[i-1][j];
		x->Q[i][j]= x->Q[i-1][j];
		x->f1[i][j]= fd1(x->C[i][j],z);
		x->f2[i][j]= fd2(x->C[i][j],z);*/
	}
			//print_C_screen(i,x);

	/*for(j=0;j<=i-z->Nc-1;j++){
		x->dCh[i][j] = I3*x->C[i-2][j]+I2*x->C[i-1][j]+I1*x->C[i][j];
		x->dQh[i][j] = I3*x->Q[i-2][j]+I2*x->Q[i-1][j]+I1*x->Q[i][j];
		x->df1h[i][j]= I3*x->f1[i-2][j]+I2*x->f1[i-1][j]+I1*x->f1[i][j];
		x->df2h[i][j]= I3*x->f2[i-2][j]+I2*x->f2[i-1][j]+I1*x->f2[i][j];
		x->dCv[i][j] = I3*x->C[i][j+2]+I2*x->C[i][j+1]+I1*x->C[i][j];
		x->dQv[i][j] = I3*x->Q[i][j+2]+I2*x->Q[i][j+1]+I1*x->Q[i][j];
		x->df1v[i][j]= I3*x->f1[i][j+2]+I2*x->f1[i][j+1]+I1*x->f1[i][j];
		x->df2v[i][j]= I3*x->f2[i][j+2]+I2*x->f2[i][j+1]+I1*x->f2[i][j];
	}*/
}

/*void extrapolate_step(struct pmct *z,struct parr *x,int i, int j) {

	x->dCh[i][j]= I3*x->C[i-2][j]+I2*x->C[i-1][j]+I1*x->C[i][j];
	x->dQh[i][j]= I3*x->Q[i-2][j]+I2*x->Q[i-1][j]+I1*x->Q[i][j];
	x->df1h[i][j]= I3*x->f1[i-2][j]+I2*x->f1[i-1][j]+I1*x->f1[i][j];
	x->df2h[i][j]= I3*x->f2[i-2][j]+I2*x->f2[i-1][j]+I1*x->f2[i][j];
	x->dCv[i][j]= I3*x->C[i][j+2]+I2*x->C[i][j+1]+I1*x->C[i][j];
	x->dQv[i][j]= I3*x->Q[i][j+2]+I2*x->Q[i][j+1]+I1*x->Q[i][j];
	x->df1v[i][j]= I3*x->f1[i][j+2]+I2*x->f1[i][j+1]+I1*x->f1[i][j];
	x->df2v[i][j]= I3*x->f2[i][j+2]+I2*x->f2[i][j+1]+I1*x->f2[i][j];
}*/

void write_C_0(struct parr *x, struct psys *w, int ini, int ifi){
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
			(double)i*x->dt,x->C[i][0].v->x1,x->C[i][0].v->x_,x->C[i][0].h->x_,
			x->Q[i][0].v->x1,x->Q[i][0].v->x_,x->Q[i][0].h->x_,x->mu[i].x1,x->E[i].x1);
  	}
	fclose(fout);
}

void write_C(struct pmct *z, struct parr *x, struct psys *w, int j){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%d_%.2e.dat",w->dir,j,j*x->dt);
	if((fout=fopen(fn, "w"))==NULL){
		fprintf(stderr," write_C: Cannot open a outfile\n");
		exit(1);
	}

	for(i=2;i<=z->Nt;i++) {
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",
			(double)(i-j)*x->dt,x->C[i][j].v->x1,x->C[i][j].v->x_,x->C[i][j].h->x_,
			x->Q[i][j].v->x1,x->Q[i][j].v->x_,x->Q[i][j].h->x_);
  	}
	fclose(fout);
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

	for (i=0; i<z->Nt+1; i++) {
		for (j=0; j<z->Nt; j++) {
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

void print_C_screen(int i,struct parr *x) {
	printf("line %d\n",i);
	int j;
	for(j=1;j<i;j++) {
		printf("C[%d]%.12f %.12f ",j,x->f1[i][j].v->x1,x->f1[i][j-1].v->x2);
	}
	printf("\n");	
	for(j=0;j<i;j++) {
		printf("v[%d]%.12f ",j,x->f1[i][j].v->x_);
	}
	printf("\n");
	for(j=0;j<i;j++) {
		printf("h[%d]%.12f ",j,x->f1[i][j].h->x_);
	}
	printf("\n");
	for(j=0;j<=i;j++) {
		printf("[%d]%lf ",j,x->mu[j].x1);
	}
	//getchar();
}

