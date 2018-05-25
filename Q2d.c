#define thisfile "Q2d.c"
#include "libQ2d.h"

/*****************************************************
p-spin model. Algorithm for chi_4. Similar to the aging (Kim&Latz, 2000)
written:  08/13/05
upgraded: 03/30/06
******************************************************
Reference article (appendix C): Spontaneous and induced dynamic correlations in glass formers. II. - L.Berthier, G.Biroli, G.-P.Bouchaud, W.Kob, K.Miyazaki, D.R.Reichman 
modified: 22/03/17
******************************************************/

#define Pi  3.141592653589793
#define Ntmax 1 << 20

int KEY;

void mct(struct pmct *z,struct parr *x,struct psys *w);
int step(int i,struct pmct *z,struct parr *x);

double SC1(double *gC,double *gQ,double D,double mu,struct pmct *z,struct parr *x,int i,int j);
double SC2(double *gC,double *gQ,double D,double mu,struct pmct *z,struct parr *x,int i,int j);
double I1C(struct parr *x,int i,int j,int m);
double I1Cs(struct parr *x,int i,int j,int m);
double I2C(struct parr *x,int i,int j,int m);
double I3C(struct parr *x,int i,int j);
double I4C(struct parr *x,int i,int j);
double I1Q(struct parr *x,int i,int j,int m);
double I1Qs(struct parr *x,int i,int j,int m);
double I2Q(struct parr *x,int i,int j,int m);

double mu_t(struct pmct *z,struct parr *x,int i);
double E_t(struct pmct *z,struct parr *x,int i);


int main(int argc, char *argv[]){

	struct pmct z;
	struct parr x;
	struct psys w;
	
	parameters_initialization(&z,&x,&w,argc,argv);

	write_parameters(&z,&x,&w);

	mct(&z,&x,&w);

	printf("\n-------------------------------------------------------END-------------------------------------------------------\n");

	return(0);
}

void mct(struct pmct *z,struct parr *x, struct psys *w){

	array_initialization(z,x);

	int i,rpt=0,itr=0;

/*------------------------------------------------------------*/
	initialarray(z,x); 	// prepare the array btwn 0 <= i,j <= Nt/2
/*------------------------------------------------------------*/

	write_C_0(x,w,0,z->Nt2);

	while(itr <= z->itr && rpt<z->rpt){
		
		i=z->Nt2+1;
		while(i<=z->Nt && rpt<z->rpt){
/*------------------------------------------------------------*/
			rpt = step(i,z,x);	// propagate the solution the array btwn Nt/2 <= i,j <= Nt
			i++;
/*------------------------------------------------------------*/
		}

		write_C_0(x,w,z->Nt2+1,z->Nt);
		write_C(z,x,w,1);

		if(itr%3==0) { snap_config(x->C,x->dt*z->Nt,z,w); }

/*------------------------------------------------------------*/
		//printf("\nmu_before = %f",mu_t(z,x,z->Nt));
		contract(z,x,&(x->dt),&(x->dmu));	// double the size of the system
		//printf("\nmu_after = %f\n",mu_t(z,x,z->Nt2));
/*------------------------------------------------------------*/
		printf(" %dth/%d cycle - window_time: %.2e \n",itr,z->itr,x->dt*z->Nt);
		itr++;
	}

	int j=2; for(i=1;i<z->Ntexp;i++) { write_C(z,x,w,j); j = power(2,i); } 
}

int step(int i,struct pmct *z,struct parr *x){

	int j,scmax;
	double D,err2,err2_temp;

	double * gC = (double *) calloc ((z->Nt),sizeof(double));
	double * gQ = (double *) calloc ((z->Nt),sizeof(double));

	extrapolate(z,x,i);

	// (3) Go to the SC (self-consistence) loop
	scmax = 0;
	err2 =1.0;

	while(err2 >= z->eps*z->eps && scmax < z->rpt){
		err2 = 0.0;

		//****** PART ----> j<i-1
		for(j=i-z->Nc-1;j>=0;j--){

			x->mu[i] = mu_t(z,x,i);
			D = D1/x->dt  + x->mu[i] + x->df1v[i][i-1];

			err2_temp = SC2(gC,gQ,D,x->mu[i],z,x,i,j);
			if (err2_temp>err2) { err2=err2_temp; }
			
			// renew all variable
			x->C[i][j] += gC[j];
			x->Q[i][j] += gQ[j];
			x->f1[i][j] = fd1(x->C[i][j],z);
			x->f2[i][j] = fd2(x->C[i][j],z);
		
			extrapolate_step(z,x,i,j);
		}

		scmax++;	
	}

	free(gC);
	free(gQ);

	x->E[i] = E_t(z,x,i);

	z->f_err = err2;

	printf("\r\t%4d\t%2.1e:%2.1e\t  %d\t",i,z->f_err,z->eps*z->eps,scmax); fflush(stdout);
	return scmax;
}

double SC2(double *gC, double *gQ, double D, double mu, struct pmct *z, struct parr *x,int i, int j){
  int m;
  double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q;
  m = (int)(0.5*(i+j));
  i1C = I1C(x,i,j,m);
  i2C = I2C(x,i,j,m);
  i3C = I3C(x,i,j);
  i4C = I4C(x,i,j);
  i1Q = I1Q(x,i,j,m);
  i2Q = I2Q(x,i,j,m);
  i3Q = i3C;
  i4Q = i4C;

  gC[j] = -D3/x->dt*x->C[i-2][j]-D2/x->dt*x->C[i-1][j]-i1C+i2C+i3C+i4C;
  gC[j]-= (1-z->T*z->beta)*fd1(x->C[i][0],z)*x->C[j][0];
  gC[j]/= D;
  gC[j]-= x->C[i][j];
  gQ[j] = -z->T+mu-D3/x->dt*x->Q[i-2][j]-D2/x->dt*x->Q[i-1][j]-i1Q-i2Q-i3Q-i4Q;
  gQ[j]+= (1-z->T*z->beta)*fd1(x->C[i][0],z)*x->C[j][0];
  gQ[j]/= D;
  gQ[j]-= x->Q[i][j];
	
  return gC[j]*gC[j]+gQ[j]*gQ[j];
}

double I1C(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x->f1[i][m]*x->C[m][j]-x->f1[i][j]*x->C[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x->df1v[i][l-1]*(x->C[l][j]-x->C[l-1][j]);
  }
  sum += x->df1v[i][i-1]*(-x->C[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x->f1[i][l]-x->f1[i][l-1])*x->dCh[l][j];
  }
  return sum;
}

double I2C(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x->df2v[i][l-1]*(x->Q[i][l]-x->Q[i][l-1])*(x->C[l][j]+x->C[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x->f2[i][l]+x->f2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCh[l][j];
  return 0.5*sum;
}

double I3C(struct parr *x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  sum += x->f1[i][j]*x->Q[j][j]-x->f1[i][0]*x->Q[j][0];
  for(l=1;l<=j;l++) sum -= (x->f1[i][l]-x->f1[i][l-1])*x->dQv[j][l-1];
  return sum;
}

double I4C(struct parr *x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  for(l=1;l<=j;l++) sum += (x->f2[i][l]+x->f2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCv[j][l-1];
  return 0.5*sum;
}

double I1Q(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x->f1[i][m]*x->Q[m][j]-x->f1[i][j]*x->Q[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x->df1v[i][l-1]*(x->Q[l][j]-x->Q[l-1][j]);
  }
  sum += x->df1v[i][i-1]*(-x->Q[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x->f1[i][l]-x->f1[i][l-1])*x->dQh[l][j];
  }
  return sum;
}

double I2Q(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x->df2v[i][l-1]*(x->Q[i][l]-x->Q[i][l-1])*(2.0-x->Q[l][j]-x->Q[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x->f2[i][l]+x->f2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*(1.0-x->dQh[l][j]);
  return 0.5*sum;
}

double mu_t(struct pmct *z,struct parr *x,int i){
	int l;
	double mu;
	mu = z->T + x->dmu;
	for(l=1;l<=i-z->Nc;l++){
		mu+=(x->Q[i][l]-x->Q[i][l-1])*
			(I3*(x->f1[i][l+1]+x->f2[i][l+1]*x->C[i][l+1])
			+I2*(x->f1[i][l  ]+x->f2[i][l  ]*x->C[i][l  ])
			+I1*(x->f1[i][l-1]+x->f2[i][l-1]*x->C[i][l-1])
		);
	}
	mu -= (1-z->T*z->beta)*fd1(x->C[i][0],z)*x->C[i][0];
	return mu;
}

double E_t(struct pmct *z,struct parr *x,int i){
	int k;
	double E = 0.;
	for(k=0;k<i;k++){
		E -= x->df1v[i][k]*(x->Q[i][k+1]-x->Q[i][k]);
	}
	E -= (z->beta*z->T-1.)*f(x->C[i][0],z)+f(1.,z);
	return E;
}
