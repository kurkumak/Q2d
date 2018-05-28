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

void CQ_Integration(struct pmct *z,struct parr *x,struct psys *w);
int single_iteration(struct pmct *z,struct parr *x, struct psys *w);
int time_step(int i,struct pmct *z,struct parr *x);

double self_consistence_loop(double *gC,double *gQ,double D,double mu,struct pmct *z,struct parr *x,int i,int j);
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

	CQ_Integration(&z,&x,&w);

	printf("\n-------------------------------------------------------END-------------------------------------------------------\n");

	return(0);
}

void CQ_Integration(struct pmct *z,struct parr *x, struct psys *w){

	array_initialization(z,x);

	int rpt=0,itr=0;

/*------------------------------------------------------------*/
	initialarray(z,x); 	// prepare the array btwn 0 <= i,j <= Nt/2
/*------------------------------------------------------------*/

	write_C_0(z,x,w,0,z->Nt2);

	while(itr <= z->itr && rpt<z->rpt){
		rpt = single_iteration(z,x,w);
		printf(" %dth/%d cycle - window_time: %.2e \n",itr,z->itr,x->dt*z->Nt);

		itr++;
	}

	final_write_C(z,x,w);
}

int single_iteration(struct pmct *z,struct parr *x, struct psys *w) {

	int i=z->Nt2+1;
	int rpt=0;

	while(i<=z->Nt && rpt<z->rpt){
/*------------------------------------------------------------*/
		rpt = time_step(i,z,x);	// propagate the solution the array btwn Nt/2 <= i,j <= Nt
		i++;
/*------------------------------------------------------------*/
	}


//------------------PRINTING--------------------------
	write_C_0(z,x,w,z->Nt2+1,z->Nt);
	write_C(z,x,w,1);

	snap_config(x->C,x->dt*z->Nt,z,w);
//----------------------------------------------------


/*------------------------------------------------------------*/
	//printf("\nmu_before = %f",mu_t(z,x,z->Nt));
	contract(z,x,&(x->dt),&(x->dmu));	// double the size of the system
	//printf("\nmu_after = %f\n",mu_t(z,x,z->Nt2));
/*------------------------------------------------------------*/

	return rpt;
}

int time_step(int i,struct pmct *z,struct parr *x){

	int j,rpt;
	double D,err2,err2_tmp;

	double gC;
	double gQ;

	extrapolate(z,x,i);

	// (3) Go to the SC (self-consistence) loop
	rpt = 0;
	err2 =1.0;

	while(err2 >= z->eps*z->eps && rpt < z->rpt){
		
		err2 = 0.0;

		for(j=i-z->Nc-1;j>=0;j--){

			err2_tmp = self_consistence_loop(&gC,&gQ,D,x->mu[i],z,x,i,j);
			if (err2_tmp>err2) { err2=err2_tmp; }
			
			// renew all variable
			x->C[i][j] += gC;
			x->Q[i][j] += gQ;
			x->f1[i][j] = fd1(x->C[i][j],z);
			x->f2[i][j] = fd2(x->C[i][j],z);
		
			extrapolate_step(z,x,i,j);
		}

		rpt++;	
	}

	x->E[i] = E_t(z,x,i);

	printf("\r\t%4d\t%2.1e:%2.1e\t  %d\t",i,err2,z->eps*z->eps,rpt); fflush(stdout);
	return rpt;
}

double self_consistence_loop(double *gC, double *gQ, double D, double mu, struct pmct *z, struct parr *x,int i, int j){

	//evaluate Lagrange multiplier
	x->mu[i] = mu_t(z,x,i);

	//factor that multiplies C[i][j] and Q[i][j] in the self-loop consistence
	D = D0/x->dt  + x->mu[i] + x->df1v[i][i-1];


	int m;
	double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q,bf1C;
	double C2,C1,Q2,Q1;
	double T2;

	m = (int)(0.5*(i+j));

	i1C = I1C(x,i,j,m);
	i2C = I2C(x,i,j,m);
	i3C = I3C(x,i,j);
	i4C = I4C(x,i,j);
	i1Q = I1Q(x,i,j,m);
	i2Q = I2Q(x,i,j,m);
	i3Q = i3C;
	i4Q = i4C;
	bf1C = (1-z->T*z->beta)*fd1(x->C[i][0],z)*x->C[j][0];

	T2 = z->T*z->T;

	C2 = x->C[i-2][j]/x->dt;
	C1 = x->C[i-1][j]/x->dt;
	Q2 = x->Q[i-2][j]/x->dt;
	Q1 = x->Q[i-1][j]/x->dt;



	*gC  = -D2*C2 -D1*C1; 
	*gC += -i1C+i2C+i3C+i4C;
	*gC -= bf1C;
	*gC /= D;
	*gC -= x->C[i][j];

	*gQ  = -D2*Q2 -D1*Q1;
	*gQ += -T2 + mu; 
	*gQ += -i1Q-i2Q-i3Q-i4Q;
	*gQ += bf1C;
	*gQ /= D;
	*gQ -= x->Q[i][j];
	
	return power(*gC,2)+power(*gQ,2);
}

double I1C(struct parr *x,int i,int j,int m){
	int l;
	double sum;
	sum = 0.0;

	for(l=m+1;l<=i-1;l++){
		sum += x->df1v[i][l-1]*(x->C[l][j]-x->C[l-1][j]);
	}

	sum += x->f1[i][m]*x->C[m][j]-x->f1[i][j]*x->C[j][j];

	for(l=j+1;l<=m;l++){
		sum -= (x->f1[i][l]-x->f1[i][l-1])*x->dCh[l][j];
	}

    sum += x->df1v[i][i-1]*(-x->C[i-1][j]); //Term of the step expansion

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

	for(l=m+1;l<=i-1;l++){
		sum += x->df1v[i][l-1]*(x->Q[l][j]-x->Q[l-1][j]);
	}

	sum += x->f1[i][m]*x->Q[m][j]-x->f1[i][j]*x->Q[j][j];

	for(l=j+1;l<=m;l++){
		sum -= (x->f1[i][l]-x->f1[i][l-1])*x->dQh[l][j];
	}

    sum += x->df1v[i][i-1]*(-x->Q[i-1][j]); //Term of the step expansion

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
	mu = z->T*z->T + x->dmu;
	for(l=1;l<=i-z->Nc;l++){
		mu+=(x->Q[i][l]-x->Q[i][l-1])*
			(I2*(x->f1[i][l+1]+x->f2[i][l+1]*x->C[i][l+1])
			+I1*(x->f1[i][l  ]+x->f2[i][l  ]*x->C[i][l  ])
			+I0*(x->f1[i][l-1]+x->f2[i][l-1]*x->C[i][l-1])
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
