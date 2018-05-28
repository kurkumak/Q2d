#define thisfile "AdaptQ2d.c"
#include "libAdaptQ2d.h"

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
void step(int i,struct pmct *z,struct parr *x);

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
void SpecialE_t(struct pmct *z,struct parr *x,int i);


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

	int i,itr=0;

/*------------------------------------------------------------*/
	initialarray(z,x); 	// prepare the array btwn 0 <= i,j <= Nt/2
/*------------------------------------------------------------*/

	write_C_0(x,w,1,z->Nt2);
	//write_C(z,x,w,1);

	while(itr <= z->itr){

				//printf("\nmu_before = %f\n\n",mu_t(z,x,z->Nt2)); 
		for(i=z->Nt2+1;i<=z->Nt;i++){
/*------------------------------------------------------------*/
			step(i,z,x);	// propagate the solution the array btwn Nt/2 <= i,j <= Nt
			//getchar();
/*------------------------------------------------------------*/
		}

		//snap_config(x->C,x->dt*z->Nt,z,w);
		//if(itr==20) { snap_vector(x->f1,x->df1v,"df1v",x->dt*z->Nt,z,w); snap_vector(x->Q,x->dQv,"Q",x->dt*z->Nt,z,w); SpecialE_t(z,x,z->Nt); }

		write_C_0(x,w,z->Nt2+1,z->Nt);
		write_C(z,x,w,1);

/*------------------------------------------------------------*/
		printf("\nmu_before = %.12f",mu_t(z,x,z->Nt));  //getchar();
		contract(z,x,&(x->dt),&(x->dmu));	// double the size of the system
		printf("\nmu_after = %.12f\n",mu_t(z,x,z->Nt2));
/*------------------------------------------------------------*/
		printf(" %dth/%d cycle - window_time: %.2e \n",itr,z->itr,x->dt*z->Nt);
		itr++;
	}
}

void step(int i,struct pmct *z,struct parr *x){

	int j,scmax;
	double D,err2,err2_temp;
	double New_C,New_Q;

	double * gC = (double *) calloc ((z->Nt),sizeof(double));
	double * gQ = (double *) calloc ((z->Nt),sizeof(double));

		
	//printf("BEFORE EXTRAPOLATE"); getchar();
	extrapolate(z,x,i);

	// (3) Go to the SC (self-consistence) loop
	scmax = 0;
	err2 =1.0;

	while( err2 >= z->eps*z->eps && scmax < z->rpt){
		err2 = 0.0;

		//****** PART ----> j<i-1
		for(j=i-z->Nc-1;j>=0;j--){

			x->mu[i].x1 = mu_t(z,x,i);
			D = D1/x->dt  + x->mu[i].x1 + x->f1[i][i-1].v->x_;
			    //printf("MU(%d,%d)%.12e %.12e %.12e",i,j,D1/x->dt,x->mu[i].x1,x->f1[i][i-1].v->x_); getchar();


			err2_temp = SC2(gC,gQ,D,x->mu[i].x1,z,x,i,j); ///
			if (err2_temp>err2) { err2=err2_temp; }
		// renew all variable

			//if(j==0) printf("\n\nprima: %.12f, difference: %.12e\n\n",x->C[i][j].v->x1,gC[j]);
			New_C = x->C[i][j].v->x1+gC[j];
			New_Q = x->Q[i][j].v->x1+gQ[j];
			extrapolate_corner(&(x->C[i-1][j]), &(x->C[i][j+1]), New_C, &(x->C[i][j]));
			extrapolate_corner(&(x->Q[i-1][j]), &(x->Q[i][j+1]), New_Q, &(x->Q[i][j]));
			extrapolate_corner(&(x->f1[i-1][j]), &(x->f1[i][j+1]), fd1(New_C,z), &(x->f1[i][j]));
			extrapolate_corner(&(x->f2[i-1][j]), &(x->f2[i][j+1]), fd2(New_C,z), &(x->f2[i][j]));

			if(i<z->Nt) {
				match_segment_right(x->C[i][j].h,x->C[i+1][j].h);
				match_segment_right(x->Q[i][j].h,x->Q[i+1][j].h);
				match_segment_right(x->f1[i][j].h,x->f1[i+1][j].h);
				match_segment_right(x->f2[i][j].h,x->f2[i+1][j].h);
			}
			if(j>0) {
				match_segment_down(x->C[i][j].v,x->C[i][j-1].v);
				match_segment_down(x->Q[i][j].v,x->Q[i][j-1].v);
				match_segment_down(x->f1[i][j].v,x->f1[i][j-1].v);
				match_segment_down(x->f2[i][j].v,x->f2[i][j-1].v);
			}			



							//print_C_screen(i,x);


			/*x->C[i][j] += gC[j];
			x->Q[i][j] += gQ[j];
			x->f1[i][j] = fd1(x->C[i][j],z);
			x->f2[i][j] = fd2(x->C[i][j],z);
		
			extrapolate_step(z,x,i,j);*/
		}

		scmax++;	
	}

	free(gC);
	free(gQ);

	x->E[i].x1 = E_t(z,x,i);

	z->f_err = err2;
	z->scmaxx[i]=scmax;
	printf("\r\t%4d\t%2.1e:%2.1e\t  %d\t",i,z->f_err,z->eps,z->scmaxx[i]); fflush(stdout);
}

/*double SC1(double *gC, double *gQ, double D, double mu, struct pmct *z, struct parr *x,int i,int j){
	int m = (int)(0.5*(i+j));
	double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q;
	i1C = I1Cs(x,i,j,m);
	i2C = I2C(x,i,j,m); //x->df2v[i][j]*(x->Q[i][j+1]-x->Q[i][j])*x->dCh[i][j];
//i2C = 0.25*(x->df2v[i][i-1]*(x->Q[i][i]-x->Q[i][i-1])*(x->C[i][i-1]+x->C[i-1][i-1])+(x->f2[i][i]+x->f2[i][i-1])*(x->Q[i][i]-x->Q[i][i-1])*x->dCh[i][i-1]);
	i3C = I3C(x,i,j);
	i4C = I4C(x,i,j);
	i1Q = I1Qs(x,i,j,m);
	i2Q = I2Q(x,i,j,m);//x->df2v[i][j]*(x->Q[i][j+1]-x->Q[i][j])*(1.-x->dQh[i][j]);
//i2Q = 0.25*(x->df2v[i][i-1]*(x->Q[i][i]-x->Q[i][i-1])*(2.0-x->Q[i][i-1]-x->Q[i-1][i-1])+(x->f2[i][i]+x->f2[i][i-1])*(x->Q[i][i]-x->Q[i][i-1])*(1.0-x->dQh[i][i-1]));
	i3Q = i3C;
	i4Q = i4C;

	gC[j] = +1.0/x->dt*x->C[i-1][j]-i1C+i2C+i3C+i4C;
	gC[j]-= (1-z->T*z->beta)*fd1(x->C[i][0],z)*x->C[j][0];
	gC[j]/= D;
	gC[j]-= x->C[i][j];
	gQ[j] = -z->T+mu+1.0/x->dt*x->Q[i-1][j]-i1Q-i2Q-i3Q-i4Q;
	gQ[j]+= (1-z->T*z->beta)*fd1(x->C[i][0],z)*x->C[j][0];
	gQ[j]/= D;
	gQ[j]-= x->Q[i][j];
	
	return gC[j]*gC[j]+gQ[j]*gQ[j];
}*/

double SC2(double *gC, double *gQ, double D, double mu, struct pmct *z, struct parr *x,int i, int j){
  int m;
  double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q,V,C3,C2,C1,Q3,Q2,Q1,dt;
  m = (int)(0.5*(i+j));
  i1C = I1C(x,i,j,m);    //printf("1c1"); fflush(stdout);
  i2C = I2C(x,i,j,m);    //printf("1c2"); fflush(stdout);
  i3C = I3C(x,i,j);      //printf("1c3"); fflush(stdout);
  i4C = I4C(x,i,j);      //printf("1c4"); fflush(stdout);
  i1Q = I1Q(x,i,j,m);    //printf("1q1"); fflush(stdout);
  i2Q = I2Q(x,i,j,m);    //printf("1q2"); fflush(stdout);
  i3Q = i3C;
  i4Q = i4C;

  //printf("QUI%.8e %.8e %.8e %.8e %.8e %.8e\n\n",i1C,i2C,i3C,i4C,i1Q,i2Q); getchar();


  V = (1-z->T*z->beta)*fd1(x->C[i][0].v->x1,z)*x->C[j+1][0].h->x1; //DA controllare x->C[j+1][0].v->x1 (j+1)?

  C1 = x->C[i][j].v->x1; C2 = x->C[i-1][j].v->x1; C3 = x->C[i-2][j].v->x1;
  Q1 = x->Q[i][j].v->x1; Q2 = x->Q[i-1][j].v->x1; Q3 = x->Q[i-2][j].v->x1;
  dt = x->dt;

  ////printf("%f %f %f",C3,C2,-i1C+i2C+i3C+i4C); getchar();

      //printf("%.8e %.8e %.8e",gC[j],-D3*C3/dt - D2*C2/dt,-i1C+i2C+i3C+i4C); getchar();

  gC[j] = -D3*C3/dt - D2*C2/dt  -i1C+i2C+i3C+i4C;

    //printf("%.12e ",gC[j]); getchar();

  gC[j]-= V;
    //printf("%.12e ",gC[j]); getchar();

  gC[j]/= D;
    //printf("%.12e ",gC[j]); getchar();

  gC[j]-= C1;
    //printf("%.12e %.12e ",C1,gC[j]); getchar();
  gQ[j] = - z->T + mu - D3*Q3/dt - D2*Q2/dt  -i1Q-i2Q-i3Q-i4Q;
      //printf("%.8e %.8e %.8e",gQ[j],- z->T + mu - D3*Q3/dt - D2*Q2/dt, -i1Q-i2Q-i3Q-i4Q); getchar();

  gQ[j]+= V;
    //printf("%.8e ",gQ[j]); getchar();

  gQ[j]/= D;
    //printf("%.8e ",gQ[j]); getchar();

  gQ[j]-= Q1;
	
  //printf("%.8e %.8e",gC[j],gQ[j]); getchar();
  return gC[j]*gC[j]+gQ[j]*gQ[j];
}

double I1Cs(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=j;l<=i-1;l++){
    sum += product_d(x->C[l+1][j].h,x->f1[i][l].v);
    	/*x->df1v[i][l]*(x->C[l+1][j]-x->C[l][j]);*/
  }
  return sum;
}

double I1C(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
    ////printf("%d %d %d %f %f %f %f ",m,i,j,x->f1[i][m].v,x->C[m][j].v->x1,x->f1[i][j].v->x1,x->C[j+1][j].h->x1); getchar();
  sum += x->f1[i][m].v->x1*x->C[m][j].v->x1 - x->f1[i][j].v->x1*x->C[j+1][j].h->x1;
   //printf("%.8f ",sum); //getchar();
  //x->f1[i][m]*x->C[m][j]-x->f1[i][j]*x->C[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += product_d(x->C[l][j].h,x->f1[i][l-1].v);
    //x->df1v[i][l-1]*(x->C[l][j]-x->C[l-1][j]);
  }
  //printf("%.8f ",sum); //getchar();

  sum += x->f1[i][i-1].v->x_*(-x->C[i-1][j].v->x1);
  //x->df1v[i][i-1]*(-x->C[i-1][j]);

  //printf("%.8f ",sum); //getchar();

  for(l=j+1;l<=m;l++){
    sum -= product_d(x->f1[i][l-1].v,x->C[l][j].h);
    //(x->f1[i][l]-x->f1[i][l-1])*x->dCh[l][j];
       //printf("ic%d %.12f %.12f %.12f %.12f ",l,x->f1[i][l].v->x1,x->f1[i][l-1].v->x1,x->C[l][j].h->x_,sum); //getchar();

  }
 	//printf("%.8f ",sum); //getchar();

  return sum;
}

double I2C(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += product3_d(x->Q[i][l-1].v,x->f2[i][l-1].v,x->C[l][j].h);
  	/*x->df2v[i][l-1]*(x->Q[i][l]-x->Q[i][l-1])*(x->C[l][j]+x->C[l-1][j]);*/
  for(l=j+1;l<=m;l++) sum += product3_d(x->Q[i][l-1].v,x->C[l][j].h,x->f2[i][l-1].v);
  	/*(x->f2[i][l]+x->f2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCh[l][j];*/
  return sum; //0.5*
}

double I3C(struct parr *x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  sum += x->f1[i][j].v->x1*x->Q[j+1][j].h->x1 - x->f1[i][0].v->x1*x->Q[j+1][0].h->x1;
  for(l=1;l<=j;l++) sum -= product_d(x->f1[i][l-1].v,x->Q[j][l-1].v);
  	//(x->f1[i][l]-x->f1[i][l-1])*x->dQv[j][l-1];
  return sum;
}

double I4C(struct parr *x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  for(l=1;l<=j;l++) sum += product3_d(x->Q[i][l-1].v,x->C[j][l-1].v,x->f2[i][l-1].v);
  	//(x->f2[i][l]+x->f2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCv[j][l-1];
  return sum; //0.5*
}

double I1Qs(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=j;l<=i-1;l++){
    sum += product_d(x->Q[l+1][j].h,x->f1[i][l].v);
    //x->df1v[i][l]*(x->Q[l+1][j]-x->Q[l][j]);
  }
  return sum;
}


double I1Q(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
    //printf("%d %d %d %f %f %f %f ",m,i,j,x->f1[i][m].v,x->Q[m][j].v->x1,x->f1[i][j].v->x1,x->Q[j+1][j].h->x1); getchar();
  sum += x->f1[i][m].v->x1*x->Q[m][j].v->x1 - x->f1[i][j].v->x1*x->Q[j+1][j].h->x1;
  // printf("%f ",sum); getchar();
  //x->f1[i][m]*x->Q[m][j]-x->f1[i][j]*x->Q[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += product_d(x->Q[l][j].h,x->f1[i][l-1].v);
    //x->df1v[i][l-1]*(x->Q[l][j]-x->Q[l-1][j]);
  }
 // printf("%f ",sum); getchar();

  sum += x->f1[i][i-1].v->x_*(-x->Q[i-1][j].v->x1);
  //x->df1v[i][i-1]*(-x->Q[i-1][j]);

  //printf("%f ",sum); getchar();

  for(l=j+1;l<=m;l++){
    sum -= product_d(x->f1[i][l-1].v,x->Q[l][j].h);
    //(x->f1[i][l]-x->f1[i][l-1])*x->dCh[l][j];
  }
 //  printf("%f ",sum); getchar();

  return sum;
}

double I2Q(struct parr *x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += -product3_d(x->Q[i][l-1].v,x->f2[i][l-1].v,x->Q[l][j].h) + product_d(x->Q[i][l-1].v,x->f2[i][l-1].v);
  	//x->df2v[i][l-1]*(x->Q[i][l]-x->Q[i][l-1])*(2.0-x->Q[l][j]-x->Q[l-1][j]);
  for(l=j+1;l<=m;l++) sum += -product3_d(x->Q[i][l-1].v,x->Q[l][j].h,x->f2[i][l-1].v) + 
  	0.5*(x->f2[i][l].v->x1+x->f2[i][l-1].v->x1)*(x->Q[i][l].v->x1-x->Q[i][l-1].v->x1);//product_d(x->Q[i][l-1].v,x->f2[i][l-1].v);
  	//(x->f2[i][l]+x->f2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*(1.0-x->dQh[l][j]);
  return sum; //0.5*
}

double mu_t(struct pmct *z,struct parr *x,int i){
	int l;
	double mu;
	mu = z->T;// + x->dmu;
	for(l=1;l<i-1;l++){
	mu += (x->Q[i][l].v->x1-x->Q[i][l-1].v->x1)*
			(I3*(x->f1[i][l+1].v->x1+x->f2[i][l+1].v->x1*x->C[i][l+1].v->x1)
			+I2*(x->f1[i][l  ].v->x1+x->f2[i][l  ].v->x1*x->C[i][l  ].v->x1)
			+I1*(x->f1[i][l-1].v->x1+x->f2[i][l-1].v->x1*x->C[i][l-1].v->x1));
	
	/*for(l=1;l<=i-z->Nc;l++){
		mu+=(x->Q[i][l].v->x1-x->Q[i][l-1].v->x1)*
			(I3*(x->f1[i][l+1].v->x1+x->f2[i][l+1].v->x1*x->C[i][l+1].v->x1)
			+I2*(x->f1[i][l  ].v->x1+x->f2[i][l  ].v->x1*x->C[i][l  ].v->x1)
			+I1*(x->f1[i][l-1].v->x1+x->f2[i][l-1].v->x1*x->C[i][l-1].v->x1));
		//product_d(x->Q[i][l-1].v,x->f1[i][l-1].v)+product3_d(x->Q[i][l-1].v,x->f2[i][l-1].v,x->C[i][l-1].v);
		/*mu+=(x->Q[i][l]-x->Q[i][l-1])*
			(I3*(x->f1[i][l+1]+x->f2[i][l+1]*x->C[i][l+1])
			+I2*(x->f1[i][l  ]+x->f2[i][l  ]*x->C[i][l  ])
			+I1*(x->f1[i][l-1]+x->f2[i][l-1]*x->C[i][l-1])
		);*/
	}
	mu -= (1-z->T*z->beta)*fd1(x->C[i][0].v->x1,z)*x->C[i][0].v->x1;
	//printf("mu%f**%f-%f-%f::",(1-z->T*z->beta)*fd1(x->C[i][0].v->x1,z)*x->C[i][0].v->x1,(1-z->T*z->beta),fd1(x->C[i][0].v->x1,z),x->C[i][0].v->x1);
	return mu;
}

double E_t(struct pmct *z,struct parr *x,int i){
	int k;
	double E = 0.;
	for(k=0;k<i;k++){
		E -= product_d(x->Q[i][k].v,x->f1[i][k].v);
		/*E -= x->df1v[i][k]*(x->Q[i][k+1]-x->Q[i][k]);*/
	}
	E -= (z->beta*z->T-1.)*f(x->C[i][0].v->x1,z)+f(1.,z);
	return E;
}

/*void SpecialE_t(struct pmct *z,struct parr *x,int i){
	int k;
	double E = 0.;
	for(k=0;k<i;k++){
		E -= x->df1v[i][k]*(x->Q[i][k+1]-x->Q[i][k]);
	}
	printf("Prod = %.5f\n",E);
}*/
  
 //fprintf(fout,"# x->R      =z->T*z->beta=%.2e \n",   x->R);

/*#ifndef PRINT
	print_vector(z->sh,x->C,z->Nt);		//PRINT
	print_vector(z->sh2,x->Q,z->Nt);		//PRINT
#endif*/
