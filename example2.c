#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "Segments.h"

int main(int argc, char *argv[]){

	int i,j,k;
	int e = 8;
	int L = 1 << e;
	int N = L*L;

	printf("Total Size: %d\n",N);


	struct corner * Q;
	Q = (struct corner *) malloc((L+1)*sizeof(struct corner));
	for(i=0;i<=L;i++) {
		Q[i].h = malloc(sizeof(struct segment));
		Q[i].v = malloc(sizeof(struct segment));
	}
	open_vector("Q.txt",Q);

	struct corner * dfv1;
	dfv1 = (struct corner *) malloc((L+1)*sizeof(struct corner));
	for(i=0;i<=L;i++) {
		dfv1[i].h = malloc(sizeof(struct segment));
		dfv1[i].v = malloc(sizeof(struct segment));
	}
	open_vector("df1v.txt",dfv1);

	struct corner * Q_an;
	Q_an = (struct corner *) malloc((L+1)*sizeof(struct corner));
	for(i=0;i<L/2;i++) {
		Q_an[i].v = malloc(sizeof(struct segment));
		anchestor(Q_an[i].v,Q[2*i].v,Q[2*i+1].v);
	}

	struct corner * dfv1_an;
	dfv1_an = (struct corner *) malloc((L+1)*sizeof(struct corner));
	for(i=0;i<L/2;i++) {
		dfv1_an[i].v = malloc(sizeof(struct segment));
		anchestor(dfv1_an[i].v,dfv1[2*i].v,dfv1[2*i+1].v);
	}

	double P,a,b=0;
	for(i=0;i<L/2;i++) {
		P+=product_d(Q_an[i].v,dfv1_an[i].v);
	}
	printf("Prod Q*dfv1 = %.8f\n",P);

	double dP = fabs(P/L*2);

	for(i=0;i<L/2;i++) {
		a=product_d(Q_an[i].v,dfv1_an[i].v);	
		Q_an[i].v->children=0;
		b=product_d(Q_an[i].v,dfv1_an[i].v);
		if(fabs(a-b)>dP*0.1) { printf("%f %.8f %.8f \n",dP,a,b); Q_an[i].v->children=1; }
	}

	P=0;
	for(i=0;i<L/2;i++) {
		P+=product_d(Q_an[i].v,dfv1_an[i].v);
	}
	printf("Prod Q*dfv1 = %.8f\n",P);

	/*struct corner * Q_anan;
	Q_anan = (struct corner *) malloc((L+1)*sizeof(struct corner));
	for(i=0;i<L/2;i++) {
		Q_anan[i].v = malloc(sizeof(struct segment));
		anchestor(Q_anan[i].v,Q_an[2*i].v,Q_an[2*i+1].v);
	}*/

	/*struct corner * dfv1_an;
	dfv1_anan = (struct corner *) malloc((L+1)*sizeof(struct corner));
	for(i=0;i<L/2;i++) {
		dfv1_anan[i].v = malloc(sizeof(struct segment));
		anchestor(dfv1_anan[i].v,dfv1[2*i].v,dfv1[2*i+1].v);
	}

	double P=0;
	for(i=0;i<L;i++) {
		P+=product_d(Q[i].v,dfv1[i].v,1);
	}
	printf("Prod Q*dfv1 = %.8f\n",P);

	P=0;
	for(i=0;i<L/2;i++) {
		P+=product_d(Q_an[i].v,dfv1_an[i].v,0);
	}
	printf("Prod Q*dfv1 = %.8f\n",P);*/


	//print_vector(dfv1,257);



	/*for(i=0;i<L;i++) {
	G[i] = (struct corner *) malloc(L*sizeof(struct corner));
		for(j=0;j<N;j++) {
			G[i][j]->h = malloc(sizeof(struct segment));
			G[i][j]->v = malloc(sizeof(struct segment));

			G[i][j]->h->x1=tanh((i-j)/100.);
			G[i][j]->h->x2=tanh((i-j+1)/100.);

			G[i][j]->v->x1=tanh(-(i-1-j)/10.);
			G[i][j]->v->x2=tanh(-(i-j)/10.);

			G[i][j]->x_=i;
			G[i][j]->h->children=0;
			G[i][j]->v->children=0;
		}
	}*/

	/*for(k=e-1;k>0;k--) {
		Nk = 1 << k;
		printf("Size: %d\n",Nk);
				
		G[k-1] = (struct segment * *) malloc(Nk*sizeof(struct segment *));
		
		for(i=0;i<Nk;i++) {		
		G[k-1][i] = malloc(sizeof(struct segment));
		anchestor(G[k-1][i],G[k][2*i],G[k][2*i+1]);
		}
	}*/
}







