#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
#include <sys/stat.h>

#include "segments.cpp"

int main(int argc, char *argv[]){

	int i,k;
	int e = 12;
	int N = 1 << e;
	int Nk;

	printf("Total Size: %d\n",N);


	class Segment * * * G;
	G = (class Segment * * *) malloc(e*sizeof(class Segment * *));


	G[e-1] = (class Segment * *) malloc(N*sizeof(class Segment *));
	for(i=0;i<N;i++) {
			G[e-1][i] = new Segment(i,i,i+0.5);
			//G[e-1][i]->fill(i,i,i+0.5);
	}

	for(k=e-1;k>0;k--) {
		Nk = 1 << k;
		printf("Size: %d\n",Nk);
				
		G[k-1] = (class Segment * *) malloc(Nk*sizeof(class Segment *));
		
		for(i=0;i<Nk;i++) {		
		G[k-1][i] = new Segment();
		G[k-1][i]->anchestor(G[k][2*i],G[k][2*i+1]);
		}
	}
	


/*	class Segment * * G3;
	G3 = (class Segment * *) malloc(N/2*sizeof(class Segment *));
		
	for(i=0;i<N/2;i++) {
		G3[i] = malloc(sizeof(class Segment));
		anchestor(G3[i],G4[2*i],G4[2*i+1]);
	}

	class Segment * * G2;
	G2 = (class Segment * *) malloc(N/4*sizeof(class Segment *));
		
	for(i=0;i<N/4;i++) {
		G2[i] = malloc(sizeof(class Segment));
		anchestor(G2[i],G3[2*i],G3[2*i+1]);
	}

	class Segment * * G1;
	G1 = (class Segment * *) malloc(N/8*sizeof(class Segment *));
		
	for(i=0;i<N/8;i++) {
		G1[i] = malloc(sizeof(class Segment));
		anchestor(G1[i],G2[2*i],G2[2*i+1]);
	}*/


	double P,Pc;

	for(k=e;k>0;k--) {
		P = product(G[k-1][0],G[k-1][1]);
		printf("product(G[%d][0],G[%d][1]) : %f \n",k-1,k-1,P); fflush(stdout);
	}

	for(k=e;k>0;k--) {
		P = product_d(G[k-1][0],G[k-1][1]);
		printf("product_d(G[%d][0],G[%d][1]) : %f \n",k-1,k-1,P); fflush(stdout);
	}

	for(k=e;k>1;k--) {
		P = product3_d(G[k-1][0],G[k-1][1],G[k-1][2]);
		printf("product3_d(G[%d][0],G[%d][1],G[%d][2]) : %f \n",k-1,k-1,k-1,P); fflush(stdout);
	}

	getchar();

	//extinguish_dinasty(G[0][0]);

	P = product(G[0][0],G[0][1]);
	Pc = G[0][0]->get_x_()*G[0][1]->get_x_();
	printf("product(G[%d][0],G[%d][1]) : %f \n",0,0,P); fflush(stdout);
	printf("product(G[%d][0]->x_*G[%d][1]->x_) : %f \n",0,0,Pc); fflush(stdout);

	getchar();
}
