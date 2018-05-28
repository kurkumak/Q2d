#define D1 3./2.
#define D2 -2.
#define D3 1./2.
#define I1 5./12.
#define I2 8./12.
#define I3 -1./12.

class Segment{

private:  
	double * x1; 	//first point of segment
	double * x_;	//average over segment   
	double * x2; 	//last point of segment

	int children;
	class Segment *child[2];

public:

    Segment();
    Segment(double a1, double a2);
    Segment(double a1, double a_, double a2);

	void put_x1(double a1);
	void put_x_(double a);
	void put_x2(double a2); 
	double get_x1();
	double get_x_();
	double get_x2();

	void fill(double a1, double a_, double a2);
	void copy(class Segment *a);
	void extrapolate_next(class Segment a, double a3);
	void extrapolate_prev(double a1, class Segment a);

	void anchestor(class Segment *child0, class Segment *child1);
	void kill_children(class Segment *a);
	void extinguish_dinasty(class Segment *A);

	void match_segment_right(class Segment *h);
	void match_segment_down(class Segment *v);

	friend double product(class Segment *a, class Segment *b);
	friend double product_d(class Segment *a, class Segment *b); //product where the first element is taken as derivative
	friend double product3_d(class Segment *a, class Segment *b, class Segment *c);//product where the first element is taken as derivative
};

class Corner{

private:
	class Segment * h;
	class Segment * v;

public:

	Corner();
	Corner(class Segment *ph, class Segment *pv);
	Corner(double h1, double c, double v2);
	Corner(double h1, double h_, double c, double v_, double v2);

	void fill(double h1, double h_, double c, double v_, double v2);
	void renew();
	void copy(class Corner *a);
	void extrapolate(class Corner *ch, double x0, class Corner *cv);
	void contract(class Corner *ch, class Corner *cv, class Corner *c);
};

/*

void match_segment_right(struct segment *d, struct segment *ch) {
	ch->x1 = d->x2;
}
void match_segment_down(struct segment *d, struct segment *cv) {
	cv->x2 = d->x1;
}


void open_vector(char *name, struct corner *x){
	printf("OPEN VECTOR: %s\n",name);
	FILE *f = fopen(name,"r");
	double d1,d_;
	fscanf(f,"%lf %lf\n",&(x->v->x1),&(x->v->x_)); x++;
	while(fscanf(f,"%lf %lf\n",&d1,&d_)==2) { x->v->x1 = d1; x->v->x_ = d_; (x-1)->v->x2 = d1; (x-1)->v->children = 0; x++; }
	fclose(f);
}

void print_vector(struct corner *x,int L){
	int i;
	for(i=0;i<L;i++) {
		printf("%lf ",x[i].v->x1);
	}
}*/