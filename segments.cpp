#include "segments.h"

//-----------------------SEGMENTS-----------------------
Segment::Segment() 
{
	x1 = new double; x_ = new double; x2 = new double;
	children = 0;
}

Segment::Segment(double a1, double a2) 
{
	x1 = new double; x_ = new double; x2 = new double; 
    fill(a1,0.5*(a1+a2),a2);
}

Segment::Segment(double a1, double a_, double a2) 
{
	x1 = new double; x_ = new double; x2 = new double;
    fill(a1,a_,a2);
}

void Segment::put_x1(double a1) 
{
	*x1 = a1;
}
void Segment::put_x_(double a_) 
{
	*x_ = a_;
}
void Segment::put_x2(double a2) 
{
	*x2 = a2;
}
double Segment::get_x1() 
{
	return *x1;
}
double Segment::get_x_() 
{
	return *x_;
}
double Segment::get_x2()
{
	return *x2;
}

void Segment::fill(double a1, double a_, double a2)
{
	put_x1(a1); put_x_(a_); put_x2(a2); 
	children = 0;
}

void Segment::copy(class Segment *a) 
{
	put_x1( a->get_x1() ); put_x_( a->get_x_() ); put_x2( a->get_x2() );
	
	if(a->children) 
	{
		children = 1;
		child[0] = new class Segment; 		child[1] = new class Segment;
		child[0]->copy(a->child[0]);		child[1]->copy(a->child[1]);
	} 
	else 
	{
		children = 0;
	}
}

void Segment::extrapolate_next(class Segment a, double a3) 
{
	put_x1( a.get_x2() );
	put_x_( I3*a.get_x1() + I2*a.get_x2() + I1*a3 );
	put_x2( a3 );
}

void Segment::extrapolate_prev(double a1, class Segment a) 
{
	put_x1( a1 );
	put_x_( I1*a1 + I2*a.get_x1() + I3*a.get_x2());
	put_x2( a.get_x1() );
}

void Segment::anchestor(class Segment *child0, class Segment *child1) 
{
	put_x1( child0->get_x1() );
	put_x_( (child0->get_x_() + child1->get_x_())/2. );	//average over segment
	put_x2( child1->get_x2() );	
	children = 1;
	child[0] = child0;
	child[1] = child1;
}

void Segment::kill_children(class Segment *a) {
	
	if(a->children) {
		kill_children(a->child[0]);			
		kill_children(a->child[1]);
	}

	free(a);
}

void Segment::extinguish_dinasty(class Segment *A) {
	
	if(A->children) {
		kill_children(A->child[0]);			
		kill_children(A->child[1]);
	}
	A->children=0;
}


//-----------------------PRODUCTS-----------------------

double product(class Segment *a, class Segment *b) 
{	
	double P=0;
	
	if(a->children && b->children) {
		P += product(a->child[0],b->child[0]);
		P += product(a->child[1],b->child[1]);
	}
	else {
		P += a->get_x_()*b->get_x_();
	}
	
	return P;
}

double product_d(class Segment *a, class Segment *b) 
{	
	double P=0;
	
	if(a->children && b->children) {
		P += product_d(a->child[0],b->child[0]);
		P += product_d(a->child[1],b->child[1]);
	}
	else {
		P += (a->get_x2() - a->get_x1()) * b->get_x_();
	}
	
	return P;
}

double product3_d(class Segment *a, class Segment *b, class Segment *c) 
{	
	double P=0;
	
	if(a->children && b->children && c->children) {
		P += product3_d(a->child[0],b->child[0],c->child[0]);
		P += product3_d(a->child[1],b->child[1],c->child[1]);
	}
	else {
		P += (a->get_x2() - a->get_x1()) * b->get_x_() * 0.5*(c->get_x2() + c->get_x1());
	}
	
	return P;
}

//-----------------------CORNERS-----------------------

Corner::Corner() 
{
	h = new Segment; v = new Segment;
}

Corner::Corner(double h1, double c, double v2)
{
	h = new Segment; v = new Segment;
    h->fill(h1,0.5*(h1+c),c);
    v->fill(c,0.5*(c+v2),v2); 
}

Corner::Corner(class Segment *ph, class Segment *pv)
{
	h = ph; v = pv; 
}

Corner::Corner(double h1, double h_, double c, double v_, double v2)
{
	h = new Segment; v = new Segment;
    h->fill(h1,h_,c);
    v->fill(c,v_,v2); 
}

void Corner::fill(double h1, double h_, double c, double v_, double v2)
{
    h->fill(h1,h_,c);
    v->fill(c,v_,v2); 
}

void Corner::renew() 
{
	h = new Segment; v = new Segment;
}

void Corner::copy(class Corner *a) 
{
	h->copy(a->h);
	v->copy(a->v);			
}

void Corner::extrapolate(class Corner *ch, double x0, class Corner *cv) 
{
	h->extrapolate_next(*(ch->h),x0);
	v->extrapolate_prev(x0,*(cv->v));
}

void Corner::contract(class Corner *ch, class Corner *cv, class Corner *c) 
{
	h->anchestor(ch->h,c->h);
	v->anchestor(c->v,cv->v);	
}