#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>
using namespace std;

double fixed_point(double (*f)(double), double x){
//	cout << x << endl;
	double y=f(x);
	if(x==y)
		return x;
	else
		return fixed_point(f,y);
}
double recurrence(double (*f)(double), double x, int iterations){
	if(iterations>0){
		x=f(x);
		cout << x << endl;
		iterations--;
		return recurrence(f,x,iterations);
	}
	else
		return x;
}

double aitken(double (*f)(double),double x, int iterations){
		double a=recurrence(f,x,iterations);
		double b=recurrence(f,x,iterations+1);
		double c=recurrence(f,x,iterations+2);
		double d= a-(pow(b-a,2)/(c-(2*b)+a));
		cout << iterations<< ": " << a << " " << d << endl;
		return d;
}

double bisection(double a, double b, double (*f)(double), double tol, int &iterations){
	if(f(a)<=tol && f(a)>=-tol)
		return a;
	else if(f(b)<=tol && f(b)>=-tol)
		return b;
	else{
		iterations++;	
		long double x=(a+b)/2;
		if((f(a)>tol && f(x)<-tol) || (f(a)<-tol && f(x)>tol))
			return bisection(a,x,f,tol, iterations);
		else
			return bisection(x,b,f,tol, iterations);
	}
}

double bisection(double a, double b, double (*f)(double), double tol){
	if(f(a)<=tol && f(a)>=-tol)
		return a;
	else if(f(b)<=tol && f(b)>=-tol)
		return b;
	else{
		long double x=(a+b)/2;
		if((f(a)>tol && f(x)<-tol) || (f(a)<-tol && f(x)>tol))
			return bisection(a,x,f,tol);
		else
			return bisection(x,b,f,tol);
	}
}

double secant(double (*f)(double), double xzero, double xone, double tol){
	double a=f(xzero);
	double b=f(xone);
	double x_n=xone-(b*((xone-xzero)/(b-a)));
	double c=f(x_n);
	cout << x_n << " " << c << endl;
	if(c>=-tol && c<=tol)
		return x_n;
	else{
		xzero=xone;
		xone=x_n;
		return secant(f,xzero,xone,tol);
	}
}

double newton(double (*f)(double), double(*fprime)(double), double xzero, double tol){
	double x_n=xzero-(f(xzero)/fprime(xzero));
	double a=f(x_n);
//	cout << a << endl;
//	cout << xzero << " " << x_n << endl;
	if(a>=-tol && a<=tol)
		return x_n;
	else
		return newton(f,fprime,x_n,tol);
}

double falsi(double (*f)(double),double a, double b, double tol, int &iterations){
	double x=a-((b-a)/(f(b)-f(a)))*f(a);
	iterations++;
	if(f(x)>=-tol && f(x)<=tol)
		return x;
	else if(f(x)<-tol)
		return falsi(f,x,b,tol,iterations);
	else
		return falsi(f,a,x,tol,iterations);
}

double steffensen(double (*g)(double),double p0,double it){
	double p1,p2,p;
	for(int c=0;c<it;c++){
		p1=g(p0);
		p2=g(p1);
		p=p0-(pow(p1-p0,2)/(p2-(2*p1)+p0));
		cout << p << endl;
		p0=p;
	}
	return p;
}

double muller(double (*f)(double), double p0, double p1, double p2, double tol, double it){
	double h1,h2,s1,s2,d,D,E,b,h,p;
	h1=p1-p0;
	h2=p2-p1;
	s1=(f(p1)-f(p0))/h1;
	s2=(f(p2)-f(p1))/h2;
	d=(s2-s1)/(h2+h1);
	for(int c=3;c<=it;c++){
		b=s2+(h2*d);
		D=pow(pow(b,2)-(4*f(p2)*d),0.5);
		if(abs(b-D)<abs(b+D))
			E=b+D;
		else
			E=b-D;
		h=-2*f(p2)/E;
		p=p2+h;
		if(abs(h)<tol){
			cout << "found after " << c << " iterations.\n";
			return(p);
		}
		p0=p1;
		p1=p2;
		p2=p;
		h1=p1-p0;
		h2=p2-p1;
		s1=(f(p1)-f(p0))/h1;
		s2=(f(p2)-f(p1))/h2;
		d=(s2-s1)/(h2+h1);
	}
	cout << "method failed" << endl;
	return 1;
}

vector<double> lagrange(vector<double> x, vector<double> y,vector<double> L_x){
	vector<vector<double> > L;
	vector<double> L_y;
	for(int c=0;c<L_x.size();c++){
		L.push_back(L_y);
		for(int d=0;d<x.size();d++){
			L[c].push_back(y[d]);
			for(int e=0;e<x.size();e++){
				if(e!=d)
					L[c][d]*=(L_x[c]-x[e])/(x[d]-x[e]);
			}
		}
	}
	for(int c=0;c<L_x.size();c++){
		L_y.push_back(0);
		for(int d=0;d<x.size();d++)
			L_y[c]+=L[c][d];
	}
	return L_y;
}

double neville(vector<double> x_val, vector<double> y_val, double x){
	double Q[y_val.size()][y_val.size()];
	for(int c=0;c<y_val.size();c++)
		Q[c][0]=y_val[c];
	for(int c=1;c<y_val.size();c++){
		for(int d=1;d<=c;d++)
			Q[c][d]=(((x-x_val[c-d])*Q[c][d-1])-((x-x_val[c])*Q[c-1][d-1]))/(x_val[c]-x_val[c-d]);
	}
	return Q[y_val.size()-1][y_val.size()-1];
}

double** divided_difference(vector<double> x_val, vector<double> y_val){
	double **table=new double*[y_val.size()];
	for(int c=0;c<y_val.size();c++){
		table[c]=new double[y_val.size()];
		table[c][0]=y_val[c];
	}
	for(int c=1;c<y_val.size();c++){
		for(int d=1;d<=c;d++){
			table[c][d]=(table[c][d-1]-table[c-1][d-1])/(x_val[c]-x_val[c-d]);
//			cout << table[c][d] << "\t";
		}
		cout << endl;
	}
	return table;
}

double forward_difference(vector<double> x_val, vector<double> y_val, double x){
	double **table=divided_difference(x_val,y_val);
	vector<double> coef;
	double h=x_val[1]-x_val[0];
	double s=(x-x_val[0])/h;
	coef.push_back(table[0][0]);
	for(int c=1;c<y_val.size();c++){
		coef.push_back(table[c][c]);
		for(double d=s;d>=s-c+1;d--){
			coef[c]*=d;
		}
		coef[c]*=pow(h,c);	
	}
	double y=0;
	for(int c=0;c<coef.size();c++)
		y+=coef[c];
	return y;
}

double backward_difference(vector<double> x_val, vector<double> y_val, double x){
	double **table=divided_difference(x_val,y_val);
	vector<double> coef;
	double h=x_val[x_val.size()-1]-x_val[x_val.size()-2];
	double s=(x-x_val[x_val.size()-1])/h;
	coef.push_back(table[y_val.size()-1][0]);
	for(int c=1;c<y_val.size();c++){
		coef.push_back(table[y_val.size()-1][c]);
		coef[c]*=pow(-1,c);
		for(double d=-s;d>=-s-c+1;d--){
			coef[c]*=d;
		}
		coef[c]*=pow(h,c);
	}	
	double y=0;
	for(int c=0;c<coef.size();c++)
		y+=coef[c];
	return y;
}

double centered_difference(vector<double> x_val, vector<double> y_val, double x){
	double **table=divided_difference(x_val,y_val);
	vector<double> coef;
	int index;
	if(x_val.size()%2==1)
		index=x_val.size()/2;
	else
		index=x_val.size()/2-1;
	double h=x_val[index]-x_val[index-1];
	double s=(x-x_val[index])/h;
	for(int c=0;c<y_val.size();c++){
		if(c%2==0){
			coef.push_back(table[index+(c/2)][c]);
			if(c!=0)
				coef[c]*=pow(s,2);
			if(c>2){
				for(double d=1;d<=(c/2)-1;d++)
					coef[c]*=pow(s,2)-pow(d,2);
			}
		}
		else{
			coef.push_back(table[index+((c-1)/2)][c]+table[index+((c-1)/2)+1][c]);
			coef[c]/=2;
			coef[c]*=s;

			if(c>2){
				for(double d=1;d<=(c-1)/2;d++)
					coef[c]*=pow(s,2)-pow(d,2);
			}
		}
		coef[c]*=pow(h,c);
		//cout << coef[c] << endl;
	}
	double y=0;
	for(int c=0;c<coef.size();c++)
		y+=coef[c];
	return y;
}

vector<double> hermite(vector<double> x_val, vector<double> y_val, vector<double> y_prime, vector<double> x){
	double table[y_val.size()*2][y_val.size()*2+1];
	vector<double> y;
	double term;
	for(int c=0;c<y_val.size();c++){
		table[2*c][0]=x_val[c];
		table[2*c+1][0]=x_val[c];
		table[2*c][1]=y_val[c];
		table[2*c+1][1]=y_val[c];
	}
	for(int c=1;c<y_val.size()*2;c++){
		if(c%2==1)
			table[c][2]=y_prime[(c-1)/2];
		else
			table[c][2]=(table[c][1]-table[c-1][1])/(x_val[c/2]-x_val[c/2-1]);
	}
	for(int c=3;c<y_val.size()*2+1;c++){
		for(int d=c-1;d<y_val.size()*2;d++)
			table[d][c]=(table[d][c-1]-table[d-1][c-1])/(table[d][0]-table[d-(c-1)][0]);
	}
	for(int c=0;c<x.size();c++){
		y.push_back(table[0][1]);
		for(int d=2;d<y_val.size()*2+1;d++){
			term=table[d-1][d];
			for(int e=0;e<d-1;e++)
				term*=x[c]-table[e][0];
			y[y.size()-1]+=term;
		}
	}
	return y;
}

vector<vector<double> > cubic_spline(vector<double> x_val, vector<double> y_val){
	vector<vector<double> > coefs;
	vector<double> coef(4);
	double h[x_val.size()], alpha[x_val.size()], l[x_val.size()], mu[x_val.size()], z[x_val.size()];
	double b[x_val.size()], e[x_val.size()], d[x_val.size()]; 
	for(int c=0;c<x_val.size()-1;c++)
		h[c]=x_val[c+1]-x_val[c];
	for(int c=1;c<x_val.size()-1;c++)
		alpha[c]=(3/h[c]*(y_val[c+1]-y_val[c]))-(3/h[c]*(y_val[c]-y_val[c-1]));
	l[0]=1;
	mu[0]=0;
	z[0]=0;
	for(int c=1;c<x_val.size()-1;c++){
		l[c]=2*(x_val[c+1]-x_val[c-1])-(h[c-1]*mu[c-1]);
		mu[c]=h[c]/l[c];
		z[c]=(alpha[c-1]-(h[c-1]*z[c-1]))/l[c];
	}
	l[x_val.size()-1]=1;
	z[x_val.size()-1]=0;
	e[x_val.size()-1]=0;
	for(int c=x_val.size()-2;c>=0;c--){
		e[c]=z[c]-(mu[c]*e[c+1]);
		b[c]=((y_val[c+1]-y_val[c])/h[c])-(h[c]*(e[c+1]+(2*e[c]))/3);
		d[c]=(e[c+1]-e[c])/(3*h[c]);
	}
	for(int c=0;c<x_val.size()-1;c++){
		coef[0]=y_val[c];
		coef[1]=b[c];
		coef[2]=e[c];
		coef[3]=d[c];
		coefs.push_back(coef);
	}
/*	for(int c=0;c<coefs.size();c++){
		for(int f=0;f<4;f++)
			cout << coefs[c][f] << " ";
		cout << endl;
	}
*/	return coefs;
}
	
vector<vector<double> > clamped_cubic_spline(vector<double> x_val, vector<double> y_val, double FPO, double FPN){
	double h[x_val.size()], alpha[x_val.size()], a[x_val.size()], b[x_val.size()], c[x_val.size()], d[x_val.size()];
	double l[x_val.size()], mu[x_val.size()], z[x_val.size()];
	vector<double> coef(4);
	vector<vector<double> > coefs;
	for(int i=0;i<x_val.size()-1;i++){
		h[i]=x_val[i+1]-x_val[i];
		a[i]=y_val[i];
	}
	a[x_val.size()-1]=y_val[y_val.size()-1];

	alpha[0]=(3*(a[1]-a[0])/h[0])-(3*FPO);
	alpha[x_val.size()-1]=(3*FPN)-((3*(a[x_val.size()-1]-a[x_val.size()-2]))/h[x_val.size()-2]);
	for(int i=1;i<x_val.size()-1;i++)
		alpha[i]=((3/h[i])*(a[i+1]-a[i]))-((3/h[i-1])*(a[i]-a[i-1]));
	l[0]=2*h[0];
	mu[0]=0.5;
	z[0]=alpha[0]/l[0];
	for(int i=1;i<x_val.size()-1;i++){
		l[i]=(2*(x_val[i+1]-x_val[i-1]))-(h[i-1]*mu[i-1]);
		mu[i]=h[i]/l[i];
		z[i]=(alpha[i]-(h[i-1]*z[i-1]))/l[i];
	}
	l[x_val.size()-1]=h[x_val.size()-2]*(2-mu[x_val.size()-2]);
	z[x_val.size()-1]=(alpha[x_val.size()-1]-(h[x_val.size()-2]*z[x_val.size()-2]))/l[x_val.size()-1];
	c[x_val.size()-1]=z[x_val.size()-1];
	for(int j=x_val.size()-2;j>=0;j--){
		c[j]=z[j]-(mu[j]*c[j+1]);
		b[j]=((a[j+1]-a[j])/h[j])-((h[j]*(c[j+1]+(2*c[j])))/3);
		d[j]=(c[j+1]-c[j])/(3*h[j]);
	}
	for(int i=0;i<x_val.size()-1;i++){
		coef[0]=a[i];
		coef[1]=b[i];
		coef[2]=c[i];
		coef[3]=d[i];
		coefs.push_back(coef);
		for(int j=0;j<4;j++)
			cout << coef[j] << " ";
		cout << endl;
	}
	return coefs;
}

vector<vector<double> > bezier(vector<double> x_end, vector<double> y_end, vector<double> x_l_guide, vector<double> y_l_guide, vector<double> x_r_guide, vector<double> y_r_guide){
	vector<vector<double> > coefs;
	vector<double> coef(8);
	for(int i=0;i<x_end.size()-1;i++){
		coef[0]=x_end[i];
		coef[4]=y_end[i];
		coef[1]=3*(x_l_guide[i]-x_end[i]);
		coef[5]=3*(y_l_guide[i]-y_end[i]);
		coef[2]=3*(x_end[i]+x_r_guide[i+1]-(2*x_l_guide[i]));
		coef[6]=3*(y_end[i]+y_r_guide[i+1]-(2*y_l_guide[i]));
		coef[3]=x_end[i+1]-x_end[i]+(3*x_l_guide[i])-(3*x_r_guide[i+1]);
		coef[7]=y_end[i+1]-y_end[i]+(3*y_l_guide[i])-(3*y_r_guide[i+1]);
		coefs.push_back(coef);
	}
	for(int i=0;i<coefs.size();i++){
		cout << "x_i(t) = ";
		for(int j=0;j<4;j++)
			cout << coefs[i][j] << " ";
		cout << endl << "y_i(t) = ";
		for(int j=4;j<8;j++)
			cout << coefs[i][j] << " ";
		cout << endl;
	}
	return coefs;
}

void graph_bezier_curve(vector<vector<double> > coefs){
	ofstream output("bezier.dat");
	double x,y;
	for(int i=0;i<coefs.size();i++){
		for(double d=0;d<=1;d+=0.001){
			x=0;
			y=0;
			for(int j=0;j<4;j++)
				x+=coefs[i][j]*pow(d,j);
			for(int j=4;j<8;j++)
				y+=coefs[i][j]*pow(d,j-4);
			output << x << " " << y << endl;
		}
	}
	output.close();
}

vector<double> three_point(vector<double> x_val, vector<double> y_val){
	vector<double> y_prime;
	double h;
	for(int i=0;i<x_val.size();i++){
		if(i==0 || i==x_val.size()-1){
			if(i==0){
				h=x_val[1]-x_val[0];
				y_prime.push_back(1/(2*h)*((-3*y_val[0])+(4*y_val[1])-y_val[2]));
			}
			else{
				h=x_val[x_val.size()-2]-x_val[x_val.size()-1];
				y_prime.push_back(1/(2*h)*((-3*y_val[x_val.size()-1])+(4*y_val[x_val.size()-2])-y_val[x_val.size()-3]));
			}
		}
		else{
			h=abs(x_val[i+1]-x_val[i]);
			y_prime.push_back(1/(2*h)*(y_val[i+1]-y_val[i-1]));
		}
	}
	return y_prime;
}


double richardson_forward_difference(double *N_1, int size){
	double **table = new double*[size];
	for(int c=0;c<size;c++)
		table[c] = new double[size];
	for(int c=0;c<size;c++)
		table[c][0]=N_1[c];
	for(int c=1;c<size;c++){
		for(int d=c;d<size;d++)
			table[c][d]=table[c-1][d]+((table[c-1][d]-table[c-1][d-1])/(pow(4,c)-1));
	}
	return table[size-1][size-1];
}

double trapezoidal(double (*f)(double), int n, double a, double b){
	double h=(b-a)/n;
	double integral=0;
	for(double c=a;c<=b-h;c+=h){
		integral+=(h/2)*(f(a)+f(a+h));
		a+=h;
	}
	return integral;
}

//n even
double simpson(double (*f)(double), int n, double a, double b){
	double h=(b-a)/n;
	double integral=0;
	for(double c=a;c<=b-(2*h);c+=2*h){
		integral+=(h/3)*(f(a)+(4*f(a+h))+f(a+(2*h)));
		a+=2*h;
	}
	return integral;
}

// 3|n
double simpson_three_eighths(double (*f)(double), int n, double a, double b){
	double h=(b-a)/n;
	double integral=0;
	for(double c=a;c<=b-(3*h);c+=3*h){
		integral+=(3*h/8)*(f(a)+(3*f(a+h))+(3*f(a+(2*h)))+f(a+(3*h)));
		a+=3*h;
	}
	return integral;
}

//4|n
double four_rule(double (*f)(double), int n, double a, double b){
	double h=(b-a)/n;
	double integral=0;
	for(double c=a;c<=b-(4*h);c+=4*h){
		integral+=(2*h/45)*((7*f(a))+(32*f(a+h))+(12*f(a+(2*h)))+(32*f(a+(3*h)))+(7*f(a+(4*h))));
		a+=4*h;
	}
	return integral;
}

//n=0
double cotes_midpoint(double (*f)(double), int n, double a, double b){
	double h=(b-a)/(n+2);
	double integral=0;
	for(double c=a+h;c<=b-h;c+=h){
		integral+=2*h*f(a+h);
		a+=h;
	}
	return integral;
}

double cotes_one(double (*f)(double), int n, double a, double b){
	double h=(b-a)/(n+2);
	double integral=0;
	for(double c=a+h;c<=b-(2*h);c+=2*h){
		integral+=(3*h/2)*(f(a+h)+f(a+(2*h)));
		a+=2*h;
	}
	return integral;
}

//n even
double cotes_two(double (*f)(double), int n, double a, double b){
	double h=(b-a)/(n+2);
	double integral=0;
	for(double c=a+h;c<=b-(3*h);c+=3*h){
		integral+=(4*h/3)*((2*f(a+h))-f(a+(2*h))+(2*f(a+(3*h))));
		a+=3*h;
	}
	return integral;
}

//3|n
double cotes_three(double (*f)(double), int n, double a, double b){
	double h=(b-a)/(n+2);
	double integral=0;
	for(double c=a+h;c<=b-(4*h);c+=4*h){
		integral+=(5*h/24)*((11*f(a+h))+f(a+(2*h))+f(a+(3*h))+(11*f(a+(4*h))));
		a+=4*h;
	}
	return integral;
}

double composite_trapezoidal(double (*f)(double), int n, double a, double b){
	double h=(b-a)/n;
	double integral=f(a)+f(b);
	for(int j=1;j<n;j++)
		integral+=2*f(a+((double)j*h));
	return h*integral/2;
}
double composite_simpson(double (*f)(double), int n, double a, double b){
	double h=(b-a)/n;
	double XI0=f(a)+f(b);
	double XI1=0;
	double XI2=0;
	for(double i=1;i<=n-1;i++){
		if((int)i%2==0)
			XI2+=f(a+(i*h));
		else
			XI1+=f(a+(i*h));
	}
	return h*(XI0+(2*XI2)+(4*XI1))/3;
}

double composite_midpoint(double (*f)(double), int n, double a, double b){
	double h=(b-a)/(n+2);
	double integral=0;
	for(int j=0;j<=n/2;j++)
		integral+=f(a+(((double)(2*j)+1)*h));
	return 2*h*integral;
}

double romberg(double (*f)(double), int n, double a, double b){
	double h=b-a;
	double table[n][n];
	table[0][0]=h*(f(a)+f(b))/2;
	for(int c=1;c<n;c++){
		h/=2;
		table[c][0]=table[c-1][0];
		for(double d=1;d<=pow(2,c-1);d++)
			table[c][0]+=2*h*f(a+((2*d-1)*h));
		table[c][0]/=2;
	}
	for(int c=1;c<n;c++){
		for(int d=c;d<n;d++)
			table[d][c]=table[d][c-1]+((table[d][c-1]-table[d-1][c-1])/(pow(4,c)-1));
	}
	return table[n-1][n-1];
}

double quadrature(double (*f)(double), double TOL, int max_it, double a, double b){
	double app=0;
	int i=0;
	vector<double> tol, a_i, h, FA, FC, FB, S, L, v(8);
	double FD, FE, S1, S2;
	tol.push_back(10*TOL);
	a_i.push_back(a);
	h.push_back((b-a)/2);
	FA.push_back(f(a));
	FC.push_back(f(a+h[i]));
	FB.push_back(f(b));
	S.push_back(h[i]*(FA[i]+(4*FC[i])+FB[i])/3);
	L.push_back(1);
	while(i>=0){
		FD=f(a_i[i]+(h[i]/2));
		FE=f(a_i[i]+(3*h[i]/2));
		S1=h[i]*(FC[i]+(4*FE)+FB[i])/6;
		S2=h[i]*(FA[i]+(4*FD)+FC[i])/6;
		v[0]=(a_i[i]);
		v[1]=(FA[i]);
		v[2]=(FC[i]);
		v[3]=(FB[i]);
		v[4]=h[i];
		v[5]=tol[i];
		v[6]=S[i];
		v[7]=L[i];
		i--;
		tol.pop_back();
		a_i.pop_back();
		h.pop_back();
		FA.pop_back();
		FC.pop_back();
		FB.pop_back();
		S.pop_back();
		L.pop_back();
		if(abs(S1+S2-v[6])<v[5])
			app+=(S1+S2);
		else{
			if(v[7]>=max_it){
				cout << "LEVEL EXCEEDED\n";
				return pow(-1,.5);
			}
			else{
				i++;
				a_i.push_back(v[0]+v[4]);
				FA.push_back(v[2]);
				FC.push_back(FE);
				FB.push_back(v[3]);
				h.push_back(v[4]/2);
				tol.push_back(v[5]/2);
				S.push_back(S2);
				L.push_back(v[7]+1);
				i++;
				a_i.push_back(v[0]);
				FA.push_back(v[1]);
				FC.push_back(FD);
				FB.push_back(v[2]);
				h.push_back(h[i-1]);
				tol.push_back(tol[i-1]);
				S.push_back(S1);
				L.push_back(L[i-1]);
			}
		}
	}
	return app;
}
	
double gauss_quad(double (*f)(double), double a, double b, vector<double> coefs, vector<double> roots){
	double integral=0;
	for(int c=0;c<coefs.size();c++)
		integral+=coefs[c]*f((((b-a)*roots[c])+b+a)/2)*((b-a)/2);
	return integral;
}

double simpson_double_integral(double (*f)(double,double), double a, double b, double (*c)(double), double (*d)(double), int m, int n){
	double h=(b-a)/(double)n;
	double J_1=0;
	double J_2=0;
	double J_3=0;
	double x, y, Q, L;
	double HX, K_1, K_2, K_3;
	for(int i=0;i<=n;i++){
		x=a+((double)i*h);
		HX=(d(x)-c(x))/m;
		K_1=f(x,c(x))+f(x,d(x));
		K_2=0;
		K_3=0;
		for(int j=1;j<m;j++){
			y=c(x)+((double)j*HX);
			Q=f(x,y);
			if(j%2==0)
				K_2+=Q;
			else
				K_3+=Q;
		}
		L=(HX*(K_1+(2*K_2)+(4*K_3)))/3;
		if(i==0 || i==n)
			J_1+=L;
		else if(i%2==0)
			J_2+=L;
		else
			J_3+=L;
	}
	return (h*(J_1+(2*J_2)+(4*J_3)))/3;
}

double second_midpoint(double (*f)(double),double x, double h){
	return (f(x-h)-(2*f(x))+f(x+h))/pow(h,2);
}

