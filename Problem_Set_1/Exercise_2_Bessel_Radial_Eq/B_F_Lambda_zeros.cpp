#include <iostream>
#include <cmath>

// const double Lambda = 1;
const double ErrMax = 1e-12;

double f1(double r, double x1, double x2, double Lambda){
    return -pow(Lambda,2) * x2 ;
}

double f2(double r, double x1, double x2, double Lambda){
    return x1 + x2/r;
}

template <class T>
void Runge_Kutta_Coupled(double dt ,double & t0, double & x10, 
                           double & x20, T g1, T g2, double Lambda){
    double dx11, dx21, dx31, dx41;
    double dx12, dx22, dx32, dx42;

    dx11 = dt * g1(t0,x10,x20, Lambda);
    dx21 = dt * g1(t0+dt/2, x10+dx11/2, x20+dx12/2, Lambda);
    dx31 = dt * g1(t0+dt/2, x10+dx21/2, x20+dx22/2, Lambda);
    dx41 = dt * g1(t0+dt, x10+dx31, x20+dx32, Lambda);

    dx12 = dt * g2(t0,x10,x20, Lambda);
    dx22 = dt * g2(t0+dt/2, x10+dx11/2, x20+dx12/2, Lambda);
    dx32 = dt * g2(t0+dt/2, x10+dx21/2, x20+dx22/2, Lambda);
    dx42 = dt * g2(t0+dt, x10+dx31, x20+dx32, Lambda);

    t0+= dt;
    x10 += (dx11 + 2*dx21 + 2*dx31 + dx41)/6;
    x20 += (dx12 + 2*dx22 + 2*dx32 + dx42)/6;
}

double R(double r, double Lambda){
    double r_aux, x1, x2 , dr = 0.01;

    for(r_aux=0.01, x1 = 0*r_aux , x2=1*r_aux ; r_aux <= r ; ){
        Runge_Kutta_Coupled(dr, r_aux, x1, x2, f1, f2, Lambda);
    }

    return x2/r;
}

template <class T>
double Zeros_by_bisection(double a, double b, T g, double r){
    double m, ga, gm;

    ga = g(r,a);

    while(b-a > ErrMax){
        m = (b+a)/2, gm = g(r, m);
        if(ga*gm>0)
            {a = m; ga = gm;}
        else
            b = m;

    }

    return (a+b)/2;
}


int main(){
    double Lambda=1, r=1;
    double A[] = {1.5, 4, 6, 10, 13};
    int n = sizeof(A)/sizeof(A[0]);

    for(int ii = 0; ii<n-1; ii++)
        std::cout << Zeros_by_bisection(A[ii], A[ii+1], R, r) << "\t";

    // for(Lambda=0.1; Lambda <= 15.0; Lambda+=0.1){
    //     std::cout << Lambda << "\t" << R(r, Lambda) << std::endl;
    // }

    return 0;
}