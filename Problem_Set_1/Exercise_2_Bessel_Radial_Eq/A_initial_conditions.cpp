#include <iostream>
#include <cmath>

const double Lambda = 1;

double f1(double r, double x1, double x2){
    return -pow(Lambda,2) * x2 ;
}

double f2(double r, double x1, double x2){
    return x1 + x2/r;
}

template <class T>
void Runge_Kutta_Coupled(double dt ,double & t0, double & x10, 
                           double & x20, T g1, T g2){
    double dx11, dx21, dx31, dx41;
    double dx12, dx22, dx32, dx42;

    dx11 = dt * g1(t0,x10,x20);
    dx21 = dt * g1(t0+dt/2, x10+dx11/2, x20+dx12/2);
    dx31 = dt * g1(t0+dt/2, x10+dx21/2, x20+dx22/2);
    dx41 = dt * g1(t0+dt, x10+dx31, x20+dx32);

    dx12 = dt * g2(t0,x10,x20);
    dx22 = dt * g2(t0+dt/2, x10+dx11/2, x20+dx12/2);
    dx32 = dt * g2(t0+dt/2, x10+dx21/2, x20+dx22/2);
    dx42 = dt * g2(t0+dt, x10+dx31, x20+dx32);

    t0+= dt;
    x10 += (dx11 + 2*dx21 + 2*dx31 + dx41)/6;
    x20 += (dx12 + 2*dx22 + 2*dx32 + dx42)/6;
}


int main(){

    double r, x1, x2 , dr = 0.01;
    double a = 1, r_max=10.0;

    for(r=0.01, x1 = 0*r , x2=1*r ; r <= r_max ; ){
        std::cout << r << "\t" << x2/r << std::endl;
        Runge_Kutta_Coupled(dr, r, x1, x2, f1, f2);
    }

    return 0;
}