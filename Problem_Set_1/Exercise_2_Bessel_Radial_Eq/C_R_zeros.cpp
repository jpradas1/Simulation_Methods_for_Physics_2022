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

    for(r_aux=0.01, x1 = 1*r_aux , x2=0*r_aux ; r_aux <= r ; ){
        Runge_Kutta_Coupled(dr, r_aux, x1, x2, f1, f2, Lambda);
    }

    return x2/r;
}

template <class T>
double Zeros_by_bisection(double a, double b, T g, double lambda){
    double m, ga, gm;

    ga = g(a, lambda);

    while(b-a > ErrMax){
        m = (b+a)/2, gm = g(m, lambda);
        if(ga*gm>0)
            {a = m; ga = gm;}
        else
            b = m;

    }

    return (a+b)/2;
}

double F(double lambda){
    double r = 0.01, Raux = 1, a, b;
    bool flag = true;
    double R_p, R1 = R(r, lambda);
    
    while(Raux > 0){
        R_p = R1;
        a = r;
        r += 0.01;
        R1 = R(r, lambda);
        b = r;
        Raux = (R_p * R1)/abs(R_p * R1 );
    }

    r = Zeros_by_bisection(a, b, R, lambda);
    return r;
}

double Lambda(double r){
    double lambda = 2.7, raux = F(lambda);
    while(raux - r > 0){
        lambda += 1e-5;
        raux = F(lambda);
    }
    return lambda;
}


int main(){
    double r = 1;
    double lambda = Lambda(r);

    std::cout.precision(8);
    std::cout.setf(std::ios::scientific);

    for(r=0.001; r<=1.1; r+=0.01)
        std::cout << r << "\t" << R(r, lambda) << "\t" << R(r, 2*lambda)
                  << "\t" << R(r, 3*lambda) << "\t" << R(r, 4*lambda)
                  << "\t" << 0 << std::endl;

    return 0;
}