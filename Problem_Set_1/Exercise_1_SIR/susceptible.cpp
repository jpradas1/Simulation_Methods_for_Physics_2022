#include <iostream>
#include <cmath>

// const double BETA = 0.35;
const double GAMMA = 0.08;

double f1(double t, double ss, double ii, double rr, double BETA){
    return -BETA * ss * ii;
}

double f2(double t, double ss, double ii, double rr, double BETA){
    return BETA * ss * ii - GAMMA * ii;
}

double f3(double t, double ss, double ii, double rr, double BETA){
    return GAMMA * ii;
}

template <class T>
void Runge_Kutta_Coupled(double dt ,double & t0, double & x10, double & x20,
                         double & x30, T g1, T g2, T g3, double BETA){
    double dx11, dx21, dx31, dx41;
    double dx12, dx22, dx32, dx42;
    double dx13, dx23, dx33, dx43;

    dx11 = dt * g1(t0,x10,x20,x30, BETA);
    dx21 = dt * g1(t0+dt/2, x10+dx11/2, x20+dx12/2, x30+dx13/2, BETA);
    dx31 = dt * g1(t0+dt/2, x10+dx21/2, x20+dx22/2, x30+dx23/2, BETA);
    dx41 = dt * g1(t0+dt, x10+dx31, x20+dx32, x30+dx33, BETA);

    dx12 = dt * g2(t0,x10,x20, x30, BETA);
    dx22 = dt * g2(t0+dt/2, x10+dx11/2, x20+dx12/2, x30+dx13/2, BETA);
    dx32 = dt * g2(t0+dt/2, x10+dx21/2, x20+dx22/2, x30+dx23/2, BETA);
    dx42 = dt * g2(t0+dt, x10+dx31, x20+dx32, x30+dx33, BETA);

    dx13 = dt * g3(t0,x10,x20, x30, BETA);
    dx23 = dt * g3(t0+dt/2, x10+dx11/2, x20+dx12/2, x30+dx13/2, BETA);
    dx33 = dt * g3(t0+dt/2, x10+dx21/2, x20+dx22/2, x30+dx23/2, BETA);
    dx43 = dt * g3(t0+dt, x10+dx31, x20+dx32, x30+dx33, BETA);

    t0+= dt;
    x10 += (dx11 + 2*dx21 + 2*dx31 + dx41)/6;
    x20 += (dx12 + 2*dx22 + 2*dx32 + dx42)/6;
    x30 += (dx13 + 2*dx23 + 2*dx33 + dx43)/6;
}


int main(){

    double t, ss, ii , rr = 0, dt = 0.1;
    double BETA; double max = 0.40;

    for(BETA = 0; BETA <= max; BETA+=0.01){
        for(t=0, ss=0.999, ii=0.001 ; t <= 70 ; ){
        Runge_Kutta_Coupled(dt, t, ss, ii, rr, f1, f2, f3, BETA);
        }
        std::cout << BETA/GAMMA << "\t" << ss << std::endl;
    }

    return 0;
}