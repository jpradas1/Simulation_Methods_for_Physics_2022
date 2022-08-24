#include <iostream>
#include <cmath>

const double ErrMax = 1e-15;

double f(double t, double x){
    return x;
}

void RK4_Step(double & t0, double & x0, double dt){
    double dx1, dx2, dx3, dx4;
    dx1 = dt * f(t0,x0);
    dx2 = dt * f(t0+dt/2, x0+dx1/2);
    dx3 = dt * f(t0+dt/2, x0+dx2/2);
    dx4 = dt * f(t0+dt, x0+dx3);

    x0 += (dx1 + 2*dx2 + 2*dx3 + dx4)/6 ; t0+= dt;
}

int main(){
    double t, x , dt = 0.01;

    for(t=0, x=1 ; t <= 2 + dt ; ){
        std::cout << t << "\t" << x << "\t" << exp(t) << std::endl;
        RK4_Step(t, x, dt);
    }

    return 0;
}