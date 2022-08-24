#include <iostream>
#include <cmath>

const double ErrMax = 1e-15;

double f(double t, double x){
    return x;
}

void Euler_Step(double & t, double & x, double dt){
    double dx;
    dx = dt * f(t,x);
    x += dx ; t+= dt;
}

int main(){
    double t, x , dt = 0.01;

    for(t=0, x=1 ; t <= 2 + dt ; ){
        std::cout << t << "\t" << x << "\t" << exp(t) << std::endl;
        Euler_Step(t, x, dt);
    }

    return 0;
}