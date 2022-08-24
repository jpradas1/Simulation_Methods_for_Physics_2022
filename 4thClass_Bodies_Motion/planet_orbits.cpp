#include <iostream>
#include <cmath>

const double g=9.8; // m/s²
// const double AU = 149.6e6 * 1000; // m
// const double G = 6.67428e-11; // N*m² / kg²

const double GM = 1;

class Body{
    public:
        void init(double x0, double y0, double Vx0, double Vy0, 
                  double m0, double R0);
        void compute_Forces();
        void movement(double dt);

        double get_x(){return x;};
        double get_y(){return y;};

        double ratio(){return sqrt(x*x+y*y);};
    private:
        double x, y, Vx, Vy, Fx, Fy, m, R;
};

void Body::init(double x0, double y0, double Vx0, double Vy0, 
                double m0, double R0){
    x= x0; y= y0; Vx= Vx0; Vy= Vy0; m= m0; R= R0;
}

void Body::compute_Forces(){
    double F = GM * m/(pow(ratio(), 3));
    Fx= -F * x; Fy= -F * y;
}

void Body::movement(double dt){
    x+= Vx * dt;    y+= Vy * dt;
    Vx+= Fx/m*dt;   Vy+= Fy/m*dt;
}

int main(){
    double t, dt=0.001;
    double omega, T, V0, r0=100; 
    double m = 1;
    Body Planet;

    omega = sqrt(GM*pow(r0,-3));
    V0 = omega*r0;
    T = 2*M_PI / omega;

    // double x0, double y0, double Vx0, double Vy0, 
    // double m0, double R0

    Planet.init(r0, 0 , 0, V0/2, m, 0.15);

    for(t=0; t<1.1*T; t+= dt){
        std::cout << Planet.get_x() << "\t" << Planet.get_y() << std::endl;
        Planet.compute_Forces();
        Planet.movement(dt);
    }

    return 0;
}