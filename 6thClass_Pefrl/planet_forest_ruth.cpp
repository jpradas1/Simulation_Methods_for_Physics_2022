#include <iostream>
#include <cmath>

#include "./../Requirements/vector.h"

const double GM = 1.0;

const double Theta = 1 /(2 - pow(2.0, 1/3.0));
const double coeff1 = Theta/2;
const double coeff2 = (1- Theta)/2;
const double coeff3 = 1 - 2 * Theta;

class Body{
    public:
        void init(double x0, double y0, double z0, double Vx0, 
                  double Vy0, double Vz0, double m0, double R0);
        void compute_Forces();
        void move_r(double dt, double coeff);
        void move_v(double dt, double coeff);

        double get_x(){return r.x();};
        double get_y(){return r.y();};
        double get_z(){return r.z();};
    private:
        vector3D r, V, F;
        double m, R;
};

void Body::init(double x0, double y0, double z0, double Vx0, 
                double Vy0, double Vz0, double m0, double R0){
    r.load(x0, y0, z0); V.load(Vx0, Vy0, Vz0);
    m= m0; R= R0;
}

void Body::compute_Forces(){
    double F_aux = GM * m/(pow(r.norm(), 3));
    F = (-F_aux) * r;
}

void Body::move_r(double dt, double coeff){
    r += V * (dt * coeff);
}

void Body::move_v(double dt, double coeff){
    V += F * (dt * coeff/m);
}

int main(){
    double t, dt=1.0;
    double omega, T, V0, r0=100; 
    double m = 1;
    Body Planet;

    omega = sqrt(GM*pow(r0,-3));
    V0 = omega*r0;
    T = 2*M_PI / omega;

    // double x0, double y0, double z0, double Vx0, 
    // double Vy0, double Vz0, double m0, double R0

    Planet.init(r0, 0 , 0, 0, V0/2, 0, m, 0.15);

    for(t=0; t<1.1*T; t+= dt){
        std::cout << Planet.get_x() << "\t" << Planet.get_y() << std::endl;
        Planet.move_r(dt, coeff1);

        Planet.compute_Forces();
        Planet.move_v(dt, Theta);

        Planet.move_r(dt, coeff2);

        Planet.compute_Forces();
        Planet.move_v(dt, coeff3);

        Planet.move_r(dt, coeff2);

        Planet.compute_Forces();
        Planet.move_v(dt, Theta);

        Planet.move_r(dt, coeff1);
    }

    return 0;
}