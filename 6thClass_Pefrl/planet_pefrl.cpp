#include <iostream>
#include <cmath>

#include "./../Requirements/vector.h"

const double GM = 1.0;

const double Zeta = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e0;
const double Chi = -0.6626458266981849e-1;

const double Coeff1 = (1-2*Lambda)/2;
const double Coeff2 = 1 - 2 * (Chi + Lambda);

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
    double omega, T, V0, r0=10; 
    double m = 1;
    Body Planet;

    omega = sqrt(GM*pow(r0,-3));
    V0 = omega*r0;
    T = 2*M_PI / omega;

    // double x0, double y0, double z0, double Vx0, 
    // double Vy0, double Vz0, double m0, double R0

    Planet.init(r0, 0 , 0, 0, V0/2, 0, m, 0.5);

    for(t=0; t<1.1*T; t+= dt){
        std::cout << Planet.get_x() << "\t" << Planet.get_y() << std::endl;
        Planet.move_r(dt, Zeta);

        Planet.compute_Forces();
        Planet.move_v(dt, Coeff1);

        Planet.move_r(dt, Chi);

        Planet.compute_Forces();
        Planet.move_v(dt, Lambda);

        Planet.move_r(dt, Coeff2);

        Planet.compute_Forces();
        Planet.move_v(dt, Lambda);

        Planet.move_r(dt, Chi);

        Planet.compute_Forces();
        Planet.move_v(dt, Coeff1);

        Planet.move_r(dt, Zeta);
    }

    return 0;
}