#include <iostream>
#include <cmath>

#include "./../Requirements/vector_old.h"

const double G = 1.0;
const int N = 2;

const double Zeta = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e0;
const double Chi = -0.6626458266981849e-1;

const double Coeff1 = (1-2*Lambda)/2;
const double Coeff2 = 1 - 2 * (Chi + Zeta);

class Body{
    public:
        void init(double x0, double y0, double z0, double Vx0, 
                  double Vy0, double Vz0, double m0, double R0);
        void delete_forces();
        void add_forces(vector3D F0);
        void move_r(double dt, double coeff);
        void move_v(double dt, double coeff);

        double get_x(){return r.x();};
        double get_y(){return r.y();};
        double get_z(){return r.z();};

        friend class Collider;
    private:
        vector3D r, V, F;
        double m, R;
};

void Body::init(double x0, double y0, double z0, double Vx0, 
                double Vy0, double Vz0, double m0, double R0){
    r.load(x0, y0, z0); V.load(Vx0, Vy0, Vz0);
    m= m0; R= R0;
}

void Body::delete_forces(){
    F.load(0, 0, 0);
}

void Body::add_forces(vector3D F0){
    F += F0;
}

void Body::move_r(double dt, double coeff){
    r += V * (dt * coeff);
}

void Body::move_v(double dt, double coeff){
    V += F * (dt * coeff/m);
}

class Collider {
    public:
        void compute_forces(Body * Planet);
        void compute_forces_between(Body & Planet1, Body & Planet2);
    private:

};

void Collider::compute_forces(Body * Planet){
    int i, j;

    for(i=0; i<N; i++)
        Planet[i].delete_forces();

    for(i=0; i<N; i++)
        for(j=i+1; j<N; j++)
            compute_forces_between(Planet[i], Planet[j]);
}

void Collider::compute_forces_between(Body & Planet1, Body & Planet2){
    vector3D r21, n, F1;  double d21, F;
    r21 = Planet2.r - Planet1.r; d21 = r21.norm();
    n = r21/d21;
    F = G * Planet1.m  * Planet2.m * pow( d21, -2.0);
    F1 = F * n; 
    Planet1.add_forces(F1); Planet2.add_forces(F1 * (-1)); 
}

int main(){
    double t, dt=1.0; int ii;
    double omega, T, V0, V1; 
    double m0 = 10, m1= 1, r = 11;
    double M = m0 + m1, x0 = -m1*r /M, x1 = m0*r /M;
    Body Planet[N];
    Collider Newton;

    omega = sqrt(G*M*pow(r,-3));
    V0 = omega*x0; V1 = omega * x1;
    T = 2*M_PI / omega;

    // double x0, double y0, double z0, double Vx0, 
    // double Vy0, double Vz0, double m0, double R0

    Planet[0].init(x0, 0 , 0, 0, V0, 0, m0, 0.5);
    Planet[1].init(x1, 0 , 0, 0, V1, 0, m1, 0.5);

    for(t=0; t<1.1*T; t+= dt){
        std::cout << Planet[1].get_x() << "\t" << Planet[1].get_y() << std::endl;
        for(ii = 0; ii < N; ii++) Planet[ii].move_r(dt, Zeta);

        Newton.compute_forces(Planet);
        for(ii = 0; ii < N; ii++) Planet[ii].move_v(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Planet[ii].move_r(dt, Chi);

        Newton.compute_forces(Planet);
        for(ii = 0; ii < N; ii++) Planet[ii].move_v(dt, Lambda);
        for(ii = 0; ii < N; ii++) Planet[ii].move_r(dt, Coeff2);

        Newton.compute_forces(Planet);
        for(ii = 0; ii < N; ii++) Planet[ii].move_v(dt, Lambda);
        for(ii = 0; ii < N; ii++) Planet[ii].move_r(dt, Chi);

        Newton.compute_forces(Planet);
        for(ii = 0; ii < N; ii++) Planet[ii].move_v(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Planet[ii].move_r(dt, Zeta);
    }

    return 0;
}