#include <iostream>
#include <cmath>

#include "./../Requirements/vector.h"
#include "./../Requirements/Random64.h"

const double g = 9.8;
const int Nx = 1, Ny = 1, N = Nx * Ny;
const double K = 1e4;
const double Lx = 60, Ly = 60;
const double Gamma = 10;

const double Zeta = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e0;
const double Chi = -0.6626458266981849e-1;

const double Coeff1 = (1-2*Lambda)/2;
const double Coeff2 = 1 - 2 * (Chi + Zeta);

class Body;
class Collider;

class Body{
    public:
        void init(double x0, double y0, double Vx0, double Vy0, 
                  double theta0, double omega0, double m0, double R0);
        void delete_forces() {F.load(0, 0, 0);};
        void add_forces(vector3D F0) {F += F0;};
        void move_r(double dt, double coeff);
        void move_v(double dt, double coeff);

        double get_x(){return r.x();};
        double get_y(){return r.y();};

        void Draw();

        friend class Collider;
    private:
        vector3D r, V, F;
        double m, R;
        double theta, omega, tau, I;
};

void Body::init(double x0, double y0, double Vx0, double Vy0, 
                double theta0, double omega0, double m0, double R0){
    r.load(x0, y0, 0); V.load(Vx0, Vy0, 0);
    m= m0; R= R0;
    theta = theta0; omega = omega0; 
    I = 2.0/5 * m * R * R; 
}

void Body::move_r(double dt, double coeff){
    r += V * (dt * coeff);
    theta += omega + (dt * coeff);
}

void Body::move_v(double dt, double coeff){
    V += F * (dt * coeff/m);
    omega += tau * (dt * coeff/I);
}

void Body::Draw(){
    std::cout << " , " << r.x() << "+" << R << "*cos(t),"
              << r.y() << "+" << R << "*sin(t), "
              << r.x() << "+" << R*cos(theta)/7.0 << "*t,"
              << r.y() << "+" << R*sin(theta)/7.0 << "*t";
}

class Collider {
    public:
        void compute_forces(Body * Grain);
        void compute_forces_between(Body & Grain1, Body & Grain2);
    private:
};

void Collider::compute_forces(Body * Grain){
    int i, j; vector3D W; double w;

    for(i=0; i<N+4; i++)
        Grain[i].delete_forces();

    for(i=0; i<N ; i++){
        w = Grain[i].m * g;
        W.load(0, -w, 0);
        Grain[i].add_forces(W);
    }

    for(i=0; i<N; i++)
        for(j=i+1; j<N+4; j++)
            compute_forces_between(Grain[i], Grain[j]);
}

void Collider::compute_forces_between(Body & Grain1, Body & Grain2){
    vector3D r21, n, F1, V21;  
    double d21, F, s, HK, m12;
    r21 = Grain2.r - Grain1.r; d21 = r21.norm();
    s = Grain2.R + Grain1.R - d21;
    if (s > 0){
        n = r21/d21; F = K * pow(s, 1.5); //Hertz's force
        m12 = (Grain1.m * Grain2.m)/(Grain1.m + Grain2.m);
        HK = Gamma*sqrt(s)*m12; // Hertz Kuwabara-Kono's force
        F1 = n * F - HK * V21;
        V21 = Grain2.V - Grain1.V;
        Grain2.add_forces(F1);
        Grain1.add_forces(F1*(-1));
    }
}

void StartAnimation(){
    // std::cout << "set terminal gif animate" << std::endl;
    // std::cout << "set output 'gas_2d.gif'" << std::endl;
    std::cout << "unset key" << std::endl;
    std::cout << "set xrange[-10:"<<Lx+10<<"]" << std::endl;
    std::cout << "set yrange[-10:"<<Ly+10<<"]" << std::endl;
    std::cout << "set size ratio -1" << std::endl;
    std::cout << "set parametric" << std::endl;
    std::cout << "set trange [0:7]" << std::endl;
    std::cout << "set isosamples 12" << std::endl;
}

void StartDomain(){
    std::cout << "plot 0,0 ";
    std::cout<<" , "<<Lx/7<<"*t,0";        
    std::cout<<" , "<<Lx/7<<"*t,"<<Ly;     
    std::cout<<" , 0,"<<Ly/7<<"*t";        
    std::cout<<" , "<<Lx<<","<<Ly/7<<"*t";
}

void FinishDomain(){
    std::cout << std::endl;
}

int main(){
    Body Grain[N+4];
    Collider Hertz;
    Crandom rand64(1);
    double m0 = 1, R0 = 8, kT = 10, theta;
    // double V0 = sqrt(2 * kT / m0);
    double t, tmax = 20*sqrt(Ly/g), dt=1e-3;
    double tdraw, tdomain = tmax/1000;
    int ii, ix, iy;
    double dx = Lx /(Nx+1), dy = Ly /(Ny+1);
    double theta0 = 0, omega0 = 0, omega_max = 1.0*10;

    // double x0, double y0, double Vx0, 
    // double Vy0, double m0, double R0

    StartAnimation();

    double RWall = 100 * Lx, MWall = 10000 * m0;

    Grain[N+0].init(Lx/2, Ly+RWall, 0, 0, theta0, omega0, MWall, RWall); // Up Wall
    Grain[N+1].init(Lx/2, -RWall, 0, 0, theta0, omega0, MWall, RWall); // Down Wall
    Grain[N+2].init(Lx+RWall, Ly/2, 0, 0, theta0, omega0, MWall, RWall); // Right Wall
    Grain[N+3].init(-RWall, Ly/2, 0, 0, theta0, omega0, MWall, RWall); //  Left Wall

    for(ix = 0; ix < Nx; ix++)
        for(iy = 0; iy < Ny; iy++){
            // theta = 2 * M_PI * rand64.r();
            Grain[Nx*iy+ix].init((ix+1)*dx, (iy+1)*dy, 0, 0, theta0, omega_max, m0, R0);
        }

    for(t=0 , tdraw = 0; t < tmax; t+= dt , tdraw+= dt){
        if(tdraw > tdomain){
            StartDomain();
            for(ii=0 ;ii<N; ii++) Grain[ii].Draw();
            FinishDomain();
            tdraw = 0;
        }
        // std::cout << Grain[1].get_x() << "\t" << Grain[1].get_y() << std::endl;
        for(ii = 0; ii < N; ii++) Grain[ii].move_r(dt, Zeta);

        Hertz.compute_forces(Grain);
        for(ii = 0; ii < N; ii++) Grain[ii].move_v(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Grain[ii].move_r(dt, Chi);

        Hertz.compute_forces(Grain);
        for(ii = 0; ii < N; ii++) Grain[ii].move_v(dt, Lambda);
        for(ii = 0; ii < N; ii++) Grain[ii].move_r(dt, Coeff2);

        Hertz.compute_forces(Grain);
        for(ii = 0; ii < N; ii++) Grain[ii].move_v(dt, Lambda);
        for(ii = 0; ii < N; ii++) Grain[ii].move_r(dt, Chi);

        Hertz.compute_forces(Grain);
        for(ii = 0; ii < N; ii++) Grain[ii].move_v(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Grain[ii].move_r(dt, Zeta);
    }

    return 0;
}