#include <iostream>
#include <cmath>

#include "./../Requirements/vector_old.h"
#include "./../Requirements/Random64.h"

const double g = 9.8;
const int Nx = 5, Ny = 5, N = Nx * Ny;
const double K = 1e4;
const double Lx = 60, Ly = 60;

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
                  double m0, double R0);
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
};

void Body::init(double x0, double y0, double Vx0, double Vy0, 
                double m0, double R0){
    r.load(x0, y0, 0); V.load(Vx0, Vy0, 0);
    m= m0; R= R0;
}

void Body::move_r(double dt, double coeff){
    r += V * (dt * coeff);
}

void Body::move_v(double dt, double coeff){
    V += F * (dt * coeff/m);
}

void Body::Draw(){
    std::cout << " , " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t)";
}

class Collider {
    public:
        void compute_forces(Body * Molecule);
        void compute_forces_between(Body & Molecule1, Body & Molecule2);
    private:
};

void Collider::compute_forces(Body * Molecule){
    int i, j;

    for(i=0; i<N+4; i++)
        Molecule[i].delete_forces();

    for(i=0; i<N; i++)
        for(j=i+1; j<N+4; j++)
            compute_forces_between(Molecule[i], Molecule[j]);
}

void Collider::compute_forces_between(Body & Molecule1, Body & Molecule2){
    vector3D r21, n, F1;  double d21, F, s;
    r21 = Molecule2.r - Molecule1.r; d21 = r21.norm();
    s = Molecule2.R + Molecule1.R - d21;
    if (s > 0){
        n = r21/d21; F = K * pow(s, 1.5);
        F1 = n * F;
        Molecule2.add_forces(F1);
        Molecule1.add_forces(F1*(-1));
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
    Body Molecule[N+4];
    Collider Hertz;
    Crandom rand64(1);
    double m0 = 1, R0 = 2, kT = 10, theta;
    double V0 = sqrt(2 * kT / m0);
    double t, tmax = 100*(Lx/V0), dt=1e-3;
    double tdraw, tdomain = tmax/1000;
    int ii, ix, iy;
    double dx = Lx /(Nx+1), dy = Ly /(Ny+1);

    // double x0, double y0, double Vx0, 
    // double Vy0, double m0, double R0

    StartAnimation();

    double RWall = 100 * Lx, MWall = 10000 * m0;

    Molecule[N+0].init(Lx/2, Ly+RWall, 0, 0, MWall, RWall); // Up Wall
    Molecule[N+1].init(Lx/2, -RWall, 0, 0, MWall, RWall); // Down Wall
    Molecule[N+2].init(Lx+RWall, Ly/2, 0, 0, MWall, RWall); // Right Wall
    Molecule[N+3].init(-RWall, Ly/2, 0, 0, MWall, RWall); //  Left Wall

    for(ix = 0; ix < Nx; ix++)
        for(iy = 0; iy < Ny; iy++){
            theta = 2 * M_PI * rand64.r();
            Molecule[Nx*iy+ix].init((ix+1)*dx, (iy+1)*dy, V0*cos(theta), V0*sin(theta), m0, R0);
        }

    for(t=0 , tdraw = 0; t < tmax; t+= dt , tdraw+= dt){
        if(tdraw > tdomain){
            StartDomain();
            for(ii=0 ;ii<N; ii++) Molecule[ii].Draw();
            FinishDomain();
            tdraw = 0;
        }
        // std::cout << Molecule[1].get_x() << "\t" << Molecule[1].get_y() << std::endl;
        for(ii = 0; ii < N; ii++) Molecule[ii].move_r(dt, Zeta);

        Hertz.compute_forces(Molecule);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_v(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_r(dt, Chi);

        Hertz.compute_forces(Molecule);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_v(dt, Lambda);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_r(dt, Coeff2);

        Hertz.compute_forces(Molecule);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_v(dt, Lambda);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_r(dt, Chi);

        Hertz.compute_forces(Molecule);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_v(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Molecule[ii].move_r(dt, Zeta);
    }

    return 0;
}