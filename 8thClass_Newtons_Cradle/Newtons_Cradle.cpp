#include <iostream>
#include <cmath>

#include "./../Requirements/vector_old.h"

const double K = 1e9;
const int N = 5;
const double g = 980;

const double Zeta = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e0;
const double Chi = -0.6626458266981849e-1;

const double Coeff1 = (1-2*Lambda)/2;
const double Coeff2 = 1 - 2 * (Chi + Zeta);

class Body;
class Collider;

class Body{
    public:
        void init(double theta0, double omega0, double m0, 
                  double R0, double l0, double x0);
        void delete_torsion() {tau = 0;};
        void add_torsion(double tau0) {tau += tau0;};
        void move_theta(double dt, double coeff);
        void move_omega(double dt, double coeff);

        double get_x(){return x0+l*sin(theta);};
        double get_y(){return -l*cos(theta);};
        double get_tau(){return tau;};

        void Draw();

        friend class Collider;
    private:
        double theta, omega, tau;
        double m, R, l, x0, I;
};

void Body::init(double theta0, double omega0, double m0, 
                double R0, double l0, double x00){
    theta = theta0; omega = omega0; m= m0; R= R0;
    l = l0; x0 = x00; I = m*l*l;
}

void Body::move_theta(double dt, double coeff){
    theta += omega * (dt * coeff);
}

void Body::move_omega(double dt, double coeff){
    omega += tau * (dt * coeff/I);
}

void Body::Draw(){
    std::cout << " , " << get_x() << "+" << R << "*cos(t)," << get_y() << "+" << R << "*sin(t)";
    std::cout << " , " << x0 << "+"<<l/7<<"*t*sin("<< theta << "),-" << l/7 << "*t*cos(" << theta<< ")";
}

class Collider {
    public:
        void compute_torsion(Body * Pendulum);
        void compute_torsion_between(Body & Pendulum1, Body & Pendulum2);
    private:
};

void Collider::compute_torsion(Body * Pendulum){
    int i; double tau0;

    for(i=0; i<N; i++){
        Pendulum[i].delete_torsion();
        tau0 = -Pendulum[i].l * Pendulum[i].m * g * sin(Pendulum[i].theta);
        Pendulum[i].add_torsion(tau0);
    }

    for(i=N-1; i>0; i--)
        compute_torsion_between(Pendulum[i], Pendulum[i-1]);
}

void Collider::compute_torsion_between(Body & Pendulum1, Body & Pendulum2){
    double F=0, s;
    s = (Pendulum2.get_x() + Pendulum2.R)-(Pendulum1.get_x()-Pendulum1.R);
    if (s>0) F = K * pow(s, 1.5);  // Hertz Force
    Pendulum1.add_torsion(F*Pendulum1.l);
    Pendulum2.add_torsion(-1*F*Pendulum2.l);
}

void StartAnimation(){
    // std::cout << "set terminal gif animate" << std::endl;
    // std::cout << "set output 'Two_Pendulums.gif'" << std::endl;
    std::cout << "unset key" << std::endl;
    std::cout << "set xrange[-5:20]" << std::endl;
    std::cout << "set yrange[-15:2]" << std::endl;
    std::cout << "set size ratio -1" << std::endl;
    std::cout << "set parametric" << std::endl;
    std::cout << "set trange [0:7]" << std::endl;
    std::cout << "set isosamples 12" << std::endl;
}

void StartDomain(){
    std::cout << "plot 0,0 ";
}

void FinishDomain(){
    std::cout << std::endl;
}

int main(){
    Body Pendulum[N];
    Collider Newton;
    double m0 = 50, l0 = 12, R0 = 2;
    double theta0 = - M_PI / 30; double x0 = 0; double omega0 = 0;
    double omega, T;
    omega = sqrt(g / l0);
    T = 2*M_PI / omega;
    double t, tmax = 11*T, dt=0.0001;
    double tdraw, tdomain = T/500;
    int ii;

    // double theta0, double omega0, double m0, 
    // double R0, double l0, double x0

    Pendulum[0].init(theta0, omega0 , m0, R0, l0, x0);
    for(ii=1; ii<N; ii++) Pendulum[ii].init(0, 0 , m0, R0, l0, 2*R0*ii);

    StartAnimation();

    for(t=0 , tdraw = 0; t < tmax; t+= dt , tdraw+= dt){
        if(tdraw > tdomain){
            StartDomain();
            for(ii=0 ;ii<N; ii++) Pendulum[ii].Draw();
            FinishDomain();
            tdraw = 0;
        }
        // std::cout << Pendulum[1].get_x() << "\t" << Pendulum[1].get_y() << std::endl;
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_theta(dt, Zeta);

        Newton.compute_torsion(Pendulum);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_omega(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_theta(dt, Chi);

        Newton.compute_torsion(Pendulum);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_omega(dt, Lambda);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_theta(dt, Coeff2);

        Newton.compute_torsion(Pendulum);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_omega(dt, Lambda);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_theta(dt, Chi);

        Newton.compute_torsion(Pendulum);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_omega(dt, Coeff1);
        for(ii = 0; ii < N; ii++) Pendulum[ii].move_theta(dt, Zeta);
    }

    return 0;
}