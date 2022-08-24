#include <iostream>
#include <cmath>

const double g=9.8; // m/s²


class Body{
    public:
        void init(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
        void compute_Forces();
        void movement(double dt);

        double get_x(){return x;};
        double get_y(){return y;};
    private:
        double x, y, Vx, Vy, Fx, Fy, m, R;
};

void Body::init(double x0, double y0, double Vx0, double Vy0, double m0, double R0){
    x= x0; y= y0; Vx= Vx0; Vy= Vy0; m= m0; R= R0;
}

void Body::compute_Forces(){
    Fx= 0; Fy= -m*g;
}

void Body::movement(double dt){
    x+= Vx * dt;    y+= Vy * dt;
    Vx+= Fx/m*dt;   Vy+= Fy/m*dt;
}

int main(){
    double t, dt=0.1;
    Body Ball;

    // x0, y0, Vx0 ,Vy0, m0, R0

    Ball.init(0, 0 , 20, 15, 0.453, 0.15);

    for(t=0; t<3.5; t+= dt){
        std::cout << Ball.get_x() << "\t" << Ball.get_y() << std::endl;
        Ball.compute_Forces();
        Ball.movement(dt);
    }

    return 0;
}