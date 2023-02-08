#include <cmath>
#include "./../Requirements/Random64.h"

const int Lx = 1024;
const double p = 0.5;

const int Q = 2;

class LatticeGas{
public:
    LatticeGas();
    void Delete();
    void Init(int N, double mu, double sigma, Crandom & rand64);
    void Show();
    void Display_rho();
    void Collision(Crandom & rand64);
    void Advection();
    double rho(int ix){return n[ix][0]+n[ix][1];};
    double Variance();
private:
    int n[Lx][Q], nnew[Lx][Q];
    int V[Q]; // V[i] i= 0 (right) i= 1 (left)
};

LatticeGas::LatticeGas(){
    V[0] = 1; V[1] = -1;
}

void LatticeGas::Delete(){
    for(int ix = 0; ix < Lx; ix++)
        for(int i=0; i<Q; i++)
            n[ix][i] = nnew[ix][i] = 0;
}

void LatticeGas::Init(int N, double mu, double sigma, Crandom & rand64){
    int ix, i;
    while (N>0){
        ix = (int) rand64.gauss(mu, sigma); // choose a random cell
        if(ix<0) ix = 0; if(ix>Lx-1) ix = Lx-1; // fit boundaries
        i = (int) Q * rand64.r(); // choose a random direction
        if(n[ix][i]==0)
            n[ix][i] = 1; N--; // fill that chosen cell
    }
    
}

void LatticeGas::Show(){
    for(int i = 0; i<Q; i++){
        for(int ix = 0; ix<Lx; ix++)
         std::cout << n[ix][i];
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void LatticeGas::Display_rho(){
    for(int ix = 0; ix < Lx; ix++)
        std::cout << ix << "\t" << rho(ix) << std::endl;
}

void LatticeGas::Collision(Crandom & rand64){
    for(int ix=0; ix<Lx; ix++){
        if(rand64.r() > p){
            // switch the content
            nnew[ix][0] = n[ix][1];
            nnew[ix][1] = n[ix][0];
        }
        else{ // no switch
            nnew[ix][0] = n[ix][0];
            nnew[ix][1] = n[ix][1];
        }
    }
}

void LatticeGas::Advection(){
    for(int ix=0; ix<Lx; ix++)
        for(int i=0; i<Q; i++)
            n[(ix+V[i]+Lx)%Lx][i] = nnew [ix][i];
}

double LatticeGas::Variance(){
    int ix; double N, Xprom, Sigma2;
    // compute N
    for(N=0, ix=0; ix<Lx; ix++)
        N+=rho(ix);
    // compute Xprom
    for(Xprom=0, ix=0; ix<Lx; ix++)
        Xprom += ix * rho(ix);
    Xprom /= N;
    // compute Sigma2
    for(Sigma2=0, ix=0; ix<Lx; ix++)
        Sigma2 += pow(ix-Xprom, 2.0) * rho(ix);
    Sigma2 /= (N-1);

    return Sigma2;
}

int main(){
    LatticeGas Diffusion;
    Crandom rand64(1);
    int N = 400; double mu = Lx/2, sigma = Lx/8;
    int t, tmax = 400;

    Diffusion.Delete();
    Diffusion.Init(N, mu, sigma, rand64);
    for(t=0; t<tmax;t++){
        std::cout << t << "\t" << Diffusion.Variance() << std::endl;
        Diffusion.Collision(rand64);
        Diffusion.Advection();
    }
    // Diffusion.Display_rho();

    return 0;
}