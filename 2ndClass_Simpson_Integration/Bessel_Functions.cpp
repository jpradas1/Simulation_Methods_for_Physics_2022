#include <iostream>
#include <cmath>

double f(double alpha, double x, double t){
    return cos(alpha*x-x*sin(t));
}

double Simpson(double alpha, double x, double a, double b, int n){
    double t, h, sum = 0;
    n *= 2; h = (b-a)/n;

    for(int ii = 0 ; ii<=n ; ii++){
        t = a + ii*h;
        if(ii==0 || ii==n)
            sum += f(alpha, x, t);
        else if(ii%2==0)
            sum += 2*f(alpha, x, t);
        else
            sum += 4*f(alpha, x, t);
    }

    return sum*h/3;
}

double Bessel(double alpha, double x){
    double a = 0, b = M_PI; int n = 50;

    return 1.0/M_PI*Simpson(alpha, x, a, b, n);
}

int main(){
    double alpha = 0, x;

    for(x=0; x<= 10; x+=0.1)
        std::cout << x << "\t" << Bessel(alpha, x) << std::endl;


    // std::cout << "The integral by Simpson Methon is "<<Simpson(alpha, x, a,b,n)<<std::endl;

    return 0;
}