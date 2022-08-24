#include <iostream>
#include <cmath>

const double ErrMax = 1e-7;

double f(double alpha, double x, double t){
    return cos(alpha*t-x*sin(t));
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

double Zeros_by_bisection(double alpha, double a, double b){
    double m, fa, fm;

    fa = Bessel(alpha, a);

    while(b-a > ErrMax){
        m = (b+a)/2, fm = Bessel(alpha, m);
        if(fa*fm>0)
            {a = m; fa = fm;}
        else
            b = m;

    }

    return (a+b)/2;
}

int main(){
    double alpha = 0;
    double a=2, b=4;

    // for(x=0; x<= 10; x+=0.1)
    //     std::cout << x << "\t" << Bessel(alpha, x) << std::endl;


    std::cout << "The Bessel Fucntion, with alpha = "<<alpha<<", vanish at "<<Zeros_by_bisection(alpha, a, b)<<std::endl;

    return 0;
}