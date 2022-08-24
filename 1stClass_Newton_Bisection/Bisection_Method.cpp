#include <iostream>
#include <cmath>

const double ErrMax = 1e-15;

double f(double x){
    return sin(x)/x;
}

double Zeros_by_bisection(double a, double b){
    double m, fa, fm;

    fa = f(a);

    while(b-a > ErrMax){
        m = (b+a)/2, fm = f(m);
        if(fa*fm>0)
            {a = m; fa = fm;}
        else
            b = m;

    }

    return (a+b)/2;
}

int main(){
    double a = 2; double b = 4;

    // for(x=0.1; x<=10; x+=0.1){
    //     std::cout<<x<<"\t"<<f(x)<<std::endl;
    // }

    std::cout << "The function vanish at "<<Zeros_by_bisection(a,b)<<std::endl;

    return 0;
}