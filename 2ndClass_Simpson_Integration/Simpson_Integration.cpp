#include <iostream>
#include <cmath>

const double ErrMax = 1e-15;

double f(double x){
    return cos(x);
}

double Simpson(double a, double b, int n){
    double x, h, sum = 0;
    n *= 2; h = (b-a)/n;

    for(int ii = 0 ; ii<=n ; ii++){
        x = a + ii*h;
        if(ii==0 || ii==n)
            sum += f(x);
        else if(ii%2==0)
            sum += 2*f(x);
        else
            sum += 4*f(x);
    }

    return sum*h/3;
}

int main(){
    double a =  0; double b = M_PI/2;
    int n = 50;

    std::cout << "The integral by Simpson Methon is "<<Simpson(a,b,n)<<std::endl;

    return 0;
}