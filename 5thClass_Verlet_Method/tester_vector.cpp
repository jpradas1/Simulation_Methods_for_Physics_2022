#include <iostream>
#include <cmath>

#include "vector.h"

void modify(vector3D & a){
    a += a; 
}

int main(){
    vector3D a,b,c;
    double y;

    a.load(0,1,2);
    b.load(8,3,3);
    c.load(5,2,1);

    y = a.y();

    a.show();
    std::cout << y << "\t" << &y << std::endl;
    modify(a);
    a.show();


    return 0;
}