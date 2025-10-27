#ifndef function_hpp
#define function_hpp

#include <cmath>

class functions{
public:
 functions(){}

 double func(double a,double r){
  return 1.0/(1.0+exp(-r/a))-a;
 }

 double gradfun(double a,double r){
  return -r*exp(-r/a)/(pow(1.0+std::exp(-r / a),2)*a*a)-1;
 }

 void cusp(double& a,double r){
    double deltaar;
    //initguess
    a=0.4;
  do{
  deltaar=func(a,r)/gradfun(a,r);
  a=a-deltaar;
  }while (abs(func(a,r))>1.0e-14);
 }

};

#endif