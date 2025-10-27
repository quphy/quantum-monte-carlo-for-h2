#ifndef function_hpp
#define function_hpp

#include <cmath>

class functions{
public:
 functions(){}
 
 template<int x1>
 double lifecount(int (&life)[x1])
 {
    int nlife=0;
    for(int i=0;i<x1;i++){
        if(life[i]>0){
            nlife+=1;
        }
    }
    return nlife;
 }

 template<int x1,int x2>
 double potential(double (&x)[x1][x2],const double& R){
 double R1L,R1R,R2L,R2R,R12,po;
 R1L=x[0][1]*x[0][1]+x[0][0]*x[0][0]+(x[0][2]-R/2)*(x[0][2]-R/2);
 R1L=sqrt(R1L);
 R1R=x[0][1]*x[0][1]+x[0][0]*x[0][0]+(x[0][2]+R/2)*(x[0][2]+R/2);
 R1R=sqrt(R1R);
 R2L=x[1][1]*x[1][1]+x[1][0]*x[1][0]+(x[1][2]-R/2)*(x[1][2]-R/2);
 R2L=sqrt(R2L);
 R2R=x[1][1]*x[1][1]+x[1][0]*x[1][0]+(x[1][2]+R/2)*(x[1][2]+R/2);
 R2R=sqrt(R2R);
 R12=(x[1][2]-x[0][2])*(x[1][2]-x[0][2])+(x[1][0]-x[0][0])*(x[1][0]-x[0][0])
 +(x[1][1]-x[0][1])*(x[1][1]-x[0][1]);
 R12=sqrt(R12);
 po=-1.0/R1L-1.0/R1R-1.0/R2L-1.0/R2R+1.0/R12+1.0/R;
 return po;
 }



};

#endif
