#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

#include "param.hpp"
#include "func.hpp"

int main(int argc, char* argv[]){

int walker_test, life[Parameters::max_walker], 
 test[Parameters::max_walker]={0};
double potentail_avg, W, nlife_before, nlife_now, energy=0;
double x[Parameters::max_walker][Parameters::num_ele][3]={0},xtemp[Parameters::num_ele][3];
double potential_value,Etemp=0;

functions fun;


for(int i=0;i<Parameters::num_walker;i++){
    life[i]=1;
}
for(int i=Parameters::num_walker;i<Parameters::max_walker;i++){
    life[i]=-1;
}


    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<double> normal_dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);





for(int i_iter=0;i_iter<Parameters::niter;i_iter++){
     
    for(int i=0;i<Parameters::max_walker;i++){
        if(life[i]>0){
            for(int j=0;j<Parameters::num_ele;j++){
                for(int k=0;k<3;k++){
                    x[i][j][k]=x[i][j][k]+sqrt(Parameters::dt)*normal_dist(gen);
                }
            }
        }
    }
 
 nlife_now=fun.lifecount(life);
 potentail_avg=0.0;
 for(int i=0;i<Parameters::max_walker;i++){
     for(int j=0;j<Parameters::num_ele;j++){
                for(int k=0;k<3;k++){
                    xtemp[j][k]=x[i][j][k];
                }
            }
    if(life[i]>0){
        potentail_avg+=fun.potential(xtemp,Parameters::R);
    }
 }
potentail_avg=potentail_avg/nlife_now;
nlife_before=nlife_now;

int count=0;

for (int i=0;i<Parameters::max_walker;i++){
   if (life[i]<0) continue;
        bool found = false;
        for (int j = 0; j < count+1; ++j) {
            if (test[j] == i) {
                found = true;
                break;
            }
        }
        if (found) continue;
    for(int j=0;j<Parameters::num_ele;j++){
        for(int k=0;k<3;k++){
        xtemp[j][k]=x[i][j][k];
        }
    }
    W = exp(-(fun.potential(xtemp, Parameters::R) - potentail_avg) 
       * Parameters::dt); 
     
    double rho=W+uniform_dist(gen);
    walker_test=std::min(int(rho),3);
    if(walker_test==0){
        life[i]=-1;
    }
     else if(walker_test==1){
        life[i]=1;
     }
    else if(walker_test==2){
        
int life_index = -1;
for (int k = 0; k < Parameters::max_walker; ++k) {
    if (life[k] < 0) {
        life_index = k;
        break;  
    }
}


    for (int j = 0; j < Parameters::num_ele; ++j) {
        for (int k = 0; k< 3; ++k) {
            x[life_index][j][k] = x[i][j][k];
        }
    }

    life[life_index] = 1;
    test[count] = life_index;
    count = count + 1;


}
    else if(walker_test==3){
        int life_index = -1;
  for (int iii=0;iii<2;++iii){

        
int life_index ;
for (int k = 0; k < Parameters::max_walker; ++k) {
    if (life[k] < 0) {
        life_index = k;
        break;  
    }
}


    for (int j = 0; j < Parameters::num_ele; ++j) {
        for (int k = 0; k< 3; ++k) {
            x[life_index][j][k] = x[i][j][k];
        }
    }

    life[life_index] = 1;
    test[count] = life_index;
    count = count + 1;
}
  
  }

    }
  
    
   
if(i_iter>Parameters::eq_iter){
    potential_value=0.0;
    for (int i=0;i<Parameters::max_walker;i++){
        if(life[i]<0) continue;
        for(int j=0;j<Parameters::num_ele;j++){
        for(int k=0;k<3;k++){
        xtemp[j][k]=x[i][j][k];
        }
        potential_value+=fun.potential(xtemp,Parameters::R);
    }
}
potential_value=potential_value/nlife_now;
Etemp=Etemp+potential_value;
}

}
double v,t;
v=Etemp/(Parameters::niter-Parameters::eq_iter);
energy=potentail_avg;
t=energy-v;

std::cout<<"Total energy of DMC of H2 is "<<energy<<" Hartree"<<std::endl;
return 0;

}