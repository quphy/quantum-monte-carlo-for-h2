#ifndef metropolis_energy_local
#define metropolis_energy_local

#include <iostream>
#include <cmath>

#include "random.hpp"
#include "param.hpp"

class metropolis{
public:
metropolis(){}


template <int x1, int x2, int x3,
 int e1, int gb, int ge>
void sampling(double (&x)[x1][x2][x3],double& beta, const double& a ,
const double& R,double (&E)[e1], double(&grad_beta)[gb], double (&grad_E)[ge])
{
 double x_new[x1][x2][x3],rnd_uni[x1][x2][x3],rnd_uni2[x2],delta, transition_prob[x2],
        tran1[x2],tran2[x2];
 int count,mark[x2],mark_3d[x1][x2][x3];
 delta=1.0;
 count=0;
 double N = static_cast<double>(x2);
 for (int i=0;i<e1;i++){
  fill_random_uniform_3d(&rnd_uni[0][0][0],x1,x2,x3,-0.5,0.5);
  fill_random_uniform_1d(rnd_uni2,rnd_uni2+x2,0.0,1.0);
  for (int j = 0; j < x1; ++j){
    for (int k = 0; k < x2; ++k){
        for (int l = 0; l < x3; ++l){
            x_new[j][k][l] = x[j][k][l] + rnd_uni[j][k][l];
        }
    }
  }
 tran_prob(x,a,R,beta,tran1);
 tran_prob(x_new,a,R,beta,tran2);

 for (int j=0;j<x2;++j){
  transition_prob[j]=( (tran2[j]*tran2[j]) /(tran1[j]*tran1[j]) );
   if(rnd_uni2[j]<transition_prob[j]){
    mark[j]=1;
    count=count+1;
   }
   else{
    mark[j]=0;
   }
}
 
for (int j=0;j<x2;++j){
    if(mark[j]==1){
    for (int k=0;k<x1;++k){
        for ( int l=0;l<x3;++l){
                x[k][j][l]=x_new[k][j][l];
            }
        }
    }
}

if(i>Parameters::eq_step){
E_local(x,a,R,beta,E[i],grad_beta[i],grad_E[i]);
}

if((i+1)%10000==0){
    
    delta=delta*count/(50.0*N);
    count=0;
}
}
}

template <int x1, int x2, int x3>
void tran_prob(double (&x)[x1][x2][x3], const double& a, const double& R, double& beta,
               double (&prob)[x2])
{
double r1[x1][x2],r1L[x1][x2],r1R[x1][x2],r2[x1][x2],r2L[x1][x2],r2R[x1][x2],r12[x1][x2];
double phi1[x2],phi1_L[x2],phi1_R[x2],phi2[x2],phi2_L[x2],phi2_R[x2],r1L_abs[x2],
       r1R_abs[x2],r2L_abs[x2],r2R_abs[x2],r12_abs[x2];

value(x,a,R,r1,r1L,r1L_abs,r1R,r1R_abs,r2,r2L,r2L_abs,r2R,r2R_abs,r12,r12_abs,phi1,phi1_L,
      phi1_R,phi2,phi2_L,phi2_R);

for (int i=0;i<x2;i++){
prob[i]=phi1[i]*phi2[i]*exp(r12_abs[i]/(Parameters::alpha*(1.0+beta*r12_abs[i])));
}      
}

template <int x1, int x2, int x3>
void value(double (&x)[x1][x2][x3], const double& a , const double& R,
     double (&r1)[x1][x2], double (&r1L)[x1][x2], double (&r1L_abs)[x2],
     double (&r1R)[x1][x2], double (&r1R_abs)[x2], double (&r2)[x1][x2],
     double (&r2L)[x1][x2], double (&r2L_abs)[x2], double (&r2R)[x1][x2],
     double (&r2R_abs)[x2], double (&r12)[x1][x2], double (&r12_abs)[x2],
     double (&phi1)[x2], double (&phi1_L)[x2], double (&phi1_R)[x2],
     double (&phi2)[x2], double (&phi2_L)[x2], double (&phi2_R)[x2]){
   
 double nucle_cor[3][Parameters::num_nuc];
 
 for (int i=0;i<3;i++){
  for (int j=0;j<Parameters::num_nuc;j++){
    nucle_cor[i][j]=0;
  }
 }
 nucle_cor[0][0]=-0.5*R;
 nucle_cor[0][1]=0.5*R;
 
 for (int j=0;j<x2;j++){
    for (int i=0; i<x1;i++){
        r1[i][j]=x[i][j][0];
        r2[i][j]=x[i][j][1];
        r1L[i][j]=r1[i][j]-nucle_cor[i][0];
        r1R[i][j]=r1[i][j]-nucle_cor[i][1];
        r2L[i][j]=r2[i][j]-nucle_cor[i][0];
        r2R[i][j]=r2[i][j]-nucle_cor[i][1];
        r12[i][j]=r1[i][j]-r2[i][j];
    }
    r1L_abs[j]=std::sqrt(r1L[0][j]*r1L[0][j]+r1L[1][j]*r1L[1][j]+r1L[2][j]*r1L[2][j]);
    r1R_abs[j]=std::sqrt(r1R[0][j]*r1R[0][j]+r1R[1][j]*r1R[1][j]+r1R[2][j]*r1R[2][j]);
    r2L_abs[j]=std::sqrt(r2L[0][j]*r2L[0][j]+r2L[1][j]*r2L[1][j]+r2L[2][j]*r2L[2][j]);
    r2R_abs[j]=std::sqrt(r2R[0][j]*r2R[0][j]+r2R[1][j]*r2R[1][j]+r2R[2][j]*r2R[2][j]);
    r12_abs[j]=std::sqrt(r12[0][j]*r12[0][j]+r12[1][j]*r12[1][j]+r12[2][j]*r12[2][j]);
    phi1_L[j]=exp(-r1L_abs[j]/a);
    phi1_R[j]=exp(-r1R_abs[j]/a);
    phi2_L[j]=exp(-r2L_abs[j]/a);
    phi2_R[j]=exp(-r2R_abs[j]/a);
    phi1[j]=phi1_L[j]+phi1_R[j];
    phi2[j]=phi2_L[j]+phi2_R[j];
}
 }

template <int x1,int x2, int x3>
void E_local(double(&x)[x1][x2][x3] , const double& a, const double& R, double& beta,
             double& E_st, double& grad_betast,double& grad_Est){
double r1[x1][x2],r1L[x1][x2],r1R[x1][x2],r2[x1][x2],r2L[x1][x2],r2R[x1][x2],r12[x1][x2];
double phi1[x2],phi1_L[x2],phi1_R[x2],phi2[x2],phi2_L[x2],phi2_R[x2],r1L_abs[x2],
       r1R_abs[x2],r2L_abs[x2],r2R_abs[x2],r12_abs[x2];
double r1L_hat[x1][x2],r1R_hat[x1][x2],r2L_hat[x1][x2],r2R_hat[x1][x2],r12_hat[x1][x2];  
double E1[x2],E2[x2],E3[x2],E4[x2],E51[x2],E52[x2],E53[x2],E5[x2],E6[x2],E_tot[x2],grad_beta_array[x2];
double L1[x2],R1[x2],L2[x2],R2[x2]; 
double N = static_cast<double>(x2);
value(x,a,R,r1,r1L,r1L_abs,r1R,r1R_abs,r2,r2L,r2L_abs,r2R,r2R_abs,r12,r12_abs,phi1,phi1_L,
      phi1_R,phi2,phi2_L,phi2_R);
E_st=0;
grad_betast=0;
grad_Est=0;
 for(int i=0; i<x2; i++){
    for(int j=0; j<x1; j++){
        r1L_hat[j][i]=r1L[j][i]/r1L_abs[i];
        r1R_hat[j][i]=r1R[j][i]/r1R_abs[i];
        r2L_hat[j][i]=r2L[j][i]/r2L_abs[i];
        r2R_hat[j][i]=r2R[j][i]/r2R_abs[i];
        r12_hat[j][i]=r12[j][i]/r12_abs[i];
    }
  E1[i]=-1.0/(a*a);
  E2[i]=(phi1_L[i]/r1L_abs[i]+ phi1_R[i]/r1R_abs[i])/(a*phi1[i]);
  E3[i]=(phi2_L[i]/r2L_abs[i]+ phi2_R[i]/r2R_abs[i])/(a*phi2[i]);
  E4[i]=1.0/r12_abs[i]- 1.0/r1L_abs[i]- 1.0/r1R_abs[i]- 1.0/r2L_abs[i]- 1.0/r2R_abs[i];
  E53[i]=1.0/(2.0*a*(1.0+beta*r12_abs[i]) * (1.0+beta*r12_abs[i]) );
  L1[i]=0;
  R1[i]=0;
  L2[i]=0;
  R2[i]=0;
 }

 for(int i=0; i<x2; i++){
    for(int j=0;j<x1;j++){
        L1[i]+=r1L_hat[j][i]*r12_hat[j][i];
        R1[i]+=r1R_hat[j][i]*r12_hat[j][i];
        L2[i]+=r2L_hat[j][i]*r12_hat[j][i];
        R2[i]+=r2R_hat[j][i]*r12_hat[j][i];
    }
    E51[i]=(phi1_L[i]*L1[i]+phi1_R[i]*R1[i])/phi1[i];
    E52[i]=(phi2_L[i]*L2[i]+phi2_R[i]*R2[i])/phi2[i];
    E5[i]=(E51[i]-E52[i])*E53[i];
    E6[i]= - ((4.0*beta+1.0)*r12_abs[i]+4.0 )/(4.0*r12_abs[i]*std::pow((1+beta*r12_abs[i]),4));
    E_tot[i]=E1[i]+E2[i]+E3[i]+E4[i]+E5[i]+E6[i]+1.0/R;
    E_st+=E_tot[i];
    grad_beta_array[i]=-r12_abs[i]*r12_abs[i]/(Parameters::alpha*(1+beta*r12_abs[i])*(1+beta*r12_abs[i]));   
    grad_betast+=grad_beta_array[i];
    grad_Est+=E_tot[i]*grad_beta_array[i];
}
E_st=E_st/N;
grad_betast=grad_betast/N;
grad_Est=grad_Est/N;
}

template<int gb>
void newbeta(double& beta ,double (&E)[gb],double (&gradb)[gb],double (&gradE)[gb]){
double Esum,gradbsum,gradEsum,Eavg,gradbavg,gradEavg;

Esum=0;
gradbsum=0;
gradEsum=0;

double N = static_cast<double>(Parameters::step-Parameters::eq_step);
for (int i=Parameters::eq_step-1;i<gb;++i){
 Esum+=E[i];
 gradbsum+=gradb[i];
 gradEsum+=gradE[i];
}
Eavg=Esum/N;
gradbavg=gradbsum/N;
gradEavg=gradEsum/N;
beta=beta-Parameters::gamma*2*(gradEavg-Eavg*gradbavg);
}




};

#endif