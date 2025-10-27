#include <iostream>
#include <random>
#include <algorithm>

#include "param.hpp"
#include "func.hpp"
#include "random.hpp"
#include "metro_el.hpp"

int main(int argc,char* agrv[])
{
functions fun;
metropolis mtsp;

double beta, a;
double E[Parameters::step],x[3][Parameters::num_walker][Parameters::num_ele];
double phi_grad[Parameters::step],E_grad[Parameters::step];
double energy;

beta=0.4;

energy=0;

fun.cusp(a, Parameters::R);

 for (int i=0;i<Parameters::beta_step;i++){


fill_random_uniform_3d(&x[0][0][0],3,Parameters::num_walker,
    Parameters::num_ele,-0.5,0.5);

mtsp.sampling(x,beta,a,Parameters::R,E,phi_grad,E_grad);
mtsp.newbeta(beta,E,phi_grad,E_grad);
}
for(int i=Parameters::eq_step-1;i<Parameters::step;i++){
    energy+=E[i];
}
double N = static_cast<double>(Parameters::step-Parameters::eq_step);
energy=energy/N;

std::cout<<"Total energy of VMC of H2 is "<<energy<<" Hartree"<<std::endl;

return 0;
}