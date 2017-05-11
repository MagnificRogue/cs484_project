#ifndef SIMULATE_CPP
#define SIMULATE_CPP

#include "Candidate.h"
#include <iostream>
#include "armadillo"
#include <cmath>

using namespace arma;

double simulate(Candidate candidate) {
  int n = candidate.numberOfLinks;
  vec Phase = zeros<vec>(n-1);

  for (int i=0; i<n-1; i++) {
    Phase(i) = pi*i/(n-2);
  }

  vec MASS = 1.0*ones<vec>(n);
  vec Length = 1.0*ones<vec>(n);

  for(int i=0; i < n; i++){
    MASS(i) = candidate.weight[i]; 
    Length(i) = candidate.length[i];
  }

  mat A = zeros<mat>(n-1, n);
  for (int i=0; i<(n-1); i++){
    A(i,i) = 1.0;
    A(i,i+1) = 1.0;
  }

  mat D = zeros<mat>(n-1, n);
  for (int i=0; i<(n-1); i++){
    D(i,i) = 1.0;
    D(i,i+1) = -1.0;
  }

  mat M = diagmat(MASS);
  double m = trace(M);
  mat L = diagmat(Length);
  mat J = M*L*L/3.0;
  
  mat H = L*A.t()*inv(D*(inv(M))*D.t())*A*L;
  mat N = inv(M)*D.t()*inv(D*(inv(M))*D.t())*A*L;

  double f_n = 10.0;
  double f_t = 0.1;

  mat Cn = M*f_n;
  mat Ct = M*f_t;
  mat Dt = Cn*J;

  vec e = ones<vec>(n);
  mat E = zeros<mat>(2*n,2);
  E.submat(0,0,n-1,0) = e;
  E.submat(n,1,2*n-1,1) = e;

  mat Df = zeros<mat>(2*n,2*n);
  Df.submat(0,0,n-1,n-1) = Ct*M;
  Df.submat(n,n,2*n-1,2*n-1) = Cn*M;

  mat AA = zeros<mat>(2*(n+2),2*(n+2));
  mat BB = zeros<mat>(2*(n+2),2*(n+2));
  mat CC = zeros<mat>(2*(n+2),n-1);
  mat DD = zeros<mat>(2*(n+2),2*(n+2)); 

  CC.submat(n+2,0,2*n+1,n-2) = D.t();

  double dt = 0.0001;
  double time = 0.0;
  int step = 0;
  double TotalTime = 5.0;
  int PrintStep = 100; 


  vec y = zeros<vec>(2*(n+2));
  
  while (time < TotalTime) {
    
    vec SinY = zeros<vec>(n);
    for (int i=0; i<n; i++){
       SinY(i) = sin(y(i));
    }
    mat So = diagmat(SinY);


    vec CosY = zeros<vec>(n);
    for (int i=0; i<n; i++){
       CosY(i) = cos(y(i));
    }
    mat Co = diagmat(CosY);


    mat Jie = J + So*H*So + Co*H*Co;
    mat C = So*H*Co - Co*H*So;

    mat Lie = zeros<mat>(2*n,n);
    Lie.submat(0,0,n-1,n-1) = So*N.t();
    Lie.submat(n,0,2*n-1,n-1) = -Co*N.t(); 

    mat Omi = zeros<mat>(2*n,2*n);
    Omi.submat(0,0,n-1,n-1) = Co;
    Omi.submat(0,n,n-1,2*n-1) = -So;
    Omi.submat(n,0,2*n-1,n-1) = So;
    Omi.submat(n,n,2*n-1,2*n-1) = Co;

    mat RSQ1 = zeros<mat>(n+2,n+2);
    RSQ1.submat(0,0,n-1,n-1) = Dt;
    mat LieE = zeros<mat>(2*n,n+2);
    LieE.submat(0,0,2*n-1,n-1) = Lie;
    LieE.submat(0,n,2*n-1,n+1) = E;
    mat RSQ = RSQ1 + LieE.t()*Omi*Df*Omi.t()*LieE; 

    mat R = RSQ.submat(0,0,n-1,n-1);
    mat S = RSQ.submat(0,n,n-1,n+1);
    mat Q = RSQ.submat(n,n,n+1,n+1);

    AA.submat(0,0,n+1,n+1) = eye<mat>(n+2,n+2);
    AA.submat(n+2,n+2,2*n+1,2*n+1) = Jie;
    AA.submat(2*n+2,2*n+2,2*n+3,2*n+3) = m*eye<mat>(2,2);

    BB.submat(0,n+2,n+1,2*n+3) = eye<mat>(n+2,n+2);
    BB.submat(n+2,n+2,2*n+1,2*n+1) = -R;
    BB.submat(n+2,2*n+2,2*n+1,2*n+3) = -S;
    BB.submat(2*n+2,n+2,2*n+3,2*n+1) = -S.t();
    BB.submat(2*n+2,2*n+2,2*n+3,2*n+3) = -Q;

    DD.submat(n+2,n+2,2*n+1,2*n+1) = C;
  
    vec U = zeros<vec>(n-1);
    for (int i=0; i<n-1; i++){
      U(i) = candidate.torque[i]*sin(time+Phase(i));
    }

    vec squareY = zeros<vec>(2*(n+2));
    for (int i=n+2; i<2*n+2; i++){
        squareY(i) = y(i)*y(i); 
    }

    vec dy = inv(AA)*(BB*y+CC*U+DD*(squareY));

    y += dy*dt;
    time += dt;
    step ++;
  }

  
  double outX = y(n);
  double outY = y(n+1);
 
  double answer = sqrt(pow(outX,2) + pow(outY,2));
  
  if(isnan(answer))
    answer = 0;

  return answer;

}
#endif
