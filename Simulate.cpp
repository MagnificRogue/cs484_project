#ifndef SIMULATE_CPP
#define SIMULATE_CPP

#include "Candidate.h"
#include <iostream>
#include "armadillo"
#include <cmath>

//using namespace std;
using namespace arma;

const mat A = { {1, 1, 0, 0, 0},
                {0, 1, 1, 0, 0}, 
                {0, 0, 1, 1, 0},
                {0, 0, 0, 1, 1}};

const mat D = { {1, -1, 0, 0, 0},
                {0, 1, -1, 0, 0}, 
                {0, 0, 1, -1, 0},
                {0, 0, 0, 1, -1}};


double simulate(Candidate candidate) {
  mat M(5,5, fill::eye);
  mat L = M;
  mat J = M/3.0;

  
  mat H = L*A.t()*inv(D*(inv(M))*D.t())*A*L;
  mat N = inv(M)*D.t()*inv(D*(inv(M))*D.t())*A*L;

  mat Cn = M*10.0;
  mat Ct = M*0.1;
  mat Dt = Cn*J;

  vec e = {5,fill::ones};
  mat E = zeros<mat>(10,2);
  E.submat(0,0,4,0) = e;
  E.submat(5,1,9,1) = e;

  mat Df = zeros<mat>(10,10);
  Df.submat(0,0,4,4) = Ct*M;
  Df.submat(5,5,9,9) = Cn*M;

  mat AA = zeros<mat>(14,14);
  mat BB = zeros<mat>(14,14);
  mat CC = zeros<mat>(14,4);
  mat DD = zeros<mat>(14,14);

  CC.submat(7,0,11,3) = D.t();


  double dt = 0.001;
  double time = 0.0;
  double TotalTime = 10.0;
  

// target: theta01 to theta05, x and y, theta01_dot to theta05_dot, x_dot and y_dot
  vec y = {0.0, 0.0, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  
  while (time < TotalTime) {
    
    vec SinY = {sin(y(0)), sin(y(1)), sin(y(2)), sin(y(3)), sin(y(4))};
    mat So = diagmat(SinY);

    vec CosY = {cos(y(0)), cos(y(1)), cos(y(2)), cos(y(3)), cos(y(4))};
    mat Co = diagmat(CosY);

    mat Jie = J + So*H*So + Co*H*Co;
    mat C = So*H*Co - Co*H*So;

    mat Lie = zeros<mat>(10,5);
    Lie.submat(0,0,4,4) = So*N.t();
    Lie.submat(5,0,9,4) = -Co*N.t(); 

    mat Omi = zeros<mat>(10,10);
    Omi.submat(0,0,4,4) = Co;
    Omi.submat(0,5,4,9) = -So;
    Omi.submat(5,0,9,4) = So;
    Omi.submat(5,5,9,9) = Co;

    mat RSQ1 = zeros<mat>(7,7);
    RSQ1.submat(0,0,4,4) = Dt;
    mat LieE = zeros<mat>(10,7);
    LieE.submat(0,0,9,4) = Lie;
    LieE.submat(0,5,9,6) = E;
    mat RSQ = RSQ1 + LieE.t()*Omi*Df*Omi.t()*LieE; 

    mat R = RSQ.submat(0,0,4,4);
    mat S = RSQ.submat(0,5,4,6);
    mat Q = RSQ.submat(5,5,6,6);

    AA.submat(0,0,6,6) = eye<mat>(7,7);
    AA.submat(7,7,11,11) = Jie;
    AA.submat(12,12,13,13) = 5*eye<mat>(2,2);

    BB.submat(0,7,6,13) = eye<mat>(7,7);
    BB.submat(7,7,11,11) = -R;
    BB.submat(7,12,11,13) = -S;
    BB.submat(12,7,13,11) = -S.t();
    BB.submat(12,12,13,13) = -Q;

    DD.submat(7,7,11,11) = C;
    
    vec U = {candidate.torque[0]*sin(time), candidate.torque[1]*sin(time), candidate.torque[2]*sin(time), candidate.torque[3]*sin(time)};


    vec squareY = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,y(7)*y(7), y(8)*y(8),y(9)*y(9),y(10)*y(10), y(11)*y(11), 0.0, 0.0};
    vec dy = inv(AA)*(BB*y+CC*U+DD*(squareY));

    y += dy*dt;
   
    time += dt;


    std::cout << "x=" << y(5) << std::endl;
    std::cout << "y=" << y(6) << std::endl;

  }
  
  return std::sqrt( std::pow(2.5 - y(5), 2) + std::pow(0-y(6), 2) );
}

#endif
