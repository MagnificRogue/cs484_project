#include <iostream>
#include "armadillo"

using namespace std;
using namespace arma;

int main(){

  int n = 5;
  contant double pi = 3.141592653;
  vec Candidate = 30*ones<vec>(n-1);
  Candidate(0) = 30;
  Candidate(1) = 30;
  Candidate(2) = 30;
  Candidate(3) = 30;
  vec Phase = zeros<vec>(n-1);
  for (int i=0; i<n-1; i++) {
    Phase(i) = pi*i/(n-2);
  }

  vec MASS = 1.0*ones<vec>(n);
  vec Length = 1.0*ones<vec>(n);

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

  double dt = 0.001;
  double time = 0.0;
  int step = 0;
  double TotalTime = 20.0;
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
      U(i) = Candidate(i)*sin(time+Phase(i));
    }

    vec squareY = zeros<vec>(2*(n+2));
    for (int i=n+2; i<2*n+2; i++){
        squareY(i) = y(i)*y(i); 
    }

    vec dy = inv(AA)*(BB*y+CC*U+DD*(squareY));

    y += dy*dt;
   
    if(step%PrintStep == 0){
       vec xout = zeros<vec>(n+1);
       vec yout = zeros<vec>(n+1);
       for (int i=1; i<n+1; i++){    
          xout(i) = xout(i-1) + L(i-1,i-1)*cos(y(i-1));
          yout(i) = yout(i-1) + L(i-1,i-1)*sin(y(i-1));
       }

       vec xcenter = zeros<vec>(n);
       vec ycenter = zeros<vec>(n);      
       for (int i=0; i<n; i++){
          xcenter(i) = 0.5*(xout(i)+xout(i+1));
          ycenter(i) = 0.5*(yout(i)+yout(i+1));
       }
       
       mat XCENTER = (1.0/m)*e.t()*M*xcenter;
       mat YCENTER = (1.0/m)*e.t()*M*ycenter;

       for (int i=0; i<n+1; i++){
          xout(i) += y(n)-XCENTER(0,0);
          yout(i) += y(n+1) - YCENTER(0,0);
       }

       ofstream myfile;
       myfile.open ("coord.txt", fstream::app);
       myfile << fixed << time << "\n";
       for (int i=0; i<n+1; i++){
       myfile << fixed << xout(i) << "\t" << fixed << yout(i) << "\n";
       }
       myfile << "\n";
       myfile.close(); 

    }

    time += dt;
    step ++;
  }

  
  
  return 0;
}
