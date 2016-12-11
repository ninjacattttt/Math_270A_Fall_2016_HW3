#ifndef _SYMMETRIC_TRIDIAGONAL_
#define _SYMMETRIC_TRIDIAGONAL_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "ImplicitQRSVD.h"

namespace JIXIE {

/**
    Symmetric tridiagonal matrix class.
*/
template <class T>
class SymmetricTridiagonal {

typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVect;
typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMat;
TVect alpha,beta;
int n;

//LDLT
TVect d,mu;

//QR
std::vector<GivensRotation<T> > Q;//the ith entry in the array introduces a zero in the i+1,i entry
TVect r1,r2,r3;

public:
  SymmetricTridiagonal(const int n_input):n(n_input){alpha.resize(n_input);beta.resize(n_input-1);
    alpha=TVect::Zero(n);
    beta=TVect::Zero(n-1);
  }
  SymmetricTridiagonal(const TMat& A):n(A.rows()){
    assert(A.rows()==A.cols());
    alpha.resize(A.rows());
    beta.resize(A.rows()-1);
    *this=A;
  }

  void SetToZero(){
    alpha(0)=(T)0;
    for(int i=1;i<n;i++){
      alpha(i)=(T)0;
      beta(i-1)=(T)0;}
  }

  void Multiply(const TVect& x,TVect& result)const {
    if(result.size()!=n || n<2 || x.size()!=n){
      std::cout<<"SymmetricTridiagonal of inappropriate size" << std::endl;
      exit(1);}

    result(0)=alpha(0)*x(0)+beta(0)*x(1);
    for(int i=1;i<n-1;i++) result(i)=beta(i-1)*x(i-1)+alpha(i)*x(i)+beta(i)*x(i+1);
    result(n-1)=beta(n-2)*x(n-2)+alpha(n-1)*x(n-1);
  }

  T& operator()(const int i,const int j){
    if(abs(i-j)>1 || i<0 || i>n-1 || j<0 || j>n-1){
      exit(1);
      return alpha(0);}
    else if(i==j)
      return alpha(i);
    else
      return beta(std::min(i,j));
  }

  void operator=(const TMat& A){
    assert(A.rows()==n && A.cols()==n);
    alpha(0)=A(0,0);
    for(int i=1;i<n;i++){
      alpha(i)=A(i,i);
      beta(i-1)=A(i-1,i);}
  }

  void Set(TMat& A)const{
    A.resize(n,n);
    A=TMat::Zero(n,n);
    A(0,0)=alpha(0);
    A(0,1)=beta(0);
    A(1,0)=beta(0);
    for(int i=1;i<n-1;i++){
      A(i,i)=alpha(i);
      A(i,i+1)=beta(i);
      A(i+1,i)=beta(i);
    }
    A(n-1,n-1)=alpha(n-1);
  }

  void LDLT(){
    d.resize(n);
    mu.resize(n-1);

    d(0)=alpha(0);
    for(int i=1;i<n;i++){
      mu(i-1)=beta(i-1)/d(i-1);
      d(i)=alpha(i)-mu(i-1)*mu(i-1)*d(i-1);
    }
  }

  void Set_L(TMat& L){
    if(d.size()!=n || mu.size()!=n-1) return;//LDLT is stale
    L.resize(n,n);
    L=TMat::Zero(n,n);
    L(0,0)=(T)1;
    for(int i=1;i<n;i++){
      L(i,i-1)=mu(i-1);
      L(i,i)=(T)1;
    }
  }

  void Set_D(TMat& D)const{
    if(d.size()!=n || mu.size()!=n-1) return;//LDLT is stale
    D.resize(n,n);
    D=TMat::Zero(n,n);
    for(int i=0;i<n;i++)
      D(i,i)=d(i);
  }

  void Set_R(TMat& R)const{
    R.resize(n,n);
    R=TMat::Zero(n,n);
    for(int i=0;i<n-2;i++){
      R(i,i)=r1(i);R(i,i+1)=r2(i);R(i,i+2)=r3(i);
    }
    R(n-2,n-2)=r1(n-2);R(n-2,n-1)=r2(n-2);
    R(n-1,n-1)=r1(n-1);
  }

  void Set_Q(TMat& Qm)const{
    Qm.resize(n,n);
    Qm=TMat::Identity(n,n);
    for(int i=n-2;i>=0;i--)
      Q[i].transposeRowRotation(Qm);
  }

  void QRowRotation(TVect& x){
    for(int i=n-2;i>=0;i--)
      Q[i].transposeRowRotation(x);
  }

  void QTransposeRowRotation(TVect& x){
    for(int i=0;i<n-1;i++)
      Q[i].rowRotation(x);
  }

  void QRSolve(TVect& x,const TVect& b){
    QR();
    TVect QTb;
    QTb.resize(n);
    QTb=b;
    QTransposeRowRotation(QTb);
    x(n-1)=QTb(n-1)/r1(n-1);
    x(n-2)=(QTb(n-2)-r2(n-2)*x(n-1))/r1(n-2);
    for(int i=n-3;i>=0;i--)
      x(i)=(QTb(i)-r2(i)*x(i+1)-r3(i)*x(i+2))/r1(i);
  }

  void QR(){
    typedef Eigen::Matrix<T,2,1> TV2;

    assert(n>2);
    Q.resize(n-1);
    r1.resize(n);
    r2.resize(n-1);
    r3.resize(n-2);

    T alpha_bar=alpha(0);
    T beta_bar=beta(0);

    for(int j=0;j<n-2;j++){
      //compute Givens that defines rjj,rj+1,rjj+2
      Q[j].SetIK(j,j+1);
      Q[j].compute(alpha_bar,beta(j));
      r1(j)=std::sqrt(alpha_bar*alpha_bar+beta(j)*beta(j));
      //compute rjj+1 and alpha bar
      TV2 temp(beta(j),alpha(j+1));
      Q[j].rowRotation(temp);
      r2(j)=temp(0);
      alpha_bar=temp(1);
      //compute rjj+2 and beta_bar
      temp(0)=(T)0;temp(1)=beta(j+1);
      Q[j].rowRotation(temp);
      r3(j)=temp(0);
      beta_bar=temp(1);
    }
    //the last givens only effects two columns of r
    //compute Givens that defines rjj,rj+1,rjj+2
    Q[n-2].SetIK(n-2,n-1);
    Q[n-2].compute(alpha_bar,beta(n-2));
    r1(n-2)=std::sqrt(alpha_bar*alpha_bar+beta(n-2)*beta(n-2));
    //compute rjj+1 and alpha bar
    TV2 temp(beta_bar,alpha(n-1));
    Q[n-2].rowRotation(temp);
    r2(n-2)=temp(0);
    r1(n-1)=temp(1);
  }
  
  void Givens(const T x, const T y, T& c, T& s){
  	T d=x*x+y*y,t;
  	c=1;
  	s=0;
  	if(d!=0){
  		t=std::sqrt(d);
  		c=x/t;
  		s=-y/t;
  	}
  }
  
  void SetColumn(TMat& A, const TVect& q, int k, int n){			//set column k, in the first n rows
	  for(int i=0;i<n;i++)
		  A(i,k)=q(i);
  }
  
  void SetVector(const TMat& A, TVect& x, int k, int n){			//set vector x to be column k of A
	  for(int i=0;i<n;i++)
		  x(i)=A(i,k);
  }
  
  void norm(const TVect& b, T& norm){						//norm of a vector with length n
	  norm=(T) 0;
	  for(int i=0;i<n;i++)
		  norm+=b(i)*b(i);
	  norm=sqrt(norm);
  }
  
  void BackSub(const int& k, const TMat& A, TVect& b, TVect& x){
	  x(k-1)=b(k-1)/A(k-1,k-1);
	  for(int i=k-2;i>=0;i--){
		  for(int j=0;j<i;j++)
			  b(j)-=x(i)*A(j,i);
		  x(i)=b(i)/A(i,i);
	  }
  }
  
  void exit_GMRES(const int& k, const TMat& A, TVect& b, TVect& x, TVect& lambda, const TMat& Q){
	  BackSub(k, A, b, lambda);
	  x=Q*lambda;
  }
  
  void GMRES(TVect& x, const TVect& b){
	  T tol=1e-15;
	  TVect lambda;
	  lambda.resize(n);
	  lambda.setZero();
	  
	  TMat A;
	  Set(A);
	  
	  TMat Q, H,G;
	  Q.resize(n,n);
	  Q=TMat::Zero(n,n);
	  H.resize(n+1,n);
	  H=TMat::Zero(n+1,n);
	  G.resize(n,n);
	  G=TMat::Identity(n,n);
	  
	  TVect q,vtemp,htemp;
	  T normtemp,c,s;
	  
	  int k=0;
	  q.resize(n);
	  vtemp.resize(n);
	  htemp.resize(n);
	  
	  norm(b,normtemp);
	  q=b/normtemp;						//q_1=b/|b|
	  SetColumn(Q,q,0,n);
	  
	  TVect t;
	  t.resize(n);
	  t.setZero();
	  t(0)=normtemp;
	  
	  H(0,0)=q.transpose()*A*q;
	  vtemp=A*q-H(0,0)*q;
	  norm(vtemp,normtemp);
	  H(1,0)=normtemp;
	  q=vtemp/normtemp;
	  SetColumn(Q,q,1,n);
	  Givens(H(0,0),H(1,0),c,s);
	  G(0,0)=c;							//G_1 transpose
	  G(0,1)=-s;
	  G(1,0)=s;
	  G(1,1)=c;
	  t=G*t;
	  H=G*H;
	  
	  if(std::abs(t(1))<tol){
		  exit_GMRES(k,H,t,x,lambda,Q);
		  return;
	  }
	  
	  for(k=1;k<n;k++){
		  htemp=Q.transpose()*A*q;
		  vtemp=A*q-Q*htemp;
		  norm(vtemp,normtemp);
		  SetColumn(H,htemp,k,k);
		  H(k+1,k)=normtemp;
		  q=vtemp/normtemp;
		  SetColumn(Q,q,k+1,n);
		  Givens(H(k,k),H(k+1,k),c,s);
		  G=TMat::Identity(n,n);
		  G(k,k)=c;							//G_1 transpose
		  G(k,k+1)=-s;
		  G(k+1,k)=s;
		  G(k+1,k+1)=c;
		  t=G*t;
		  H=G*H;
		  
		  if(std::abs(t(k+1))<tol){
			  exit_GMRES(k,H,t,x,lambda,Q);
			  return;
		  }
		
		exit_GMRES(k,H,t,x,lambda,Q);
		
	  }
	  
	  
	
	
  }
  
};
}
#endif
