#ifndef LEGENDRE1DELEMENT_HPP
#define LEGENDRE1DELEMENT_HPP 

#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "NumericalIntegration.hpp"
using namespace std;

//
// Some Supplementary Functions
// needed for calculations
//
template<typename Num>
Num MyCos(Num a){
  return cos(a);
};

template<typename Num>
Num MySqrt(Num a){
  return sqrt(a);
};

//
// 1-D Legendre polynomial shape functions
// uses a single integration rule so far
//

//
// Legendre polynomial function
//
template<typename Num, typename Iter>
Num LegendrePolyFunc(Num r, Iter I, Num Ln1, Num Ln2){
  Num LP   = Num(0.0);
  Num one  = Num(1.0);
  Num two  = Num(2.0);
  Num n    = Num(I);
  if(I==0)           LP = one;
  if(I==1)           LP = r;
  if((I!=0)and(I!=1)) LP = (one/n)*( (two*n - one)*r*Ln1 - (n - one)*Ln2 );
  return LP;
};

//
// Legendre polynomial function derivative
//
template<typename Num, typename Iter>
Num LegendrePolyFuncDer(Num r, Iter I, Num Ln1, Num dLn1, Num dLn2){
  Num Zero = Num(0.0);
  Num one  = Num(1.0);
  Num two  = Num(2.0);
  Num n    = Num(I);
  Num dLP  = Zero;
  if(I==0)           dLP = Zero;
  if(I==1)           dLP = one;
  if((I!=0)and(I!=1)) dLP = (one/n)*( (two*n - one)*(Ln1 + r*dLn1) - (n - one)*dLn2 );
  return dLP;
};

//
// Legendre Shape function
//
template<typename Num, typename Iter>
Num LegendreShapeFunc(Num r, Iter I, Num Ln0, Num Ln2){
  Num Ni   = Num(0.0);
  Num half = Num(0.5);
  Num one  = Num(1.0);
  Num two  = Num(2.0);
  Num four = Num(4.0);
  Num n    = Num(I);
  if(I==0)           Ni = half*(one + r);
  if(I==1)           Ni = half*(one - r);
  if((I!=0)and(I!=1)) Ni = (Ln0 - Ln2)/MySqrt<Num>(four*n-two);
  return Ni;
};

//
// Legendre Shape function derivative
//
template<typename Num, typename Iter>
Num LegendreShapeFuncDer(Iter I, Num dLn0, Num dLn2){
  Num dNi  = Num(0.0);
  Num half = Num(0.5);
  Num two  = Num(2.0);
  Num four = Num(4.0);
  Num n    = Num(I);
  if(I==0) dNi =  half;
  if(I==1) dNi = -half;
  if((I!=0)and(I!=1)) dNi = (dLn0 - dLn2)/MySqrt<Num>(four*n-two);
  return dNi;
};


//
//
//  Legendre 1-D element
//  (Line) on the reference coordinate (r)
//
//
template<typename Num, typename Iter>
class LegendreElement1D
{
  private:
    const Num Zero = Num(0.0);//Used for zeroing values etc..
    Iter pOrder=0, nSample=0; //Maximal pOrder and nSamples
    vector<Num> wi;           //Integration Weights at sample points
    vector<Num> SamplePtns;   //Sample points for numerical integration
    vector<vector<Num>> Lp;   //The legendre polynomials evaluated sample ptns
    vector<vector<Num>> dLp;  //The legendre polynomials derivatives evaluated sample ptns
    vector<vector<Num>> Phi;  //Shape function evals at sample ptns
    vector<vector<Num>> dPhi; //Shape function derivative evals at sample ptns
    GaussLegendreIntRules<Num,Iter> IntRule;

    //Generates the Legendre Polynomial evals
    //on reference element
    void Lp_eval();

    //Generates the Legendre Polynomial derivative evals
    //on reference element
    void dLp_eval();

    //Generates the Shape function evals
    //on reference element
    void Phi_eval();

    //Generates the Shape function derivative
	//evals on reference element
    void dPhi_eval();


  public:
    //Initialise the 1-D element
    LegendreElement1D(Iter pOrder_, Iter nSample_){
      SetIntegrationRule(nSample_);
      SetpOrder(pOrder_);
    };

    //Set the maximal p-Order (expensive routine)
    void SetpOrder(Iter pOrder_);

    //Set integration rule
    void SetIntegrationRule(Iter nSample_);

    //Retrieval Routines
    Num Get_Phi (Iter I, Iter Igauss){return Phi[Igauss][I];}
    Num Get_dPhi(Iter I, Iter Igauss){return dPhi[Igauss][I];};
    Num Get_xi(Iter Igauss){return SamplePtns[Igauss];};
    Num Get_wi(Iter Igauss){return wi[Igauss];};
};


//Generates the Legendre Polynomial evals
//on reference element
template <typename Num, typename Iter>
void LegendreElement1D<Num,Iter>::Lp_eval(){
  Lp.clear();
  for(Iter Isample=0; Isample<nSample; Isample++){
    Num r = SamplePtns[Isample];
    vector<Num> Lp_tmp; Lp_tmp.clear();
    for(Iter Order=0; Order<pOrder; Order++){
      Num Ln1=Zero, Ln2=Zero;
      if(Order > 1) Ln1 = Lp_tmp[Order-1];
      if(Order > 1) Ln2 = Lp_tmp[Order-2];
      Lp_tmp.push_back( LegendrePolyFunc<Num,Iter>(r, Order, Ln1, Ln2) );
    }
    Lp.push_back(Lp_tmp);
  }
};

//Generates the Legendre Polynomial derivative evals
//on reference element
template <typename Num, typename Iter>
void LegendreElement1D<Num,Iter>::dLp_eval(){
  dLp.clear();
  for(Iter Isample=0; Isample<nSample; Isample++){
    Num r = SamplePtns[Isample];
    vector<Num> dLp_tmp;
    for(Iter Order=0; Order<pOrder; Order++){
      Num Ln1=Zero, dLn1=Zero, dLn2=Zero;
      if(Order > 1) Ln1  = Lp[Isample][Order-1];
      if(Order > 1) dLn1 = dLp_tmp[Order-1];
      if(Order > 1) dLn2 = dLp_tmp[Order-2];
      dLp_tmp.push_back( LegendrePolyFuncDer<Num,Iter>(r, Order, Ln1, dLn1, dLn2) );
    }
    dLp.push_back(dLp_tmp);
  }
};

//Generates the Shape function evals
//on reference element
template <typename Num, typename Iter>
void LegendreElement1D<Num,Iter>::Phi_eval(){
  Phi.clear();
  for(Iter Isample=0; Isample<nSample; Isample++){
    Num r = SamplePtns[Isample];
    vector<Num> Phi_tmp;
    for(Iter Order=0; Order<pOrder; Order++){
      Num Ln2=Zero, Ln0=Zero;
      if(Order > 1) Ln0 = Lp[Isample][Order];
      if(Order > 1) Ln2 = Lp[Isample][Order-2];
      Phi_tmp.push_back( LegendreShapeFunc<Num,Iter>(r, Order, Ln0, Ln2) );
    }
    Phi.push_back(Phi_tmp);
  }
};

//Generates the Shape function derivative
//evals on reference element
template <typename Num, typename Iter>
void LegendreElement1D<Num,Iter>::dPhi_eval(){
  dPhi.clear();
  for(Iter Isample=0; Isample<nSample; Isample++){
    vector<Num> dPhi_tmp;
    for(Iter Order=0; Order<pOrder; Order++){
      Num dLn2=Zero, dLn0=Zero;
      if(Order > 1) dLn0 = dLp[Isample][Order];
      if(Order > 1) dLn2 = dLp[Isample][Order-2]; 
      dPhi_tmp.push_back( LegendreShapeFuncDer<Num,Iter>(Order, dLn0, dLn2) );
    }
    dPhi.push_back(dPhi_tmp);
  }
};

//Set the maximal p-Order (expensive routine)
template <typename Num, typename Iter>
void LegendreElement1D<Num,Iter>::SetpOrder(Iter pOrder_){
  pOrder = pOrder_;
  Lp_eval();
  dLp_eval();
  Phi_eval();
  dPhi_eval();
};

//Set integration rule
template <typename Num, typename Iter>
void LegendreElement1D<Num,Iter>::SetIntegrationRule(Iter nSample_){
  nSample = nSample_;
  for(Iter I=0; I<nSample; I++){
    wi.push_back(IntRule.GL_weight(I,nSample));
    SamplePtns.push_back(IntRule.GL_point(I,nSample));
  }
};
#endif
