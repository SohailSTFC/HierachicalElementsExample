#ifndef LEGENDRETENSORDELEMENT_HPP
#define LEGENDRETENSORDELEMENT_HPP 
#include "Legendre1DElement.hpp"
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//
//
//  Legendre n-D element
//  (Aribitrary dimension) on the reference coordinate (r^n)
//  for now order in all dimensions is the same
//
template<typename Num, typename Iter>
class LegendreTensorElement
{
  private:
    const Num One  = Num(1.0);    //Used for one values etc..
    const Num Zero = Num(0.0);    //Used for zeroing values etc..
    Iter ndim=0;                  //The problem dimension
    Iter pOrder1D=0, nSample1D=0; //Maximal pOrder and nSamples
    Iter nSample=0;               //Total sample points
    LegendreElement1D<Num,Iter> *Ref1DElm;  //Reference 1-D element
  public:
    //Initialise the n-D element
    LegendreTensorElement(Iter pOrder1D_, Iter nSample1D_, Iter ndim_)
    {
      ndim = ndim_;
      pOrder1D=pOrder1D_; 
	  nSample1D=nSample1D_;
      nSample = nSample1D; 
      if(ndim>1) for(Iter I=1; I<ndim; I++) nSample = nSample*nSample1D; //(Maybe use symmetry to reduce?) 
      Ref1DElm = new LegendreElement1D<Num,Iter>(pOrder1D_, nSample1D_);
    };

    //Set the maximal p-Order (expensive routine)
    void SetpOrder(Iter pOrder_);

    //Set integration rule
    void SetIntegrationRule(Iter nSample_);

    //Retrieval Routines
    Num Get_wi  (Iter Igauss[]);
    Num Get_Phi (Iter I[], Iter Igauss[]);
    Num Get_dPhi(Iter I[], Num J[], Iter Igauss[]); //J is an zero array with only 1 entry equal to 1
};

//Set the maximal p-Order (expensive routine)
template <typename Num, typename Iter>
void LegendreTensorElement<Num,Iter>::SetpOrder(Iter pOrder_){
  Ref1DElm->SetpOrder(pOrder_);
};

//Set integration rule
template <typename Num, typename Iter>
void LegendreTensorElement<Num,Iter>::SetIntegrationRule(Iter nSample_){
  Ref1DElm->SetIntegrationRule(nSample_);
};

//Calculate the n-D gauss integration weight
template <typename Num, typename Iter>
Num LegendreTensorElement<Num,Iter>::Get_wi(Iter Igauss[]){
  Num wi = Ref1DElm->Get_wi(Igauss[0]);
  if(ndim != 1) for(Iter J=1; J<ndim; J++) wi = wi*(Ref1DElm->Get_wi(Igauss[0]));
  return wi;
};

//Calculate the n-D shape function weight
template <typename Num, typename Iter>
Num LegendreTensorElement<Num,Iter>::Get_Phi (Iter I[], Iter Igauss[]){
  Num Phi = (Ref1DElm->Get_Phi(I[0], Igauss[0]));
  if(ndim != 1) for(Iter J=1; J<ndim; J++) Phi = Phi*( Ref1DElm->Get_Phi(I[J],Igauss[J]) );
  return Phi;
};

//Calculate the n-D shape function derivative weight (this avoids if statements)
template <typename Num, typename Iter>
Num LegendreTensorElement<Num,Iter>::Get_dPhi(Iter I[], Num J[], Iter Igauss[]){
  Num dPhi = J[0]*(Ref1DElm->Get_dPhi(I[0], Igauss[0])) + (One - J[0])*( Ref1DElm->Get_Phi(I[0],Igauss[0]) );
  if(ndim != 1){ for(Iter K=1; K<ndim; K++){
    Num dPhi_tmp = J[0]*(Ref1DElm->Get_dPhi(I[0], Igauss[0])) + (One - J[0])*( Ref1DElm->Get_Phi(I[0],Igauss[0]) );
    dPhi = dPhi*dPhi_tmp;
  }}
  return dPhi;
};
#endif
