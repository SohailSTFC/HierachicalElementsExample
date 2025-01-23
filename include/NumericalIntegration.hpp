#ifndef NUMERICALINTEGRATION_HPP
#define NUMERICALINTEGRATION_HPP 
#include <iostream>
using namespace std;

//
//
// Returns the 1-D sample points
// and weights for a Gauss-Legendre
// rule of Order N
//
//
template<typename Num, typename Iter>
class GaussLegendreIntRules{
  private: 
    const Num  Xi_g1[1] = { Num(0.000000)};
    const Num  Wi_g1[1] = { Num(2.000000)};

    const Num  Xi_g2[2] = { Num(-0.577350), Num(0.577350)};
    const Num  Wi_g2[2] = { Num( 1.000000), Num(1.000000)};

    const Num  Xi_g3[3] = { Num(-0.774597), Num( 0.000000), 0.774597};
    const Num  Wi_g3[3] = { Num( 0.555556), Num( 0.888889), 0.555556};

    const Num  Xi_g4[4] = {-0.861136,-0.339981, 0.339981, 0.861136};
    const Num  Wi_g4[4] = { 0.347855, 0.652145, 0.652145, 0.347855};

    const Num  Xi_g5[5] = {-0.906180,-0.538469, 0.000000, 0.538469, 0.906180};
    const Num  Wi_g5[5] = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927};

  public:
    GaussLegendreIntRules(){};

    Num GL_point(Iter I, Iter nSamples){
      if((I>=nSamples)or(I<0)) cout << "Error out I: of bounds" << endl;
      if((I>=nSamples)or(I<0)) return -1.0;
      if(nSamples==1) Xi_g1[I];
      if(nSamples==2) Xi_g2[I];
      if(nSamples==3) Xi_g3[I];
      if(nSamples==4) Xi_g4[I];
      if(nSamples==5) Xi_g5[I];
      return -1.0;
    };

    Num GL_weight(Iter I, Iter nSamples){
      if((I>=nSamples)or(I<0)) cout << "Error out I: of bounds" << endl;
      if((I>=nSamples)or(I<0)) return -1.0;
      if(nSamples==1) return Wi_g1[I];
      if(nSamples==2) return Wi_g2[I];
      if(nSamples==3) return Wi_g3[I];
      if(nSamples==4) return Wi_g4[I];
      if(nSamples==5) return Wi_g5[I];
      return -1.0;
    };
};
#endif