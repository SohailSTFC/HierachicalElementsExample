#include <vector>
#include <iomanip>
#include <iostream>
#include "Legendre1DElement.hpp"
#include "LegendreTensorElement.hpp"

using namespace std;
int main(){
  const int ndim = 1;
  int pOrder = 6;
  int nsamples = 5;
  LegendreTensorElement<float,int> MyRefElm(pOrder, nsamples, ndim);

  for(int I=0; I<nsamples; I++){
  //  for(int J=0; J<nsamples; J++){
    //  int Igauss[ndim] = {I,J};
      int Igauss[ndim] = {I};
      cout <<  setw(10) << MyRefElm.Get_wi(Igauss) << endl;
   // }
  }

  return 0;
}