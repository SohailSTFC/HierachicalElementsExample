#include <vector>
#include <iomanip>
#include <iostream>
#include "include/Legendre1DElement.hpp"
#include "include/LegendreTensorElement.hpp"

using namespace std;
int main(){
  const int ndim = 3;
  int pOrder = 4;
  int nsamples = 3;
  LegendreTensorElement<float,int> MyRefElm(pOrder, nsamples, ndim);

  for(int I=0; I<nsamples; I++){
    for(int J=0; J<nsamples; J++){
      for(int K=0; K<nsamples; K++){
        int Igauss[ndim] = {I,J,K};
        cout <<  setw(10) << MyRefElm.Get_wi(Igauss) << endl;
      }
    }
  }

  return 0;
}
