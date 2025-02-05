#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

//Standard C++ libraries
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>

//
// My libraries
//
#include "include/hierachicalFiniteElement.hpp"
#include "include/hierachicalFEMesh.hpp"
#include "include/integrationRules.hpp"

//
// Problem Operators
//
#include "sampleProblems/poissonProb.hpp"

int main(int argc, char* argv[]){

  const int pOrder_1D=3;
  const int nsampl=4;
  std::vector<std::array<amrex::Real,nsampl>> Ni(pOrder_1D+1), dNi(pOrder_1D+1);
  std::vector<std::array<amrex::Real,2>> XiWi(nsampl);
  CalcGaussLegendreXW1D<int,amrex::Real,nsampl>(&XiWi);
  HLegendreElement1D::CalcSF_SDFs1D<int,amrex::Real,nsampl,pOrder_1D>(&Ni, &dNi, XiWi);


  for(int I=0; I<nsampl; I++){
    std::cout << std::setw(12) << XiWi[I][0]
              << std::setw(12) << XiWi[I][1] << std::endl;
  }
  for(int I=0; I<nsampl; I++){
    for(int J=0; J<(pOrder_1D+1); J++){
      std::cout << std::setw(12) << Ni[J][I];
    }
    std::cout << std::endl;
  }


  amrex::Initialize(argc, argv);
  {
    amrex::Print() << "Hi AMReX version " << amrex::Version() << "\n";


    poissonProblem MyPoisson(AMREX_SPACEDIM,"Config.txt");
    hierachicalFEMeshProblemSolver<int,amrex::Real,AMREX_SPACEDIM>  MySolver(&MyPoisson);

/*

    //Data for MultiFab
    amrex::RealBox  real_box({0.,0.,0.},{1.,1.,1.});//Physical Dimension of domain
    amrex::Geometry geom(domain, &real_box);        //Geometric object
    amrex::GpuArray<amrex::Real, 3> dx = geom.CellSizeArray();


    const amrex::MultiArray4<amrex::Real>& mf_arrs = mf.arrays();
    const amrex::IntVect ngs(nGhost);
    amrex::ParallelFor(mf, ngs, [=] AMREX_GPU_DEVICE( int nbx, int I, int J, int K) noexcept {
      amrex::Real x = I * dx[0];
      amrex::Real y = J * dx[1];
      amrex::Real z = K * dx[2];
      amrex::Real rSqaured = amrex::Real(EntityMapping::EntityTypeTest(I,J,K,pOrder_1D));
      for(int L=0; L<ncomp; L++){//Iterate over components
        mf_arrs[nbx](I,J,K,L)  = rSqaured + amrex::Real(L);//1.0 + std::exp(-r_squared);
      };
    });


    //Plot the MultiFab Data
    WriteSingleLevelPlotfile("plt001",mf,{"comp0","comp1","comp2"},geom,0.,0);
*/

  }
  amrex::Finalize();
  return 0;
};