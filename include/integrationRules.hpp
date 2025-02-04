#ifndef INTEGRATIONRULES_HPP
#define INTEGRATIONRULES_HPP

//
// Standard C++ libraries
//
#include <vector>
#include <array>

//
// AMReX libraries for 
//
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>

//
// Eigen Linear algebra
// libraries for dense matrix
// arithmatic
//
//
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/**************************************\
! n-th Order Gauss-Legendre Integration
! rules
! This generates the integration sample
! points and weights of the Gauss-
! Legendre integration rules for an
! arbitrary number of sample points on
! the fly (This routine is really only
! for use once per process before the
! FE-integration starts to generate
! the integration tables.
!
! This uses the Golub-Welsch algorithm
! using libEigen for dense matrix 
! arithmetic
!
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
//
// Sample points and integration
// weights for the Gauss-Legendre
// integration scheme
//
template<typename Integer, typename RealNum, Integer nSampl>
void CalcGaussLegendreXW1D(std::vector<std::array<amrex::Real,2>> *XiWi){
  if(nSampl <= 0) throw std::invalid_argument("nsampl must be a positive integer.");
  if(XiWi->size() != nSampl) throw std::invalid_argument("XiWi are not of size nsampl");

  // Some useful constants
  // used in the function
  RealNum eps  = RealNum(1.0e-15);
  RealNum zero = RealNum(0.0);
  RealNum one  = RealNum(1.0);
  RealNum Two  = RealNum(2.0);
  RealNum four = RealNum(4.0);

  //Calculates an stores the symetric
  //tridiagonal Jacobi matrix using only
  //2-entries
  // 0:Diagonal entry
  // 1:Offdiagonal entry
  Eigen::MatrixXd JacobiMat = Eigen::MatrixXd::Zero(nSampl, nSampl);

  #pragma unroll
  for(Integer I=0; I<nSampl-1; I++){
    RealNum i       = RealNum(I+1);
    RealNum Offdiag = i/std::sqrt(four*i*i - one);
    // Create the symmetric tridiagonal Jacobi matrix
    JacobiMat(I, I+1) = Offdiag;
    JacobiMat(I+1, I) = Offdiag;
  }

  // Calculate the Nodes
  // and weights by solving
  // the eigenvalue problem
  // on the Jacobi-Matrix

  // Compute eigenvalues and eigenvectors
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(JacobiMat);
  if(eigen_solver.info() != Eigen::Success)
    throw std::runtime_error("Eigenvalue decomposition failed.");

  Eigen::VectorXd eigenvalues = eigen_solver.eigenvalues();
  Eigen::MatrixXd eigenvectors = eigen_solver.eigenvectors();

  // Calculate the Gauss points
  // and weights
  #pragma unroll
  for(Integer I=0; I<nSampl; I++){
    (*XiWi)[I][0] = eigenvalues(I);
    if(abs(eigenvalues(I))<eps) (*XiWi)[I][0] = zero;
    (*XiWi)[I][1] = Two * std::pow(eigenvectors(0, I), 2);
  }
};
#endif