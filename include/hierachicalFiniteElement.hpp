#ifndef HIERACHICALFINITEELEMENT_HPP
#define HIERACHICALFINITEELEMENT_HPP

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
! EntityType markers
!  Vertex 0  Permutation(Node-Node-Node)
!  Edge   1  Permutation(Node-Node-Edge)
!  Face   2  Permutation(Node-Edge-Edge)
!  Volume 3  Permutation(Edge-Edge-Edge)
!
!
\**************************************/
namespace EntityMapping{
  //This mapper assumes a the IJK iterator
  //starts at [0,0,0]
  AMREX_FORCE_INLINE
  int EntityTypeTest(int I, int J, int K, int pOrder){
    return (( I%(pOrder+1) )!=0)
         + (( J%(pOrder+1) )!=0)
         + (( K%(pOrder+1) )!=0);
  };

  AMREX_FORCE_INLINE
  std::array<int,3> EntityIterator(int I, int J, int K, int pOrder){
    return { (I+1)%pOrder, (J+1)%pOrder, (K+1)%pOrder };
  };

};

/**************************************\
! Tensor element N-D
!
! Takes a 1-D element basis and forms
! higher order basis functions, mappings
! supported for Tensor elements e.g.
! lines, quads, hexs, etc...
! as arguments it takes Legendre
! 1-D polynomial modal shape functions
! and derivative
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
namespace TensorElementND{
  //
  //  Form a higher dimension shape
  //  function weight using a 1-D shape
  //  function weight and the DOF position
  //
  template<typename Integer, typename RealNum, Integer nSampl, Integer DIM>
  AMREX_FORCE_INLINE
  RealNum Ni_ND(const std::vector<std::array<RealNum,nSampl>> Ni1D
              , const std::array<Integer,DIM> IJK       //The Node
              , const std::array<Integer,DIM> IJKGauss) //The Gauss-point
  {
    //Calculate shape function
    RealNum Ni = Ni1D[IJK[0]][IJKGauss[0]];
    #pragma unroll
    for(Integer I=1; I<DIM; I++) Ni = Ni*Ni1D[IJK[I]][IJKGauss[I]];
	return Ni;
  };


  //
  //  Form a higher dimension shape function
  //  derivative weight using the 1-D shape
  //  function weight, derivative weight,
  //  DOF position and the Gauss-Point
  //
  template<typename Integer, typename RealNum, Integer nSampl, Integer DIM>
  AMREX_FORCE_INLINE
  std::array<RealNum,DIM> dNi_ND(const std::vector<std::array<RealNum,nSampl>> Ni1D
                               , const std::vector<std::array<RealNum,nSampl>> dNi1D
                               , const std::array<Integer,DIM> IJK       //The Node
                               , const std::array<Integer,DIM> IJKGauss) //The Gauss-point
  {
    RealNum Zero = RealNum(0.0);
    RealNum One  = RealNum(1.0);
    std::array<RealNum,DIM> dNi, mask;

    //Calculate shape function derivative
    #pragma unroll
    for(Integer I=0; I<DIM; I++){
      #pragma unroll
      for(Integer J=0; J<DIM; J++) mask[J]=Zero;
      mask[I] = One;

      dNi[I] = dNi[I]*( (One - mask[0])*Ni1D[IJK[0]][IJKGauss[0]] 
                              + mask[0]*dNi1D[IJK[0]][IJKGauss[0]] );
      #pragma unroll
      for(Integer J=1; J<DIM; J++){
        dNi[I] = dNi[I]*( (One - mask[J])*Ni1D[IJK[J]][IJKGauss[J]] 
                               + mask[J]*dNi1D[IJK[J]][IJKGauss[J]] );
      }
    }
    return dNi;
  };

  //
  //  For effiency purposes one may want to
  //  store the entire shape-function and shape
  //  function derivative weights array for 
  //  an N-D element so it can be accessed as
  //  oppose to being recalculated on every element
  //  in these cases this function can be used
  //
  template<typename Integer, typename RealNum, Integer nSampl, Integer DIM>
  AMREX_FORCE_INLINE
  void CalcND_SF_SFD(std::vector<std::array<RealNum,nSampl>>                 NiND
                     std::vector<std::array<std::array<RealNum,DIM>,nSampl>> dNiND
                   , const std::vector<std::array<RealNum,nSampl>>           Ni1D
                   , const std::vector<std::array<RealNum,nSampl>>           dNi1D)
  {

  //TODO
/*    for(Integer I=0; I<; I++){
      NiND
      dNiND
    }*/
  };


  //
  // Assume the element shape Jacobian matrix
  // dni/dxj is only a function of the vertex
  // contributions this kind of means that the
  // that the other shape functions are just
  // affine transformations of their modal
  // coord-Shape functions
  //
  template<typename Integer, typename RealNum, Integer DIM>
  AMREX_FORCE_INLINE
  void calcElmJacDet_JacInv(Eigen::MatrixXd *JacInv
                          , RealNum         *JacDet
                          , const Integer Isample
                          , const Integer nVerts
                          , const std::vector<std::array<std::array<RealNum,DIM>,nSampl>>  dNiND
                          , const std::vector<std::array<RealNum,DIM>>  coord){
    Eigen::MatrixXd Jac = Eigen::MatrixXd::Zero(DIM, DIM);
    for(int I=0; I<DIM; I++){
      for(int J=0; J<DIM; J++){
        for(int K=0; K<nVerts; K++){
          Jac(I,J) += dNiND[I][Isample][K]*coord[J][K];
        }
      }
    }
    *JacDet = Jac.determinant();
    *JacInv = Jac.inverse();
  }
};//End of TensorElementND namespace


/**************************************\
! Hierachical Legendre element 1-D
! in modal space
!
! This generates the shape
! function (Ni) and shape function
! derivative (dNi) weights for a 1-D
! finite element and returns them
! as a pair of arrays over a span of
! input set of sample points
! it is templated class for the
! choice of precision and integer size
! so that it can be used potentially
! for other applications and libraries
!
! This documentation was used as
! reference
! http://mofem.eng.gla.ac.uk/mofem/html/hierarchical_approximation_1.html#mjx-eqn-eqElementDOF3D
!
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
namespace HLegendreElement1D{
  //
  // Calculate the Legendre
  // polynomial function
  //
  template<typename Integer, typename RealNum>
  AMREX_FORCE_INLINE
  RealNum LegendrePolyFunc(RealNum r, Integer I, RealNum Ln1, RealNum Ln2){
    RealNum Lp  = RealNum(0.0);
    RealNum One = RealNum(1.0);
    RealNum Two = RealNum(2.0);
    RealNum N   = RealNum(I);

    if(I==0)           Lp = One;
    if(I==1)           Lp = r;
    if((I!=0)and(I!=1))Lp = (One/N)*( (Two*N-One)*r*Ln1 - (N-One)*Ln2 );
    return Lp;
  };


  //
  // Calculate the Legendre
  // polynomial function derivative
  //
  template<typename Integer, typename RealNum>
  AMREX_FORCE_INLINE
  RealNum LegendrePolyFuncDer(RealNum r, Integer I, RealNum Ln1
                            , RealNum dLn1, RealNum dLn2){
    RealNum Zero = RealNum(0.0);
    RealNum one  = RealNum(1.0);
    RealNum two  = RealNum(2.0);
    RealNum n    = RealNum(I);
    RealNum dLP  = Zero;
    if(I==0)            dLP = Zero;
    if(I==1)            dLP = one;
    if((I!=0)and(I!=1)) dLP = (one/n)*( (two*n - one)*(Ln1 + r*dLn1) - (n - one)*dLn2 );
    return dLP;
  };


  //
  // Legendre type shape
  // function calculator
  //
  template<typename Integer, typename RealNum>
  AMREX_FORCE_INLINE
  RealNum LegendreShapeFunc(RealNum r, Integer I, RealNum Ln0, RealNum Ln2){
    RealNum Ni   = RealNum(0.0);
    RealNum half = RealNum(0.5);
    RealNum one  = RealNum(1.0);
    RealNum two  = RealNum(2.0);
    RealNum four = RealNum(4.0);
    RealNum n    = RealNum(I);
    if(I==0)            Ni = half*(one + r);
    if(I==1)            Ni = half*(one - r);
    if((I!=0)and(I!=1)) Ni = (Ln0 - Ln2)/std::sqrt(four*n-two);
    return Ni;
  };


  //
  // Legendre type shape
  // function derivative
  // calculator
  //
  template<typename Integer, typename RealNum>
  AMREX_FORCE_INLINE
  RealNum LegendreShapeFuncDer(Integer I, RealNum dLn0, RealNum dLn2){
    RealNum dNi  = RealNum(0.0);
    RealNum half = RealNum(0.5);
    RealNum two  = RealNum(2.0);
    RealNum four = RealNum(4.0);
    RealNum n    = RealNum(I);
    if(I==0) dNi =  half;
    if(I==1) dNi = -half;
    if((I!=0)and(I!=1)) dNi = (dLn0 - dLn2)/std::sqrt(four*n-two);
    return dNi;
  };


  //
  // Legendre-Basis Spectral Element
  //
  // Calculate the shape function
  // and derivative weights for a 1-D
  // element on a given set of sample
  // integration points, in the
  // reference element space
  //
  template<typename Integer, typename RealNum, Integer nSampl, Integer pOrder>
  AMREX_FORCE_INLINE
  void CalcSF_SDFs1D(std::vector<std::array<RealNum,nSampl>> *Ni
                   , std::vector<std::array<RealNum,nSampl>> *dNi
                   , const std::vector<std::array<RealNum,2>> XiWi)
  {
    if(nSampl <= 0)              throw std::invalid_argument("nSampl is smaller zero");
    if(XiWi.size() != nSampl)    throw std::invalid_argument("XiWi are not of size nsampl");
    if(Ni->size() != (pOrder+1)) throw std::invalid_argument("Ni is not of size pOrder+1");
    if(dNi->size() != (pOrder+1))throw std::invalid_argument("dNi is not of size pOrder+1");
    std::vector<std::array<RealNum,nSampl>> Lp, dLp;
    Lp.clear(); dLp.clear();

    //Calculate the Legendre polynomials
    //and their derivatives as well as
    //the Shape function basis and
    //derivatives
    for(Integer J=0; J<(pOrder+1); J++){
      std::array<RealNum,nSampl> arr_empt;
      Lp.push_back(arr_empt);
      dLp.push_back(arr_empt);
      for(Integer I=0; I<nSampl; I++){
        RealNum Ln1=0.0, Ln2=0.0, dLn1=0.0, dLn2=0.0;
        if(J > 1) Ln1  = Lp[J-1][I];
        if(J > 1) Ln2  = Lp[J-2][I];
        if(J > 1) dLn1 = dLp[J-1][I];
        if(J > 1) dLn2 = dLp[J-2][I];
        Lp[J][I]  = LegendrePolyFunc<Integer,RealNum>(XiWi[I][0], J, Ln1, Ln2);
        dLp[J][I] = LegendrePolyFuncDer<Integer,RealNum>(XiWi[I][0], J, Ln1, dLn1, dLn2);
        (*Ni)[J][I]  = LegendreShapeFunc<Integer,RealNum>(XiWi[I][0], J, Lp[J][I], Ln2);
        (*dNi)[J][I] = LegendreShapeFuncDer<Integer,RealNum>(J, dLp[J][I], dLn2);
      }
    }
    Lp.clear(); dLp.clear();
  };

};//End of HLegendreElement namespace
#endif
