#ifndef HIERACHICALFEMESH_HPP
#define HIERACHICALFEMESH_HPP

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


//
// My libraries needed
// for completeness
//
//
#include "abstractProblem.hpp"

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;

/**************************************\
! Hierachical FE-Mesh problem solver 
! class
!
! Generates a mesh like structure for
! Lengendre basis type finite elements
! this type of element has upto 4-entity
! types depending on the problem dimension:
!
! Vertices : (1-3) D
! Edges    : (1-3) D
! Faces    : (2-3) D
! Volumes  : 3-D
!
! This class is a bit problem specific
! in order to generate sufficient memory
! for the FE-problem being solved as
! such the problem defined must have 
! enough DOFs defined for the entirety 
! of the problem.
!
! To support hp-FEM/AMR multiple of these
! classes must be generated
!
!
! Uses element-by-element method avoiding
! avoiding global assembly
!
!
\**************************************/
template<typename Integer, typename RealNum, Integer DIM>
class hierachicalFEMeshProblemSolver{
  private:
    const unsigned nGhost=1; //Luckily for FE-analysis 
                             //only need 1-Ghost/halo-Cell
                             //thickness per processor domain

    //Number of DOFs per node for
    //each polynomial order
    unsigned MINpOrder;
    unsigned MAXpOrder;
    Integer Vert_NDOFs;
    std::array<Integer,3> EFV_NDOFs; //EFV-(Edge-Face-Volume)

    //Mappings, Iterators and distribution
    //of mesh/grid
    amrex::BoxArray            *nba;
    amrex::DistributionMapping dm;

    //Residual and solution
    //vectors/multiFabs
    amrex::MultiFab  *VertsSolution;
    amrex::MultiFab  *VertsResidual;
    std::array<amrex::MultiFab*,3> EFVsSolution; //EFV-(Edge-Face-Volume)
    std::array<amrex::MultiFab*,3> EFVsResidual; //EFV-(Edge-Face-Volume)

    abstractProblem *problem; //The abstract problem your solving
//  Solver *pAMGSolver;

    void Calc_dEFVDOFs();     //Calculates the number of DOFs
    void Set_pFEMMultiFabs(); //Sets the MultiFabs
  public:
    //Class constructor
    hierachicalFEMeshProblemSolver(abstractProblem *prob
                                 , unsigned MAXpOrder_
                                 , unsigned MINpOrder_);

    hierachicalFEMeshProblemSolver(abstractProblem *prob);

    //Class destructor
    ~hierachicalFEMeshProblemSolver();

    //Solve the equation
    void solveEquation();

    //Plot data
    void plotData();
};


// Construct the mesh class
// for a hierachical finite
// element problem
//
template<typename Integer, typename RealNum, Integer DIM>
hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::hierachicalFEMeshProblemSolver(abstractProblem *prob
                                                                                  , unsigned MAXpOrder_
                                                                                  , unsigned MINpOrder_)
{
  if(MAXpOrder_== 0)        throw("Zero order elements not supported MAXp");
  if(MINpOrder_== 0)        throw("Zero order elements not supported MINp");
  if(MINpOrder_>MAXpOrder_) throw("MAXpOrder less than MINpOrder");
  MINpOrder = MINpOrder_;
  MAXpOrder = MAXpOrder_;

  problem = prob;
  Calc_dEFVDOFs();
/******************************************\
  This is a hack to generate a mesh/grid
  quickly in AMReX just to see if the other
  functions in the class work
\******************************************/
    int n_cell=16;
    int max_grid_size=32;

    amrex::IntVect domElm_lo(0,0,0);
    amrex::IntVect domElm_hi(n_cell-1,n_cell-1,n_cell-1);

    amrex::Box domain(domElm_lo,domElm_hi);
    amrex::BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    nba = new amrex::BoxArray(amrex::convert(ba,amrex::IntVect::TheNodeVector()));
    dm = amrex::DistributionMapping(*nba);
/******************************************\
\******************************************/
  Set_pFEMMultiFabs();
};

// Construct the mesh class
// for a hierachical finite
// element problem (default setting)
//
template<typename Integer, typename RealNum, Integer DIM>
hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::hierachicalFEMeshProblemSolver(abstractProblem *prob)
 :hierachicalFEMeshProblemSolver(prob, 2, 1){};

// Destruct the mesh class
// for a hierachical finite
// element problem
//
template<typename Integer, typename RealNum, Integer DIM>
hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::~hierachicalFEMeshProblemSolver()
{
  //Solution Vectors
  if(Vert_NDOFs > 0) delete VertsSolution;
  for(int I=0; I<3; I++)
    if(EFV_NDOFs[I] > 0) delete EFVsSolution[I];

  //Residual Vectors
  if(Vert_NDOFs > 0) delete VertsResidual;
  for(int I=0; I<3; I++)
    if(EFV_NDOFs[I] > 0) delete EFVsResidual[I];
};


// Calculate the number of DOFs
// for the Edges-Faces-Volumes
// at each polynomial order
// level
//
template<typename Integer, typename RealNum, Integer DIM>
void hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::Calc_dEFVDOFs()
{
  Vert_NDOFs = std::pow(2,DIM)*(problem->getNODOFVert());

  //The total dofs associated with higher order
  Integer ndofE=0, ndofF=0, ndofV=0;
  if(MAXpOrder>1){//Only quadratic or higher order have these
    for(int I=MINpOrder; I<=MAXpOrder; I++){
      //The change in nodes associated with higher order
      Integer dnodE=0, dnodF=0, dnodV=0;
      if(DIM==1) dnodE=1;
      if(DIM==2) dnodE=4;
      if(DIM==3) dnodE=12;

      if(DIM==2) dnodF=2*I - 3;
      if(DIM==3) dnodF=12*I - 18;

      if(DIM==3) dnodV=3*I*I - 9*I + 7;

      ndofE += dnodE*(problem->getNODOFEdge());
      ndofF += dnodF*(problem->getNODOFFace());
      ndofV += dnodV*(problem->getNODOFVolm());
    }
  }
  EFV_NDOFs[0] = ndofE;
  EFV_NDOFs[1] = ndofF;
  EFV_NDOFs[2] = ndofV;

  if(Vert_NDOFs<0) throw std::invalid_argument("The number of Vertex DOFs is less than 0");
};


// Constructs the MultiFabs for all
// polynomial orders being used to
// solve the problem for each entity
// type
//
template<typename Integer, typename RealNum, Integer DIM>
void hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::Set_pFEMMultiFabs()
{
  //Solution Vectors
  if(Vert_NDOFs > 0) VertsSolution = new amrex::MultiFab(*nba, dm, Vert_NDOFs, nGhost);
  for(int I=0; I<3; I++)
    if(EFV_NDOFs[I] > 0)
      EFVsSolution[I] = new amrex::MultiFab(*nba, dm, EFV_NDOFs[0], nGhost);

  //Residual Vectors
  if(Vert_NDOFs > 0) VertsResidual = new amrex::MultiFab(*nba, dm, Vert_NDOFs, nGhost);
  for(int I=0; I<3; I++)
    if(EFV_NDOFs[I] > 0)
      EFVsResidual[I] = new amrex::MultiFab(*nba, dm, EFV_NDOFs[0], nGhost);
};
#endif
