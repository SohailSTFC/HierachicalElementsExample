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
! To support h-FEM/AMR multiple of these
! classes must be generated
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
    std::vector<std::array<Integer,3>> EFV_NDOFs; //EFV-(Edge-Face-Volume)

    //Mappings, Iterators and distribution
    //of mesh/grid
    const amrex::BoxArray&     nba;
    amrex::DistributionMapping dm;

    //Residual and solution
    //vectors/multiFabs
    amrex::MultiFab  VertsSolution;
    amrex::MultiFab  VertsResidual;
    std::vector<std::array<amrex::MultiFab,3>> EFVsSolution; //EFV-(Edge-Face-Volume)
    std::vector<std::array<amrex::MultiFab,3>> EFVsResidual; //EFV-(Edge-Face-Volume)

    abstractProblem *problem; //The abstract problem your solving
//  Solver *pAMGSolver;

    void CalcEFVDOFs();
    void Set_pFEMMultiFabs();
  public:
    //Class constructor
    hierachicalFEMeshProblemSolver(abstractProblem *prob, unsigned MAXpOrder_, unsigned MINpOrder_);

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
                                                        , unsigned MAXpOrder_=2
                                                        , unsigned MINpOrder_=1)
{
  if(MAXpOrder_== 0)        throw("Zero order elements not supported MAXp");
  if(MINpOrder_== 0)        throw("Zero order elements not supported MINp");
  if(MINpOrder_>MAXpOrder_) throw("MAXpOrder less than MINpOrder");
  MINpOrder = MINpOrder_;
  MAXpOrder = MAXpOrder_;

  problem = new abstractProblem(*prob);
  CalcEFVDOFs();
  Set_pFEMMultiFabs();
};


// Calculate the number of DOFs
// for the Edges-Faces-Volumes
// at each polynomial order
// level
//
template<typename Integer, typename RealNum, Integer DIM>
void hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::CalcEFVDOFs(){
  EFV_NDOFs.clear();
  if(MAXpOrder>1){//Only quadratic or higher order have these
    for(int I=2; I<=MAXpOrder; I++){
      Integer ndofE=0, ndofF=0, ndofV=0; //The total dofs associated with higher order
      Integer dnodE=0, dnodF=0, dnodV=0; //The nodes associated with higher order
      if(DIM==1) dnodE=1;
      if(DIM==2) dnodE=4;
      if(DIM==3) dnodE=12;

      if(DIM==2) dnodF=2*I - 3;
      if(DIM==3) dnodF=12*I - 18;

      if(DIM==3) dnodV=3*I*I - 9*I + 7;

      ndofE = dnodE*(problem->getNODOFEdge());
      ndofF = dnodF*(problem->getNODOFFace());
      ndofV = dnodV*(problem->getNODOFVolm());
      EFV_NDOFs.push_back({ndofE,ndofF,ndofV});
    }
  }
};


// Constructs the MultiFabs for each
// polynomial order being used to
// solve the problem
//
template<typename Integer, typename RealNum, Integer DIM>
  void hierachicalFEMeshProblemSolver<Integer,RealNum,DIM>::Set_pFEMMultiFabs();


    amrex::MultiFab  VertsSolution;
    amrex::MultiFab  VertsResidual;
    std::vector<std::array<amrex::MultiFab,3>> EFVsSolution; 
    std::vector<std::array<amrex::MultiFab,3>> EFVsResidual;
};
#endif