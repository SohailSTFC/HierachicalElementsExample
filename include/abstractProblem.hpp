#ifndef ABSTRACTPROBLEM_HPP
#define ABSTRACTPROBLEM_HPP

//
// Standard C++ libraries
//
#include <vector>
#include <array>


/**************************************\
! Abstract Problem Class
!
! this class acts as a template for
! solving problems using Hierachical 
! Finite Elements.
!
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
class abstractProblem{
  private:


  public:
    //The number of degrees
    //of freedom on each entity
    //type
    virtual int getNODOFVert() = 0;
    virtual int getNODOFEdge() = 0;
    virtual int getNODOFFace() = 0;
    virtual int getNODOFVolm() = 0;

    virtual void calcQpResidual() = 0;
    virtual void calcQpJacobian() = 0;
};
#endif