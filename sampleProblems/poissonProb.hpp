#ifndef POISSONPROBLEM_HPP
#define POISSONPROBLEM_HPP

#include <string>

#include "../include/abstractProblem.hpp"
/**************************************\
! Poisson Problem Class
!
! A sample problem that solves the
! poisson problem in the form
!   grad(u) = f
! This is discretised using H1 elements
! with a legendre polynomial basis
! (only one offered currentl ;-) )
!
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
class poissonProblem : public abstractProblem{
  private:
    const int nodofVt=1, nodofEd=1, nodofFa=1, nodofVo=1;
    
  public:
    poissonProblem(int DIM, std::string ConfigFile);

    //The number of degrees
    //of freedom on each entity
    //type
    virtual int getNODOFVert();
    virtual int getNODOFEdge();
    virtual int getNODOFFace();
    virtual int getNODOFVolm();

    virtual void calcQpResidual();

    virtual void calcQpJacobian();
};

poissonProblem::poissonProblem(int DIM, std::string ConfigFile){};

int poissonProblem::getNODOFVert(){return nodofVt;};
int poissonProblem::getNODOFEdge(){return nodofEd;};
int poissonProblem::getNODOFFace(){return nodofFa;};
int poissonProblem::getNODOFVolm(){return nodofVo;};


void poissonProblem::calcQpResidual(){
};


void poissonProblem::calcQpJacobian(){
};
#endif