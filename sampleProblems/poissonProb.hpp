#ifndef POISSONPROBLEM_HPP
#define POISSONPROBLEM_HPP
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
    const int nodofV=1, nodofE=1, nodofF=1, nodofV=1;
    
  public:
    poissonProblem(int DIM);

    //The number of degrees
    //of freedom on each entity
    //type
    int getNODOFVert();
    int getNODOFEdge();
    int getNODOFFace();
    int getNODOFVolm();

    void calcQpResidual(vector<Vector> *ElmVecs);

    void calcQpJacobian(vector<Matrix> *ElmJacs);
};

int poissonProblem::getNODOFVert(){return nodofV;};

int poissonProblem::getNODOFEdge(){return nodofE;};

int poissonProblem::getNODOFFace(){return nodofF;};

int poissonProblem::getNODOFVolm(){return nodofV;};


void poissonProblem::calcQpResidual(vector<Vector> *ElmVecs){

};


void poissonProblem::calcQpJacobian(vector<Matrix> *ElmJacs){



};
#endif