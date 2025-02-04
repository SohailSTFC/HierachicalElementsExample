#ifndef ABSTRACTPROBLEM_HPP
#define ABSTRACTPROBLEM_HPP
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
template<Data>
class abstractProblem{
  private:


  public:
    //The number of degrees
    //of freedom on each entity
    //type
    int getNODOFVert();
    int getNODOFEdge();
    int getNODOFFace();
    int getNODOFVolm();

    virtual void calcQpResidual(vector<Vector> *ElmVecs) = 0;

    virtual void calcQpJacobian(vector<Matrix> *ElmJacs) = 0;
};
#endif