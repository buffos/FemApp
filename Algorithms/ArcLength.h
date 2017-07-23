// ArcLength.h: interface for the ArcLength class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARCLENGTH_H__72C582BB_C670_409D_A850_3A662A8C857F__INCLUDED_)
#define AFX_ARCLENGTH_H__72C582BB_C670_409D_A850_3A662A8C857F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


//#include "Fem.h"
#include "ALGORITHMS\Analysis.h"

class ControlCenter;
class ParametersCenter;
class ArcLength : public Analysis  
{
public:
	ArcLength();
	ArcLength(FEM* pfem);
	virtual ~ArcLength();
	void Solve();

	/*----------------------------------------------------------\
	|															|
	|				    	Auxilliary Procudures				|
	|															|
	\----------------------------------------------------------*/

void ArcLengthMethod();
void ldl (Matrix A,LowerTriangularMatrix &L,DiagonalMatrix &D,int &neg, int n); 
void ScaleUp();
void Accelerate(){;}
int IterateToEquilibrium();
void ArcLengthRoutine(ColumnVector &DT,ColumnVector &DELBAR,int &illfail,int iterations);
void SolveQuadraticEquation(double a,double b,double c, double &R1,double &R2,double &Rlinear,int &failcondition);
void PrepareNextIncr();



//Variables
ParametersCenter* pc;
ControlCenter* cc;
double incrloadfactor;double totalloadfactor;double oldtotalloadfactor;
double Stifi;//Initial Stiffness Parameter
double CStif;//Current Stiffness Parameter
double StrainEnergy;
double bet;
double dl;//this is <DT,DT>
double arc;
int neg;//number of negative pivot points of K in an LDL analysis
int iterations; //iterations made by iterateto equilibrium function. Not Load Steps.Load Steps = incr
int ActiveSystemDofs;

StiffnessMatrix K,Kinv;
ColumnVector LPS,LCS;//LPS = Load Previous Stage   LCS==Load Current Stage
ColumnVector DT;//DT = Total Displacement= Kinv*LCS
ColumnVector Residual;//DT = Total Displacement= Kinv*LCS

};

#endif // !defined(AFX_ARCLENGTH_H__72C582BB_C670_409D_A850_3A662A8C857F__INCLUDED_)
