// Constrains.h: interface for the Constrains class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CONSTRAINS_H__C49CCC6E_3012_414B_A225_57A310787798__INCLUDED_)
#define AFX_CONSTRAINS_H__C49CCC6E_3012_414B_A225_57A310787798__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define _MATRICES_
#define _CONSTANTS_
#include "FemIncludes.h"


class FEM;
class Constrains  
{
public:
	Constrains(FEM* fempointer){fem=fempointer;}
	Constrains(){;}
	virtual ~Constrains();

	void ConstrainMatrix(StiffnessMatrix& K,int stage);
	void PenaltyConstrainMatrix(StiffnessMatrix& K,int stage);
	void LagrangeConstrainMatrix(StiffnessMatrix& K,int stage);
	void AugmentedConstrainMatrix(StiffnessMatrix& K,int stage);

//Constrain Matrix with Non Force Creating equations
	void BConstrainMatrix(StiffnessMatrix& K,int stage);
	void BPenaltyConstrainMatrix(StiffnessMatrix& K,int stage);
	void BLagrangeConstrainMatrix(StiffnessMatrix& K,int stage);
	void BAugmentedConstrainMatrix(StiffnessMatrix& K,int stage);

	void ConstrainVector(ColumnVector& V,int stage);
	void PenaltyConstrainVector(ColumnVector& V,int stage);
	void LagrangeConstrainVector(ColumnVector& V,int stage);
	void AugmentedConstrainVector(ColumnVector& V,int stage);

	void BConstrainVector(ColumnVector& V,int stage,double incrfactor);
	void BPenaltyConstrainVector(ColumnVector& V,int stage,double incrfactor);
	void BLagrangeConstrainVector(ColumnVector& V,int stage,double incrfactor);
	void BAugmentedConstrainVector(ColumnVector& V,int stage,double incrfactor);

	
	ColumnVector GetConstrainLoadVector(int stage,double incrfactor);
	ColumnVector GetPenaltyConstrainLoadVector(int stage,double incrfactor);
	ColumnVector GetLagrangeConstrainLoadVector(int stage,double incrfactor);
	ColumnVector GetAugmentedConstrainLoadVector(int stage,double incrfactor);



	int GetConstrainMethod(){return ConstrainMethod;}
	void SetConstrainMethod(int method){ConstrainMethod=method;}
	 // This is matrix B in B*U=V equation witch defines the cinematic condition
	Matrix CreateConstrainsMatrix(int stage);
	ColumnVector CreateConstrainsVector(int stage);
	Matrix CreateBConstrainsMatrix(int stage);
private:
	FEM* fem;
	int ConstrainMethod;

};

#endif // !defined(AFX_CONSTRAINS_H__C49CCC6E_3012_414B_A225_57A310787798__INCLUDED_)
