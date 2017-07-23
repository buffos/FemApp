// PlaneStress.cpp: implementation of the PlaneStress class.
//
//////////////////////////////////////////////////////////////////////

#include "PlaneStress.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlaneStress::PlaneStress(double YoungsMeter,double poisson)
{
	this->SetDescription("Plane Stress");
	this->SetNumberOfParameters(2);
	this->SetParameterSize(2);
	this->SetParameter(1,YoungsMeter);
	this->SetParameter(2,poisson);
}

PlaneStress::~PlaneStress()
{

}

Matrix PlaneStress::C(Element *elem, int gausspointnumber)
{
	Matrix C(3,3);
	double E=this->GetParameter(1);
	double v=this->GetParameter(2);
	C=0;
	C(1,1)=1;
	C(1,2)=v;
	C(2,1)=v;
	C(2,2)=1;
	C(3,3)=0.5*(1-v);
	C=( E/(1-v*v) ) * C;
	return C;
}

ColumnVector PlaneStress::Stress(Element* elem , int gausspointnumber)
{
	Matrix C(3,3);
	double E=this->GetParameter(1);
	double v=this->GetParameter(2);
	ColumnVector strain=elem->GetStrain(gausspointnumber);
	C=0;
	C(1,1)=1;
	C(1,2)=v;
	C(2,1)=v;
	C(2,2)=1;
	C(3,3)=0.5*(1-v);
	C=( E/(1-v*v) ) * C;
	return C*strain;
}
