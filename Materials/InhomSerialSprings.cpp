// InhomSerialSprings.cpp: implementation of the InhomSerialSprings class.
//
//////////////////////////////////////////////////////////////////////

#include "InhomSerialSprings.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

InhomSerialSprings::InhomSerialSprings(double Eo,int NElat,double pincr,double Length)
{
/*-------------------------------------------------------------------\
| We have NElat number of springs in serial order. Each has a        |
| stiffness coefficient kn=En*An/Ln where An=area of beam ln L/NElat |
| and En = Eo + Bita*(1/n^2). Bita has to be calculated so the       |
| increase En-E1/E1 is exacly pincr                                  |
\-------------------------------------------------------------------*/
	this->SetDescription("Representation of series of springs");
	this->SetNumberOfParameters(5);
	this->SetParameterSize(5);
	double poisson=0;
//	l=poisson*YoungsMeter/((1+poisson)*(1-2*poisson));
//	m=YoungsMeter/(2+2*poisson);
	this->SetParameter(1,Eo);
	this->SetParameter(2,poisson);
	this->SetParameter(3,NElat);//Number of Springs
	this->SetParameter(4,pincr);//Increase percent in E value
	this->SetParameter(5,Length);//Length of beam


}

InhomSerialSprings::~InhomSerialSprings()
{

}


Matrix InhomSerialSprings::C(Element *elem, int gausspointnumber)
{
	Matrix El(3,3);

	double x=0; 

	x = elem->GetGlobalCoords(gausspointnumber,1);

	double E = GetParameter(1);
	double N= GetParameter(3);// NUmber of Springs
	double L= GetParameter(5); // Beam Length
	double v = 0.0L;
	double pinc = this->GetParameter(4);

	double step= (L/N);
	int currentspring=x/step+1;
	double offset = (x/step+1-currentspring);


	E= E*(1+offset*pinc);// This is En
	El=0.0L;
	El(1,1)=1;	El(1,2)=v;	El(2,1)=v;	El(2,2)=1;	El(3,3)=0.5L*(1.0L-v);
	El=(E/(1.0L-v*v))*El;
	return El;
}

ColumnVector InhomSerialSprings::Stress(Element* elem , int gausspointnumber)
{
	ColumnVector strain=elem->GetStrain(gausspointnumber);
	Matrix El(3,3);


	double x=0; 

	x = elem->GetGlobalCoords(gausspointnumber,1);

	double E = GetParameter(1);
	double N= GetParameter(3);// NUmber of Springs
	double L= GetParameter(5); // Beam Length
	double v = 0.0L;
	double pinc = this->GetParameter(4);

	double step= (L/N);
	int currentspring=x/step+1;
	double offset = (x/step+1-currentspring);

	E= E*(1+offset*pinc);// This is En
	El=0.0L;
	El(1,1)=1;	El(1,2)=v;	El(2,1)=v;	El(2,2)=1;	El(3,3)=0.5L*(1.0L-v);
	El=(E/(1.0L-v*v))*El;
	return El*strain;
}

