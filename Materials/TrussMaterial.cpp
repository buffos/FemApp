// TrussMaterial.cpp: implementation of the TrussMaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "TrussMaterial.h"
#include <stdio.h>
#include <stdarg.h>
#include <iostream.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TrussMaterial::TrussMaterial(double Eo)
{
	this->SetParameterSize(7);
	this->SetDescription("Yliko Gia Truss");
	this->SetNumberOfParameters(7);
	this->SetParameter(1,Eo);
	this->SetParameter(2,100);
	this->SetParameter(3,Eo);
	this->SetParameter(4,200);
	this->SetParameter(5,Eo);
	this->SetParameter(6,300);
	this->SetParameter(7,Eo);

}

TrussMaterial::TrussMaterial(double Eo,double e1,double E1,double e2,double E2, double e3, double E3)
{
	this->SetParameterSize(7);
	this->SetDescription("Yliko Gia Truss");
	this->SetNumberOfParameters(7);
	this->SetParameter(1,Eo);
	this->SetParameter(2,e1);
	this->SetParameter(3,E1);
	this->SetParameter(4,e2);
	this->SetParameter(5,E2);
	this->SetParameter(6,e3);
	this->SetParameter(7,E3);
}

TrussMaterial::~TrussMaterial()
{

}

void TrussMaterial::info()
{
	cout<<"Yliko "<<"\t"<<this->GetDescription()<<"\t"<<"Metro Elastikotitas  E ="<<this->GetParameter(1)<<"\n";

}

Matrix TrussMaterial::C(Element *elem, int gausspointnumber)
{
	Matrix t(1,1);
	double strain= elem->GetStrain(gausspointnumber)(1);
	if( fabs(strain)<GetParameter(2) )
	{
		t(1,1)=this->GetParameter(1);
	}
	else if (GetParameter(2)<=fabs(strain) && fabs(strain)<GetParameter(4) )
	{
		t(1,1)=this->GetParameter(3);
	}
	else if ( GetParameter(4)<=fabs(strain) && fabs(strain)<GetParameter(6) )
	{
		t(1,1)=this->GetParameter(5);
	}
	else if ( GetParameter(6)<=fabs(strain) )
	{
		t(1,1)=this->GetParameter(7);
	}
	return t;

}

ColumnVector TrussMaterial::Stress(Element *elem, int gausspointnumber)
{
	ColumnVector stress(1);
	double strain= elem->GetStrain(gausspointnumber)(1);
	stress=this->GetParameter(1)*strain;
	if( fabs(strain)<GetParameter(2) )  // å < å1
	{
		stress(1)=this->GetParameter(1)*strain;
	}
	else if (GetParameter(2)<=fabs(strain) && fabs(strain)<GetParameter(4) )  // å1 <å < å2
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) );
	}
	else if ( GetParameter(4)<=fabs(strain) && fabs(strain)<GetParameter(6))  // å2 < å <å3
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( this->GetParameter(4)-this->GetParameter(2) ) + sign(strain)*this->GetParameter(5)*( fabs(strain)-this->GetParameter(4) );
	}
	else if ( GetParameter(6)<=fabs(strain))  // å3 < å
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( this->GetParameter(4)-this->GetParameter(2)  ) 
			+ sign(strain)*this->GetParameter(5)*( this->GetParameter(6)-this->GetParameter(4) ) + sign(strain)*this->GetParameter(7)*( fabs(strain)-this->GetParameter(6) );
	}
	return stress;
}
