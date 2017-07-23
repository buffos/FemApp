// TrussPolyonimial.cpp: implementation of the TrussPolyonimial class.
//
//////////////////////////////////////////////////////////////////////

#include "TrussPolyonimial.h"
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <iostream.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TrussPolyonimial::TrussPolyonimial(double a0,double a1,double a2,double a3,double a4,double a5)
{
	this->SetParameterSize(6);
	this->SetDescription("Yliko Gia Truss");
	this->SetNumberOfParameters(6);
	this->SetParameter(1,a0);
	this->SetParameter(2,a1);
	this->SetParameter(3,a2);
	this->SetParameter(4,a3);
	this->SetParameter(5,a4);
	this->SetParameter(6,a5);

}

TrussPolyonimial::~TrussPolyonimial()
{

}

Matrix TrussPolyonimial::C(Element *elem, int gausspointnumber)
{
	Matrix t(1,1);
	t(1,1)=0;
	double strain= elem->GetStrain(gausspointnumber)(1);
	if(strain==0)
	{
		t(1,1)=this->GetParameter(2);
	}
	else
	{
		// Paragogos tou polyonymou
		for(int k = 2;k<=this->GetNumberOfParameters();k++)
		{
			t(1,1)+=this->GetParameter(k)*pow(strain,k-2)*(k-1);
		}
//		t(1,1)=t(1,1)/strain;
	}
	return t;

}

void TrussPolyonimial::info()
{
	cout<<"Yliko "<<"\t"<<this->GetDescription()<<"\t"<<"Metro Elastikotitas  E ="<<this->GetParameter(1)<<"\n";

}

ColumnVector TrussPolyonimial::Stress(Element *elem, int gausspointnumber)
{
	Matrix t(1,1);
	t(1,1)=0;
	double strain= elem->GetStrain(gausspointnumber)(1);
	for(int k = 1;k<=this->GetNumberOfParameters();k++)
	{
		t(1,1)+=this->GetParameter(k)*pow(strain,k-1);
	}
	return t;

}
