// TrussMaterialS.cpp: implementation of the TrussMaterialS class.
//
//////////////////////////////////////////////////////////////////////

#include "TrussMaterialS.h"
#include <stdio.h>
#include <stdarg.h>
#include <iostream.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TrussMaterialS::TrussMaterialS(double Eo)
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

TrussMaterialS::TrussMaterialS(double Eo,double e1,double E1,double e2,double E2,double e3,double E3)
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

TrussMaterialS::~TrussMaterialS()
{

}

void TrussMaterialS::info()
{
	cout<<"Yliko "<<"\t"<<this->GetDescription()<<"\t"<<"Metro Elastikotitas  E ="<<this->GetParameter(1)<<"\n";

}

Matrix TrussMaterialS::C(Element *elem, int gausspointnumber)
{
	Matrix t(1,1);//Υπολογίζω σε κάθε περίπτωση το τέμνων μέτρο δυκαμψίας
	double strain= elem->GetStrain(gausspointnumber)(1);
	//Είναι γενικα k=σ/ε
	if( fabs(strain)<GetParameter(2) )
	{
		t(1,1)=this->GetParameter(1);
	}
	else if (GetParameter(2)<=fabs(strain) && fabs(strain)<GetParameter(4) )
	{
		t(1,1)=(  this->GetParameter(1)*this->GetParameter(2)+this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) ) )/fabs(strain);
	}
	else if ( GetParameter(4)<=fabs(strain) && fabs(strain)<GetParameter(6) )
	{
		t(1,1)=( this->GetParameter(1)*this->GetParameter(2)+this->GetParameter(3)*( this->GetParameter(4)-this->GetParameter(2) ) + this->GetParameter(5)*( fabs(strain)-this->GetParameter(4) ) )/fabs(strain);
	}
	else if ( GetParameter(6)<=fabs(strain) )
	{
		t(1,1)=this->GetParameter(1)*this->GetParameter(2) + this->GetParameter(3)*( this->GetParameter(4)-this->GetParameter(2)  ) 
			+ this->GetParameter(5)*( this->GetParameter(6)-this->GetParameter(4) ) + this->GetParameter(7)*( fabs(strain)-this->GetParameter(6) );
	}
	return t;
}

ColumnVector TrussMaterialS::Stress(Element *elem, int gausspointnumber)
{
	ColumnVector stress(1);
	double strain= elem->GetStrain(gausspointnumber)(1);
	stress=this->GetParameter(1)*strain;
	if( fabs(strain)<GetParameter(2) )  // ε < ε1
	{
		stress(1)=this->GetParameter(1)*strain;
	}
	else if (GetParameter(2)<=fabs(strain) && fabs(strain)<GetParameter(4) )  // ε1 <ε < ε2
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) );
	}
	else if ( GetParameter(4)<=fabs(strain) && fabs(strain)<GetParameter(6))  // ε2 < ε
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( this->GetParameter(4)-this->GetParameter(2) ) + sign(strain)*this->GetParameter(5)*( fabs(strain)-this->GetParameter(4) );
	}
	else if ( GetParameter(6)<=fabs(strain))  // ε2 < ε
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( this->GetParameter(4)-this->GetParameter(2)  ) 
			+ sign(strain)*this->GetParameter(5)*( this->GetParameter(6)-this->GetParameter(4) ) + sign(strain)*this->GetParameter(7)*( fabs(strain)-this->GetParameter(6) );
	}
	return stress;
}
