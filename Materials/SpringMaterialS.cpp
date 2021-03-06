// SpringMaterialS.cpp: implementation of the SpringMaterialS class.
//
//////////////////////////////////////////////////////////////////////

#include "SpringMaterialS.h"
#include "iostream.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SpringMaterialS::SpringMaterialS(double Ko,double u1,double K1,double u2,double K2)
{
	this->SetParameterSize(5);
	this->SetDescription("Yliko Gia Springs");
	this->SetNumberOfParameters(5);
	this->SetParameter(1,Ko);
	this->SetParameter(2,u1);
	this->SetParameter(3,K1);
	this->SetParameter(4,u2);
	this->SetParameter(5,K2);

}

SpringMaterialS::~SpringMaterialS()
{

}

Matrix SpringMaterialS::C(Element *elem, int gausspointnumber)
{
	Matrix t(1,1);//��������� �� ���� ��������� �� ������ ����� ���������
	double strain= elem->GetStrain(gausspointnumber)(1);
	//����� ������ k=�/�
	if( fabs(strain)<GetParameter(2) )
	{
		t(1,1)=this->GetParameter(1);
	}
	else if (GetParameter(2)<=fabs(strain) && fabs(strain)<GetParameter(4) )
	{
		t(1,1)=(  this->GetParameter(1)*this->GetParameter(2)+this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) ) )/fabs(strain);
	}
	else if ( GetParameter(4)<=fabs(strain) )
	{
		t(1,1)=( this->GetParameter(1)*this->GetParameter(2)+this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) ) + this->GetParameter(5)*( fabs(strain)-this->GetParameter(4) ) )/fabs(strain);
	}
	return t;

}

void SpringMaterialS::info()
{
	cout<<"Yliko "<<"\t"<<this->GetDescription()<<"\n";

}

ColumnVector SpringMaterialS::Stress(Element *elem, int gausspointnumber)
{
	ColumnVector stress(1);
	double strain= elem->GetStrain(gausspointnumber)(1);
	if( fabs(strain)<GetParameter(2) )  // � < �1
	{
		stress(1)=this->GetParameter(1)*strain;
	}
	else if (GetParameter(2)<=fabs(strain) && fabs(strain)<GetParameter(4) )  // �1 <� < �2
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) );
	}
	else if ( GetParameter(4)<=fabs(strain) )  // �2 < �
	{
		stress(1)=sign(strain)*this->GetParameter(1)*this->GetParameter(2)+sign(strain)*this->GetParameter(3)*( fabs(strain)-this->GetParameter(2) ) + sign(strain)*this->GetParameter(5)*( fabs(strain)-this->GetParameter(4) );
	}
	return stress;

}

SpringMaterialS::SpringMaterialS(double Ko)
{
	this->SetParameterSize(5);
	this->SetDescription("Yliko Gia Springs");
	this->SetNumberOfParameters(5);
	this->SetParameter(1,Ko);
	this->SetParameter(2,100);
	this->SetParameter(3,Ko);
	this->SetParameter(4,200);
	this->SetParameter(5,Ko);


}
