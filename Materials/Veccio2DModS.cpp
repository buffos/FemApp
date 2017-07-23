// Veccio2DModS.cpp: implementation of the Veccio2DModS class.
//
//////////////////////////////////////////////////////////////////////

#include "Veccio2DModS.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Veccio2DModS::Veccio2DModS(double fc,double ec,bool confined,int lockingcondition, double ksi, double ita)
{
	this->SetParameterSize(6);
	this->SetDescription("Veccio Collins Curve 1986 Modified");
	this->SetNumberOfParameters(6);
	this->SetParameter(1,fc);//���� ��������
	this->SetParameter(2,ec);//����������� ��������
	this->SetParameter(3,(double)confined);
	this->SetParameter(4,lockingcondition);
	// �� � ���� ��� ����� 0 ���� �� ���� gausspoint ������������ ��� ������� ������ (��������� �����)
	// �� � ���� ��� ����� 1 ���� �������� ��� �� ������ ������ �� ��� �� �������� ����� �������� ���
	// ������������� ��� ������ �� ������������� (ksi,ita)
	// �� � ���� ����� 2 ���� �� ������ ������ ������������� ��� gausspoint �� ��� ���������� ����������� �����������.
	// �� � ���� ����� 3 ���� �� ������ ������ ������������� ��� gausspoint �� ��� ���������� �������� �����������.

	this->SetParameter(5,ksi);//������������ �
	this->SetParameter(6,ita);//������������ �
}

Veccio2DModS::~Veccio2DModS()
{

}

Matrix Veccio2DModS::C(Element *elem, int gausspointnumber)
{
	ColumnVector strain;
	Matrix mat(3,3);
	Matrix T(3,3);//������� ������� �������������
	mat=0;
	strain = this->GiveTheRightStrain(elem,gausspointnumber);
	// ��������� �� �������� ������������� ��� ������ ��� ������� �� ����������� � ������������ �������.
	double bita;
	double v=Poisson( Min(strain(1),strain(2)) );
	if ( (strain(1)==0 && strain(2)==0 && strain(3)==0) || (strain(1)>=0 && strain(2)>=0))
	{ 
		bita=1;
	}
	else
	{
		bita=1/( 0.85- (0.27*Max(strain(1),strain(2)) / Min(strain(1),strain(2)) ));
		if (bita>1)
		{
			bita=1;
		}
	}
	mat(1,1)=E(strain(1),bita);
	mat(2,2)=E(strain(2),bita);
	if(v==0.0)
	{
		mat(3,3)=mat(1,1)*mat(2,2)/( mat(1,1)+mat(2,2) );
	}
	else
	{
		mat(3,3)=mat(1,1)+mat(2,2)-2*sqrt(v*v*mat(1,1)*mat(2,2) )/4;
	}
	mat(1,2)=sqrt( v*v*mat(1,1)*mat(2,2) );
	mat(2,1)=mat(1,2);
	mat=(1/(1-v*v))*mat;
	T=TransposeStrain(strain(3));
	mat=T.t()*mat*T;
	return mat;
	
}

double Veccio2DModS::Poisson(double e2)
{
	double ec=GetParameter(2);
	double v;
	if (e2/ec<=0.85 && e2<=0)
	{
		v=0.23;
	}
	if (e2/ec>0.85)
	{
		v=(0.4-.23)*(fabs(e2)-0.85)/(1-0.85);
	}
	if (e2>0)
	{
		v=0;
	}
	return v;
}

double Veccio2DModS::E(double strain, double bita)
{
	double young;
	double fc=GetParameter(1);
	double ec=GetParameter(2);
	double fcr=0.1*fabs(fc);
	if(strain<0 && strain>bita*ec) //�������� ����
	{
		young= fc*bita*( 2/(bita*ec)-strain/(bita*bita*ec*ec) );
	}
	if(strain<=bita*ec && GetParameter(3)==1.0) //�������� ����
	{
		young=fc*bita/strain;
	}
	if(strain<=bita*ec && strain>2*ec && GetParameter(3)==0.0) //�������� ����
	{
		young=bita*fc*( 1- pow(strain-bita*ec,2)/pow(2*ec-bita*ec,2) )/strain;
	}
	if(strain<=2*ec && GetParameter(3)==0.0) //�������� ����
	{
		young=0.00000000001;
	}
	if(strain>=0 && strain<0.05*fabs(ec) ) //�������� ����
	{
		young=2*fabs(fc/ec);
	}
	if(strain>=0.05*fabs(ec) ) //�������� ����
	{
		young=fcr/( (1+sqrt(200*strain)) * strain );
	}
	return young;
}

ColumnVector Veccio2DModS::Stress(Element *elem, int gausspointnumber)
{
	ColumnVector strain=elem->GetStrain(gausspointnumber);
	ColumnVector stress(3);
	stress = this->C(elem,gausspointnumber)*strain;
	return stress;

}

double Veccio2DModS::StressStrainCurve(double strain,double bita)
{
	double ec=GetParameter(2);
	double fc=GetParameter(1);
	double fcr=0.1*fabs(fc);
	double stress;
	double ep=ec*bita;
	double fp=fc*bita;
	if(strain<0 && strain>bita*ec) //�������� ����
	{
		stress=fp*(2*strain/ep-strain*strain/(ep*ep));
	}
	if(strain<=bita*ec && GetParameter(3)==1.0) //�������� ����
	{
		stress=fp;
	}
	if(strain<=bita*ec && strain>2*ec && GetParameter(3)==0.0) //�������� ����
	{
		stress=fp*(1- pow(strain-ep,2)/pow(2*ec-ep,2) );
	}
	if(strain<=2*ec && GetParameter(3)==0.0) //�������� ����
	{
		stress = 0;
	}
	if(strain>=0 && strain<0.05*fabs(ec) ) //�������� ����
	{
		stress=2*fabs(fc/ec)*strain;
	}
	if(strain>=0.05*fabs(ec) ) //�������� ����
	{
		stress=fcr/(1+ sqrt(200*strain) );
	}
	return stress;

}


ColumnVector Veccio2DModS::GiveTheRightStrain(Element *elem, int gausspointnumber)
{
	ColumnVector strain;
	double ksi,ita;
	int lcase;
	ksi =this->GetParameter(5);
	ita =this->GetParameter(6);
	lcase = this->GetParameter(4);
	int nogp; //number of gausspoints
	nogp = elem->GetNumberOfGaussPoints();

	if (lcase == 0)
	{
		strain=elem->GetStrain(gausspointnumber);
		strain=PrincipalVectorStrain(strain);//����� �� �������� ������ �������������
	}
	else if (lcase == 1)
	{
		strain=elem->GetStrain(ksi,ita);
		strain=PrincipalVectorStrain(strain);//����� �� �������� ������ �������������
	}
	else if (lcase == 2)
	{
		ColumnVector strains[10]; // Mexri 10 gausspoints blepei Gia parap'anv apla prepei na ton allakso ton kodika
		int maxgausspoint;
		double maxvalue = -1e10;
		for(int i=0;i<nogp;i++)
		{
			strains[i] = elem->GetStrain(i);
			strains[i] = PrincipalVectorStrain(strains[i]);//����� �� �������� ������ �������������
		}
		for(i=0;i<nogp;i++)
			for(int j=1;j<=2;j++)
			{
				if (strains[i](j)>maxvalue)
				{ 
					maxgausspoint = i+1;
					maxvalue = strains[i](j);

				}
			}
		strain = strains[maxgausspoint];
	}
	else if (lcase == 3)
	{
		ColumnVector strains[10];
		int mingausspoint;
		double minvalue = 1e10;
		for(int i=0;i<nogp;i++)
		{
			strains[i] = elem->GetStrain(i);
			strains[i] = PrincipalVectorStrain(strains[i]);//����� �� �������� ������ �������������
		}
		for(i=0;i<nogp;i++)
			for(int j=1;j<=2;j++)
			{
				if (strains[i](j)<minvalue)
				{ 
					mingausspoint = i+1;
					minvalue = strains[i](j);
				}
			}
		strain = strains[mingausspoint];
	}
	else if (lcase == 4)
	{
		ColumnVector strains[10];
		double maxdiogosi;
		int maxgausspoint;
		maxdiogosi = -1e10;
		for(int i=0;i<nogp;i++)
		{
			strains[i] = elem->GetStrain(i);
		}
		for(i=0;i<nogp;i++)
		{
			if (strains[i](1)+strains[i](2)> maxdiogosi)
			{ 
					maxgausspoint = i+1;
					maxdiogosi = strains[i](1)+strains[i](2);
			}
		}
		strain = strains[maxgausspoint];
		strain = PrincipalVectorStrain(strain);//����� �� �������� ������ �������������
}
	return strain;
}
