// Veccio2DS.cpp: implementation of the Veccio2DS class.
//
//////////////////////////////////////////////////////////////////////

#include "Veccio2DS.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Veccio2DS::Veccio2DS(double fc,double ec,bool confined)
{
	this->SetParameterSize(3);
	this->SetDescription("Veccio Collins Curve 1986");
	this->SetNumberOfParameters(3);
	this->SetParameter(1,fc);//Ταση διαρροής
	this->SetParameter(2,ec);//Παραμόρφωση διαρροής
	this->SetParameter(3,(double)confined);
}

Veccio2DS::~Veccio2DS()
{

}

Matrix Veccio2DS::C(Element *elem, int gausspointnumber)
{
	Matrix mat(3,3);
	Matrix T(3,3);//Πίνακας αλλαγής συντεταγμένων
	mat=0;
	ColumnVector strain=elem->GetStrain(gausspointnumber);
	strain=PrincipalVectorStrain(strain);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
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

double Veccio2DS::Poisson(double e2)
{
	double ec=GetParameter(2);
	double v;
	if (e2/ec<=0.85 && e2<=0)
	{
		v=0.23;
	}
	if (e2/ec>0.85)
	{
		v=0.23+(0.4-.23)*(fabs(e2/ec)-0.85)/(1-0.85);
	}
	if (e2>0)
	{
		v=0;
	}
	return v;
}

double Veccio2DS::E(double strain, double bita)
{
	double young;
	double fc=GetParameter(1);
	double ec=GetParameter(2);
	double fcr=0.1*fabs(fc);
	if(strain<0 && strain>bita*ec) //Θλιπτική τάση
	{
		young= fc*bita*( 2/(bita*ec)-strain/(bita*bita*ec*ec) );
	}
	if(strain<=bita*ec && GetParameter(3)==1.0) //Θλιπτική τάση
	{
		young=fc*bita/strain;
	}
	if(strain<=bita*ec && strain>2*ec && GetParameter(3)==0.0) //Θλιπτική τάση
	{
		young=bita*fc*( 1- pow(strain-bita*ec,2)/pow(2*ec-bita*ec,2) )/strain;
	}
	if(strain<=2*ec && GetParameter(3)==0.0) //Θλιπτική τάση
	{
		young=0.00000000001;
	}
	if(strain>=0 && strain<0.05*fabs(ec) ) //Θλιπτική τάση
	{
		young=2*fabs(fc/ec);
	}
	if(strain>=0.05*fabs(ec) ) //Θλιπτική τάση
	{
		young=fcr/( (1+sqrt(200*strain)) * strain );
	}
	return young;
}

ColumnVector Veccio2DS::Stress(Element *elem, int gausspointnumber)
{
	ColumnVector strain=elem->GetStrain(gausspointnumber);
	strain=PrincipalVectorStrain(strain);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
	ColumnVector stress(3);
	double bita;
	if ( (strain(1)==0 && strain(2)==0 && strain(3)==0) || (strain(1)>=0 && strain(2)>=0))
	{ 
		bita=1;
	}
	else
	{
		bita=Min( 1/( 0.85- ( Max(strain(1),strain(2))*0.27 /Min(strain(1),strain(2)) )) , 1);
		if (bita>1)
		{
			bita=1;
		}
	}
	stress(1)=StressStrainCurve(strain(1),bita);
	stress(2)=StressStrainCurve(strain(2),bita);
	stress(3)=0;
	Matrix T(3,3);
	T=TransposeStrain(strain(3));
	stress=T.t()*stress;
	return stress;

}

double Veccio2DS::StressStrainCurve(double strain,double bita)
{
	double ec=GetParameter(2);
	double fc=GetParameter(1);
	double fcr=0.1*fabs(fc);
	double stress;
	double ep=ec*bita;
	double fp=fc*bita;
	if(strain<0 && strain>bita*ec) //Θλιπτική τάση
	{
		stress=fp*(2*strain/ep-strain*strain/(ep*ep));
	}
	if(strain<=bita*ec && GetParameter(3)==1.0) //Θλιπτική τάση
	{
		stress=fp;
	}
	if(strain<=bita*ec && strain>2*ec && GetParameter(3)==0.0) //Θλιπτική τάση
	{
		stress=fp*(1- pow(strain-ep,2)/pow(2*ec-ep,2) );
	}
	if(strain<=2*ec && GetParameter(3)==0.0) //Θλιπτική τάση
	{
		stress = 0;
	}
	if(strain>=0 && strain<0.05*fabs(ec) ) //Θλιπτική τάση
	{
		stress=2*fabs(fc/ec)*strain;
	}
	if(strain>=0.05*fabs(ec) ) //Θλιπτική τάση
	{
		stress=fcr/(1+ sqrt(200*strain) );
	}
	return stress;

}