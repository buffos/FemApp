// Veccio2DT.cpp: implementation of the Veccio2D class.
//
//////////////////////////////////////////////////////////////////////

#include "Veccio2DT.h"
#include "iostream.h" // for testing
#include "Newmatio.h" // for testing
#include <fstream.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Veccio2DT::Veccio2DT(double fc,double ec,bool confined)
{
	this->SetParameterSize(3);
	this->SetDescription("Veccio Collins Curve 1986");
	this->SetNumberOfParameters(3);
	this->SetParameter(1,fc);//Ταση διαρροής
	this->SetParameter(2,ec);//Παραμόρφωση διαρροής
	this->SetParameter(3,(double)confined);

}

Veccio2DT::~Veccio2DT()
{

}

Matrix Veccio2DT::C(Element *elem, int gausspointnumber)
{
	Matrix mat(3,3);
	Matrix T(3,3);//Πίνακας αλλαγής συντεταγμένων
	mat=0;
	ColumnVector strain=elem->GetStrain(gausspointnumber);
	strain=PrincipalVectorStrain(strain);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
	double bita;
	double v=Poisson( Min(strain(1),strain(2)) );
	double ec=GetParameter(2);
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
	mat(1,1)=E(strain(1),bita);
	mat(2,2)=E(strain(2),bita);
//	mat(3,3)=mat(1,1)*mat(2,2)*(1-v*v)/( mat(1,1)*(1+v)+mat(2,2)*(1+v) );
	mat(3,3)=E(0,1)*(1-v)/2;
	mat(1,2)=sqrt( v*v*fabs( mat(1,1)*mat(2,2)) );//*sign(mat(1,1)*mat(2,2));
	mat(2,1)=mat(1,2);
	mat=(1/(1-v*v))*mat;
	T=TransposeStrain(strain(3));
	mat=T.t()*mat*T;
//			fstream Stiffness2;
//			Stiffness2.open("Test.txt",ios::app);//Εξωτερικά φορτία
//			if(elem->GetIndex()==1)
//			{
//				cout<<"";
//			}
//			Stiffness2<<"Arithmos Stoixeiou "<<elem->GetIndex() << "\n"<<mat<<"\n";
//			flush(Stiffness2);
//			Stiffness2.close();
	return mat;
}

double Veccio2DT::E(double strain, double bita)
{
	double young;
	double fc=GetParameter(1);
	double ec=GetParameter(2);
	double fcr=0.1*fabs(fc);
	if(strain<0 && strain>bita*ec) //Θλιπτική τάση
	{
		young=2*fc*(bita*ec-strain)/(bita*ec*ec);
	}
	if(strain<=bita*ec && GetParameter(3)==1.0) //Θλιπτική τάση
	{
		young=0.00000001;
	}
	if(strain<=bita*ec && strain>2*ec && GetParameter(3)==0.0) //Θλιπτική τάση
	{
		young=-bita*fc*2*(strain-bita*ec)/pow(2*ec-bita*ec,2);
	}
	if(strain<=2*ec && GetParameter(3)==0.0) //Θλιπτική τάση
	{
		young=0.00000001;
	}
	if(strain>=0 && strain<0.05*fabs(ec) ) //Εφελκυστική τάση
	{
		young=2*fabs(fc/ec);
	}
	if(strain>=0.05*fabs(ec) ) //Εφελκυστική τάση
	{
			young=-fcr*100/( sqrt(200*strain)*pow( 1+sqrt(200*strain),2) );
	}
	return young;
}

double Veccio2DT::Poisson(double e2)
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

ColumnVector Veccio2DT::Stress(Element *elem, int gausspointnumber)
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

double Veccio2DT::StressStrainCurve(double strain,double bita)
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
