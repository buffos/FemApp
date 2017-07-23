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
	this->SetParameter(1,fc);//Ταση διαρροής
	this->SetParameter(2,ec);//Παραμόρφωση διαρροής
	this->SetParameter(3,(double)confined);
	this->SetParameter(4,lockingcondition);
	// Αν η τιμή του ειναι 0 τότε σε κάθε gausspoint υπολογίζουμε του κύριους άξονες (ορθότροπο υλικό)
	// Αν η τιμή του είναι 1 τότε θεωρούμε οτι οι κύριοι άξόνες σε όλο το στοιχείο είναι σταθεροί και
	// Υπολογίζονται στο σημείο με συντεταγμένες (ksi,ita)
	// Αν η τιμή είναι 2 τότε οι κύριοι αξονες υπολογίζονται στο gausspoint με την μεγαλύτερη εφελκυστική παραμόρφωση.
	// Αν η τιμή είναι 3 τότε οι κύριοι αξονες υπολογίζονται στο gausspoint με την μεγαλύτερη θλιπτική παραμόρφωση.

	this->SetParameter(5,ksi);//συντεταγμενη ξ
	this->SetParameter(6,ita);//συντεταγμενη η
}

Veccio2DModS::~Veccio2DModS()
{

}

Matrix Veccio2DModS::C(Element *elem, int gausspointnumber)
{
	ColumnVector strain;
	Matrix mat(3,3);
	Matrix T(3,3);//Πίνακας αλλαγής συντεταγμένων
	mat=0;
	strain = this->GiveTheRightStrain(elem,gausspointnumber);
	// Υπολογίζω το διάνυσμα παραμορφώσεων στο σημείο που επιθυμώ να υπολογιστεί ο καταστατικός πίνακας.
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
		strain=PrincipalVectorStrain(strain);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
	}
	else if (lcase == 1)
	{
		strain=elem->GetStrain(ksi,ita);
		strain=PrincipalVectorStrain(strain);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
	}
	else if (lcase == 2)
	{
		ColumnVector strains[10]; // Mexri 10 gausspoints blepei Gia parap'anv apla prepei na ton allakso ton kodika
		int maxgausspoint;
		double maxvalue = -1e10;
		for(int i=0;i<nogp;i++)
		{
			strains[i] = elem->GetStrain(i);
			strains[i] = PrincipalVectorStrain(strains[i]);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
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
			strains[i] = PrincipalVectorStrain(strains[i]);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
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
		strain = PrincipalVectorStrain(strain);//Βρήκα το διάνυσμα κύριων παραμορφώσεων
}
	return strain;
}
