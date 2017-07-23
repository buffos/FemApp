// Veccio2DModS.h: interface for the Veccio2DModS class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_Veccio2DModS_H__65D30646_0DD5_11D3_B4AD_8CFA05929556__INCLUDED_)
#define AFX_Veccio2DModS_H__65D30646_0DD5_11D3_B4AD_8CFA05929556__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"

class Veccio2DModS : public Material  
{
public:
	ColumnVector GiveTheRightStrain(Element* elem,int gausspointnumber);
	double StressStrainCurve(double strain,double bita);
	ColumnVector Stress(Element* elem,int gausspointnumber);
	Veccio2DModS(double fc,double ec,bool confined,int lockingcondition, double ksi, double ita);
	virtual ~Veccio2DModS();
	double Poisson(double e2);//e2 ειναι η μεγιστη θλιπτική παραμόρφωση με  προσημο (-)
	double E(double strain,double bita);
	Matrix C(Element* elem,int gausspointnumber);

};


#endif // !defined(AFX_Veccio2DModS_H__65D30646_0DD5_11D3_B4AD_8CFA05929556__INCLUDED_)
