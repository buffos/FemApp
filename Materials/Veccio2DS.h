// Veccio2DS.h: interface for the Veccio2DS class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECCIO2DS_H__65D30646_0DD5_11D3_B4AD_8CFA05929556__INCLUDED_)
#define AFX_VECCIO2DS_H__65D30646_0DD5_11D3_B4AD_8CFA05929556__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"

class Veccio2DS : public Material  
{
public:
	double StressStrainCurve(double strain,double bita);
	ColumnVector Stress(Element* elem,int gausspointnumber);
	Veccio2DS(double fc,double ec,bool confined);
	virtual ~Veccio2DS();
	double Poisson(double e2);//e2 ειναι η μεγιστη θλιπτική παραμόρφωση με  προσημο (-)
	double E(double strain,double bita);
	Matrix C(Element* elem,int gausspointnumber);

};


#endif // !defined(AFX_VECCIO2DS_H__65D30646_0DD5_11D3_B4AD_8CFA05929556__INCLUDED_)
