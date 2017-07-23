// Veccio2DT.h: interface for the Veccio2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECCIO2D_H__1AD5DA21_09F0_11D3_B4AD_85C7D747173F__INCLUDED_)
#define AFX_VECCIO2D_H__1AD5DA21_09F0_11D3_B4AD_85C7D747173F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"
#include "IsoparametricElement.h"

//Implementation υλικού σκυρόδεμα
//Καταστατικό μοντέλο Veccio-Collins 1986
//Tangent Stiffness Aproach

class Veccio2DT : public Material  
{
public:
	double StressStrainCurve(double strain,double bita);
	ColumnVector Stress(Element* elem,int gausspointnumber);
	double Poisson(double e2);//e2 ειναι η μεγιστη θλιπτική παραμόρφωση με  προσημο (-)
	double E(double strain,double bita);
	Matrix C(Element* elem,int gausspointnumber);
	Veccio2DT(double fc,double ec,bool confined);
	virtual ~Veccio2DT();

};

#endif // !defined(AFX_VECCIO2D_H__1AD5DA21_09F0_11D3_B4AD_85C7D747173F__INCLUDED_)
