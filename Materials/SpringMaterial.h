// SpringMaterial.h: interface for the SpringMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SPRINGMATERIAL_H__DEF3ADA1_18E7_11D3_B4AD_C705777D715F__INCLUDED_)
#define AFX_SPRINGMATERIAL_H__DEF3ADA1_18E7_11D3_B4AD_C705777D715F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"
//Το υλικό αυτό χρησιμοποιείται μόνο σε tangent stifness μεθόδους
//Για secant stiffness μεθόδους πρέπει να γίνεται χρήση του υλίκου SpringMaterialS

class SpringMaterial : public Material  
{
public:
	SpringMaterial(double Ko);
	ColumnVector Stress(Element* elem , int gausspointnumber);
	void info();
	Matrix C(Element* elem,int gausspointnumber);
	SpringMaterial(double Ko,double u1,double K1,double u2,double K2);
	virtual ~SpringMaterial();

};

#endif // !defined(AFX_SPRINGMATERIAL_H__DEF3ADA1_18E7_11D3_B4AD_C705777D715F__INCLUDED_)
