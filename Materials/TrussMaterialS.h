// TrussMaterialS.h: interface for the TrussMaterialS class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRUSSMATERIALS_H__4225D8E2_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
#define AFX_TRUSSMATERIALS_H__4225D8E2_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"

class TrussMaterialS : public Material  
{
public:
	ColumnVector Stress(Element* elem,int gausspointnumber);
	Matrix C(Element* elem,int gausspointnumber);
	void info();
	TrussMaterialS(double Eo);
	TrussMaterialS(double Eo,double e1,double E1,double e2,double E2, double e3,double E3);	virtual ~TrussMaterialS();

};

#endif // !defined(AFX_TRUSSMATERIALS_H__4225D8E2_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
