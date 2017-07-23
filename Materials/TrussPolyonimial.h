// TrussPolyonimial.h: interface for the TrussPolyonimial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRUSSPOLYONIMIAL_H__C30CA470_E225_4CA4_8656_85C7B18DDA55__INCLUDED_)
#define AFX_TRUSSPOLYONIMIAL_H__C30CA470_E225_4CA4_8656_85C7B18DDA55__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"

class TrussPolyonimial : public Material  
{
public:
	ColumnVector Stress(Element* elem,int gausspointnumber);
	Matrix C(Element* elem,int gausspointnumber);
	void info();
	TrussPolyonimial(double a0,double a1,double a2,double a3,double a4,double a5);
	virtual ~TrussPolyonimial();

};

#endif // !defined(AFX_TRUSSPOLYONIMIAL_H__C30CA470_E225_4CA4_8656_85C7B18DDA55__INCLUDED_)
