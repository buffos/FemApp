// PlaneStress.h: interface for the PlaneStress class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PLANESTRESS_H__AC337041_FD57_11D2_B4AD_ABBCB3C2BC44__INCLUDED_)
#define AFX_PLANESTRESS_H__AC337041_FD57_11D2_B4AD_ABBCB3C2BC44__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"

class PlaneStress : public Material  
{
public:
	Matrix C(Element *elem,int gausspointnumber);
	PlaneStress(double YoungsMeter,double poisson);
	ColumnVector Stress(Element* elem , int gausspointnumber);
	virtual ~PlaneStress();

};

#endif // !defined(AFX_PLANESTRESS_H__AC337041_FD57_11D2_B4AD_ABBCB3C2BC44__INCLUDED_)
