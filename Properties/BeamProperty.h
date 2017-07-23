// BeamProperty.h: interface for the BeamProperty class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BEAMPROPERTY_H__0F9D1002_FE61_11D2_B4AD_E262055F05D1__INCLUDED_)
#define AFX_BEAMPROPERTY_H__0F9D1002_FE61_11D2_B4AD_E262055F05D1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Property.h"

class BeamProperty : public Property  
{
public:
	BeamProperty(double CrossSection,double Inertia);
	virtual ~BeamProperty();

};

#endif // !defined(AFX_BEAMPROPERTY_H__0F9D1002_FE61_11D2_B4AD_E262055F05D1__INCLUDED_)
