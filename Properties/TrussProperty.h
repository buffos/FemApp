// TrussProperty.h: interface for the TrussProperty class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRUSSPROPERTY_H__4225D8E3_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
#define AFX_TRUSSPROPERTY_H__4225D8E3_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Property.h"

class TrussProperty : public Property  
{
public:
	TrussProperty(double CrossSection);
	virtual ~TrussProperty();

};

#endif // !defined(AFX_TRUSSPROPERTY_H__4225D8E3_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
