// PlaneElementProperty.h: interface for the PlaneElementProperty class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PLANEELEMENTPROPERTY_H__7C22AC81_FD94_11D2_B4AD_CE37F915E471__INCLUDED_)
#define AFX_PLANEELEMENTPROPERTY_H__7C22AC81_FD94_11D2_B4AD_CE37F915E471__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Property.h"

class PlaneElementProperty : public Property  
{
public:
	PlaneElementProperty(double thickness);
	virtual ~PlaneElementProperty();

};

#endif // !defined(AFX_PLANEELEMENTPROPERTY_H__7C22AC81_FD94_11D2_B4AD_CE37F915E471__INCLUDED_)
