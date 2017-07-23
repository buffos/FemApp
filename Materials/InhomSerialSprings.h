// InhomSerialSprings.h: interface for the InhomSerialSprings class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INHOMSERIALSPRINGS_H__7E9B0657_F4B7_4E21_9E93_2AA9E5EB74FC__INCLUDED_)
#define AFX_INHOMSERIALSPRINGS_H__7E9B0657_F4B7_4E21_9E93_2AA9E5EB74FC__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Material.h"

class InhomSerialSprings : public Material  
{
public:
	Matrix C(Element *elem,int gausspointnumber);
	ColumnVector Stress(Element* elem,int gausspointnumber);
	InhomSerialSprings(double Eo,int NElat,double pincr,double Length);
	virtual ~InhomSerialSprings();

};

#endif // !defined(AFX_INHOMSERIALSPRINGS_H__7E9B0657_F4B7_4E21_9E93_2AA9E5EB74FC__INCLUDED_)
