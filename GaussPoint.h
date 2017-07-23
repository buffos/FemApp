// GaussPoint.h: interface for the GaussPoint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GAUSSPOINT_H__B7B84864_6FE4_11D2_B4AC_A18C90874F72__INCLUDED_)
#define AFX_GAUSSPOINT_H__B7B84864_6FE4_11D2_B4AC_A18C90874F72__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define _MATRICES_
#define _DEFINITIONS_
#define _EXTRAFUNCTIONS_

#include "FemIncludes.h"
#include "Element.h"


class FEM;
class GaussPoint  
{
public:
	void SetStress(ColumnVector Stress);
	void SetStrain(ColumnVector Strain);
	ColumnVector GetStress();
	ColumnVector GetStrain();
	void Initialise(int size);
	GaussPoint();
	virtual ~GaussPoint();

private:
	ColumnVector strain;
	ColumnVector stress;

};

#endif // !defined(AFX_GAUSSPOINT_H__B7B84864_6FE4_11D2_B4AC_A18C90874F72__INCLUDED_)
