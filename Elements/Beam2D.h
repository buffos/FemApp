// Beam2D.h: interface for the Beam2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BEAM2D_H__0F9D1001_FE61_11D2_B4AD_E262055F05D1__INCLUDED_)
#define AFX_BEAM2D_H__0F9D1001_FE61_11D2_B4AD_E262055F05D1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Element.h"
#include "FEM.h"

class Beam2D : public Element  
{
public:
	void GeometryToFile(char* filename);
	Matrix CreateStiffnessMatrix();
	Beam2D(FEM* fem,int index,int firstnodenumber,int secondnodenumber);
	Beam2D();
	virtual ~Beam2D();
};

#endif // !defined(AFX_BEAM2D_H__0F9D1001_FE61_11D2_B4AD_E262055F05D1__INCLUDED_)
