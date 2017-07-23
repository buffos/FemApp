// ParametersCenter.h: interface for the ParametersCenter class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARAMETERSCENTER_H__713FB1DA_501C_4390_A695_840E6C603518__INCLUDED_)
#define AFX_PARAMETERSCENTER_H__713FB1DA_501C_4390_A695_840E6C603518__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define _MATRICES_
#define _CONSTANTS_
#include "FemIncludes.h"


class FEM;
class ParametersCenter  
{
public:
	ParametersCenter();
	ParametersCenter(FEM* pfem){fem=pfem;}
	virtual ~ParametersCenter();

private: 
	FEM* fem;


};

#endif // !defined(AFX_PARAMETERSCENTER_H__713FB1DA_501C_4390_A695_840E6C603518__INCLUDED_)
