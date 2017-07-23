// Analysis.h: interface for the Analysis class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ANALYSIS_H__623AA634_8E01_4B22_9EA5_4A5C20A9B25A__INCLUDED_)
#define AFX_ANALYSIS_H__623AA634_8E01_4B22_9EA5_4A5C20A9B25A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#define _MATRICES_
#define _CONSTANTS_
#define _CONSTRAINS_
#include "FemIncludes.h"

class FEM;
class Analysis  
{
public:
	Analysis();
	Analysis(FEM* pfem){fem=pfem;}
	virtual ~Analysis();

	void virtual Solve(){;}

	FEM* fem;

};

#endif // !defined(AFX_ANALYSIS_H__623AA634_8E01_4B22_9EA5_4A5C20A9B25A__INCLUDED_)
