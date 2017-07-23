// ElasticSolver.h: interface for the ElasticSolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ELASTICSOLVER_H__43EC0DD6_0636_45FA_94C3_F2C86DB54035__INCLUDED_)
#define AFX_ELASTICSOLVER_H__43EC0DD6_0636_45FA_94C3_F2C86DB54035__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ALGORITHMS\Analysis.h"

class ElasticSolver : public Analysis  
{
public:
	ElasticSolver();
	ElasticSolver(FEM* pfem){fem = pfem ;}
	virtual ~ElasticSolver();
	void Solve();

};

#endif // !defined(AFX_ELASTICSOLVER_H__43EC0DD6_0636_45FA_94C3_F2C86DB54035__INCLUDED_)
