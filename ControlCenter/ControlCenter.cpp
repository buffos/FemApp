// ControlCenter.cpp: implementation of the ControlCenter class.
//
//////////////////////////////////////////////////////////////////////

#include "ControlCenter.h"
#include "FEM.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


ControlCenter::ControlCenter()
{

}

ControlCenter::~ControlCenter()
{

}


void  ControlCenter::SetCurrentDisplacementZero()
{
	int ActiveSystemDofs=fem->GetNumberOfActiveSystemDofs();
	int nfixes;
	if(fem->ConstrainCenter->GetConstrainMethod()==PENALTY_METHOD)
	{
		nfixes=0;	
	}
	else
	{
		nfixes=fem->countfixes()+fem->countkinimatic(fem->CurrentStage);	
	}

	ColumnVector Zero(ActiveSystemDofs+nfixes);Zero=0;
	LastDisplacement=Zero;
	CurrentDisplacement=Zero;
	DisplacementPreviousStage=Zero;
	DisplacementPreviousIteration=Zero;
}
