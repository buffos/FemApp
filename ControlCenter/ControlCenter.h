// ControlCenter.h: interface for the ControlCenter class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CONTROLCENTER_H__4FA3077E_8137_4E18_A00D_6319F550DFDA__INCLUDED_)
#define AFX_CONTROLCENTER_H__4FA3077E_8137_4E18_A00D_6319F550DFDA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#define _MATRICES_
#define _CONSTANTS_
#include "FemIncludes.h"


class FEM;
class ControlCenter  
{
public:
	ControlCenter();
	ControlCenter(FEM* pfem){fem=pfem;}
	virtual ~ControlCenter();

	void  SetDisplacementPreviousStage(ColumnVector d){DisplacementPreviousStage=d;}
	void  SetDisplacementPreviousIteration(ColumnVector d){DisplacementPreviousIteration=d;}
	void  SetCurrentDisplacement(ColumnVector d){CurrentDisplacement=d;}
	void  SetLastDisplacement(ColumnVector d){LastDisplacement=d;}
	void  SetCurrentDisplacementZero();
	void  IncreaseCurrentDisplacement(ColumnVector d){CurrentDisplacement+=d;}

	void  SetExternalLoadFromPreviousStage(ColumnVector d){ExternalLoadFromPreviousStage=d;}
	void  SetExternalLoad(ColumnVector d){ExternalLoad=d;}
	void  SetCurrentExternalLoad(ColumnVector d){ExternalLoadPreviousIteration=ExternalLoadCurrent;ExternalLoadCurrent=d;}
	void  IncreaseCurrentExternalLoad(ColumnVector d){ExternalLoadPreviousIteration=ExternalLoadCurrent;ExternalLoadCurrent+=d;}

	void  SetStiffnessMatrix(StiffnessMatrix Km){K=Km;}

	ColumnVector GetDisplacementPreviousStage(){return DisplacementPreviousStage;}
	ColumnVector GetDisplacementPreviousIteration(){return DisplacementPreviousIteration;}
	ColumnVector GetCurrentDisplacement(){return CurrentDisplacement;}
	ColumnVector GetLastDisplacement(){return LastDisplacement;}
	ColumnVector GetDU(){return (CurrentDisplacement-LastDisplacement);}
	ColumnVector GetDUPreviousIteration(){return (CurrentDisplacement-DisplacementPreviousIteration);}

	ColumnVector GetExternalLoadFromPreviousStage(){return ExternalLoadFromPreviousStage;}
	ColumnVector GetExternalLoad(){return ExternalLoad;}
	ColumnVector GetCurrentExternalLoad(){return ExternalLoadCurrent;}// will always have K*DU
	ColumnVector GetExternalLoadPreviousIteration(){return ExternalLoadPreviousIteration;}

	StiffnessMatrix GetStiffnessMatrix(){return K;}

	ColumnVector CalculateExternalLoad(){IncreaseCurrentExternalLoad(K*(CurrentDisplacement-LastDisplacement));}

private: 
	FEM* fem;

	// Vectors
	ColumnVector DisplacementPreviousStage;        // Displacement at end of previous stage
	ColumnVector DisplacementPreviousIteration;    // Displacement at end of previous iteration
	ColumnVector CurrentDisplacement;              // Displacement write now
	ColumnVector LastDisplacement;                 // Displacement just before this one

	ColumnVector ExternalLoadFromPreviousStage;    // The Current External Load the structure has
	ColumnVector ExternalLoad;                     // Total External Load for the structure for current stage
	ColumnVector ExternalLoadCurrent;				//LCS
	ColumnVector ExternalLoadPreviousIteration;    


	StiffnessMatrix K;                             //StiffnessMatrix of the Structure.

};

#endif // !defined(AFX_CONTROLCENTER_H__4FA3077E_8137_4E18_A00D_6319F550DFDA__INCLUDED_)
