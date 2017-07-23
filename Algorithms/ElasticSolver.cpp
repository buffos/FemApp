// ElasticSolver.cpp: implementation of the ElasticSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "ElasticSolver.h"
#include "FEM.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ElasticSolver::ElasticSolver()
{

}

ElasticSolver::~ElasticSolver()
{

}


void ElasticSolver::Solve()
{
	int n = fem->GetNumberOfActiveSystemDofs();
	StiffnessMatrix K;									// Stiffness Matrix Vector
	ColumnVector L;										//Load Vector
	ColumnVector U;										//Displacement Vector
	K=fem->AssembleStiffnessMatrix();
	L=fem->AssembleLoads()+ fem->AssembleStageLoads(1);


	fem->CCenter->SetStiffnessMatrix(K);				//Storing the stiffness Matrix is not needed for Elastic Solve
	fem->ConstrainCenter->ConstrainMatrix(K,1);			//Impose Boundary Conditions
	fem->ConstrainCenter->ConstrainVector(L,1);

	if (fem->SolverType == 1)							//Solving
	{
		U=GaussSolver(K,L);
	}
	fem->UpdateNodes(U.Rows(1,n));						// Update Only the active dofs Ignoring Lagrange Multipliers
	fem->CCenter->SetCurrentDisplacement(U.Rows(1,n));	//Storing the Displacement vector is not needed for Elastic Solve

	/*----------------------------------------------------------\
	|															|
	|				    	Results Storage						|
	|															|
	\----------------------------------------------------------*/

	fstream loads;loads.open(fem->LoadVectorFileOut,ios::app);		//Open filename for External Loads
		printvector(loads,L,1,n);
		loads.close();
	fstream displ;displ.open(fem->DisplacementsFileOut,ios::app);	//Open filename for Displacments
		printvector(displ,U,1,n);
		displ.close();
	fem->ReactionToFile(fem->ReactionFileOut);				//Reactions
	fem->InternalForceToFile(fem->ForceFileOut);			//Internal Forces
	fem->StressAndStainToFile(fem->ElementFileOut,true);	//Principal Strain and Stresses
	fem->StressAndStainToFile("Stresses.txt",false);		//Strain and Stresses
}