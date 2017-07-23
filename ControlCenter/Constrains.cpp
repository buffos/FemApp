// Constrains.cpp: implementation of the Constrains class.
//
//////////////////////////////////////////////////////////////////////

#include "Constrains.h"
#include "FEM.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Constrains::~Constrains()
{

}


void Constrains::ConstrainMatrix(StiffnessMatrix& K,int stage)
{
	switch(ConstrainMethod)
	{
	case PENALTY_METHOD :
		{
			PenaltyConstrainMatrix(K,stage);
		}
		break;
	case LAGRANGE_MULTIPLIERS :
		{
			LagrangeConstrainMatrix(K,stage);
		}
		break;
	case AUGMENTED_LAGRANGE :
		{
			AugmentedConstrainMatrix(K,stage);
		}
		break;
	}
}


void Constrains::ConstrainVector(ColumnVector& V,int stage)
{
	switch(ConstrainMethod)
	{
	case PENALTY_METHOD :
		{
			PenaltyConstrainVector(V,stage);
		}
		break;
	case LAGRANGE_MULTIPLIERS :
		{
			LagrangeConstrainVector(V,stage);
		}
		break;
	case AUGMENTED_LAGRANGE :
		{
			AugmentedConstrainVector(V,stage);
		}
		break;
	}
}

void Constrains::BConstrainVector(ColumnVector& V,int stage,double incrfactor)
{
	switch(ConstrainMethod)
	{
	case PENALTY_METHOD :
		{
			BPenaltyConstrainVector(V,stage,incrfactor);
		}
		break;
	case LAGRANGE_MULTIPLIERS :
		{
			BLagrangeConstrainVector(V,stage,incrfactor);
		}
		break;
	case AUGMENTED_LAGRANGE :
		{
			BAugmentedConstrainVector(V,stage,incrfactor);
		}
		break;
	}
}

	/*----------------------------------------------------------\
	|	                Handling Constrains						|
	|															|
	|	Simply altering the Stiffness Matrix and Load Vector	|
	|	For more info Look at Bathe Book page 145.				|
	|	The code is straight forward and the notation is		|
	|   A matrix one and a simple to understand					|
	|															|
	\----------------------------------------------------------*/

void Constrains::PenaltyConstrainMatrix(StiffnessMatrix& K,int stage)
{
	int tempindex=0;
	double penalty = fem->GetPenalty();
	Matrix B;B=CreateConstrainsMatrix(stage);
	K=K+penalty*B.t()*B;
}


void Constrains::LagrangeConstrainMatrix(StiffnessMatrix& K,int stage)
{
	Matrix B;B=CreateConstrainsMatrix(stage);
	int m=B.Ncols(); int n=B.Nrows();
	Matrix Zero(n,n); Zero=0;
	K = (K | B.t()) & (B | Zero) ; //& concatenates vertically and | horizontally
}

void Constrains::AugmentedConstrainMatrix(StiffnessMatrix& K,int stage)
{
	double penalty=fem->GetPenalty();
	Matrix B;B=CreateConstrainsMatrix(stage);
	int m=B.Ncols(); int n=B.Nrows();
	Matrix Zero(n,n); Zero=0;
	K = K + penalty*B.t()*B;
	K = (K | B.t()) & (B | Zero) ; //& concatenates vertically and | horizontally
}



void Constrains::PenaltyConstrainVector(ColumnVector& V,int stage)
{
	double penalty = fem->GetPenalty();
	Matrix B;B=CreateConstrainsMatrix(stage);
	ColumnVector Vc;Vc=CreateConstrainsVector(stage);
	V=V+penalty*B.t()*Vc;
}


void Constrains::LagrangeConstrainVector(ColumnVector& V,int stage)
{
	Matrix B;B=CreateConstrainsMatrix(stage);
	ColumnVector Vc;Vc=CreateConstrainsVector(stage);
	V = V & Vc ; //& concatenates vertically and | horizontally
}

void Constrains::AugmentedConstrainVector(ColumnVector& V,int stage)
{
	double penalty = fem->GetPenalty();
	Matrix B;B=CreateConstrainsMatrix(stage);
	ColumnVector Vc;Vc=CreateConstrainsVector(stage);
	V = (V+penalty*B.t()*Vc) & Vc ; //& concatenates vertically and | horizontally
}







Matrix Constrains::CreateConstrainsMatrix(int stage)
{
	int m,n,l;
	m=fem->countfixes();
	n=fem->countkinimatic(stage);
	l=fem->GetNumberOfActiveSystemDofs();
	Matrix B(m+n,l); //this matrix has rows as many as kinematic conditions
	B=0;//Initiallise

	//First i process fixes.
	int tempindex=0;
	for(int i=1;i<=fem->GetNumberofNodes();i++)
	{
		for(int j=1;j<=fem->GetVariables();j++)
		{
			if( (fem->GetNode(i)->GetExist(j)==1) && (fem->GetNode(i)->GetFixed(j)==1) )
			{
				tempindex++;
				B(tempindex,fem->GetGlobalDof(i,j))=1;
			}//endif
		}//end for(j)
	}//end for(i)

	//Now proccessing the kinimatic conditions
	for(i=1;i<=fem->GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= fem->GetFixCommand(i);
		int vecsize = TempCommand->coeff.size();
		if(stage>=TempCommand->startingstage && stage<=TempCommand->endingstage)//Process Fix Command
		{
			double pivot=0;
			for(int l=1;l<=TempCommand->shiftingtimes;l++)
			{
				tempindex++;//new matrix line
				int svalue= TempCommand->shiftvalue;
				for(int j=1;j<=vecsize;j++)
				{
					int dofj = fem->GetGlobalDof(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);
					B(tempindex,dofj)= TempCommand->coeff[j-1];
				}//next coeff
			}//next  shift
		}// end if
	}//Next Command

	return B;
}


ColumnVector Constrains::CreateConstrainsVector(int stage)
{
	int m,n;
	m=fem->countfixes();n=fem->countkinimatic(stage);
	ColumnVector B(m+n); //this vector rows as many as kinematic conditions
	B=0;
	//First i process fixes.
	int tempindex=0;
	for(int i=1;i<=fem->GetNumberofNodes();i++)
	{
		for(int j=1;j<=fem->GetVariables();j++)
		{
			if( (fem->GetNode(i)->GetExist(j)==1) && (fem->GetNode(i)->GetFixed(j)==1) )
			{
				tempindex++;
				B(tempindex)=0;// dof is fixed
			}//endif
		}//end for(j)
	}//end for(i)

	//Now proccessing the kinimatic conditions
	for(i=1;i<=fem->GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= fem->GetFixCommand(i);
		int vecsize = TempCommand->coeff.size();
		if(stage>=TempCommand->startingstage && stage<=TempCommand->endingstage)//Process Fix Command
		{
			for(int l=1;l<=TempCommand->shiftingtimes;l++)
			{
				tempindex++;//new matrix line
				B(tempindex)= TempCommand->result;
			}//next  shift
		}// end if
	}//Next Command

	return B;
}


void Constrains::BConstrainMatrix(StiffnessMatrix& K,int stage)
{
	switch(ConstrainMethod)
	{
	case PENALTY_METHOD :
		{
			BPenaltyConstrainMatrix(K,stage);
		}
		break;
	case LAGRANGE_MULTIPLIERS :
		{
			BLagrangeConstrainMatrix(K,stage);
		}
		break;
	case AUGMENTED_LAGRANGE :
		{
			BAugmentedConstrainMatrix(K,stage);
		}
		break;
	}
}


Matrix Constrains::CreateBConstrainsMatrix(int stage)
{
	Matrix Bc;Bc=CreateConstrainsMatrix(stage);
	ColumnVector Vc;Vc=CreateConstrainsVector(stage);
	int countn=0;
	Matrix B(1,Bc.Ncols()); // Create B with ONLY enties when V(i)!=0
	B=0;
	for(int i=1;i<=Vc.Nrows();i++)
	{
		if(Vc(i)==0)
		{
			countn++;
			if(countn==1)
			{
				B.Row(1) = Bc.Row(i);
			}
			else
			{
				B = B & Bc.Row(i);
			}
		}
	}
	return B;
}


void Constrains::BPenaltyConstrainMatrix(StiffnessMatrix& K,int stage)
{
	int tempindex=0;
	double penalty = fem->GetPenalty();
	Matrix B;B=CreateBConstrainsMatrix(stage);
	K=K+penalty*B.t()*B;
}


void Constrains::BLagrangeConstrainMatrix(StiffnessMatrix& K,int stage)
{
	Matrix B;B=CreateBConstrainsMatrix(stage);
	int m=B.Ncols(); int n=B.Nrows();
	Matrix Zero(n,n); Zero=0;
	K = (K | B.t()) & (B | Zero) ; //& concatenates vertically and | horizontally
}

void Constrains::BAugmentedConstrainMatrix(StiffnessMatrix& K,int stage)
{
	double penalty=fem->GetPenalty();
	Matrix B;B=CreateBConstrainsMatrix(stage);
	int m=B.Ncols(); int n=B.Nrows();
	Matrix Zero(n,n); Zero=0;
	K = K + penalty*B.t()*B;
	K = (K | B.t()) & (B | Zero) ; //& concatenates vertically and | horizontally
}


ColumnVector Constrains::GetConstrainLoadVector(int stage,double incrfactor)
{
	switch(ConstrainMethod)
	{
	case PENALTY_METHOD :
		{
			return GetPenaltyConstrainLoadVector(stage,incrfactor);
		}
		break;
	case LAGRANGE_MULTIPLIERS :
		{
			return GetLagrangeConstrainLoadVector(stage,incrfactor);
		}
		break;
	case AUGMENTED_LAGRANGE :
		{
			return GetAugmentedConstrainLoadVector(stage,incrfactor);
		}
		break;
	}
	return 0;
}


ColumnVector Constrains::GetPenaltyConstrainLoadVector(int stage,double incrfactor)
{
	ColumnVector DU;
	double penalty=fem->GetPenalty();
	DU=fem->CCenter->GetCurrentDisplacement()-fem->CCenter->GetDisplacementPreviousStage();
	Matrix B;B=CreateConstrainsMatrix(stage);
	ColumnVector V;V=CreateConstrainsVector(stage);
//	printmatrix("b",B,0);
//	printvector("Vpen",B.t()*incrfactor*V,0);
//	printvector("BDU",B.t()*B*DU,0);
//	printvector("penu2",penalty*B.t()*(incrfactor*V),0);
	DU= penalty*B.t()*(incrfactor*V-B*DU);
//	DU= -penalty*B.t()*B*DU;
	return DU;
}


ColumnVector Constrains::GetLagrangeConstrainLoadVector(int stage,double incrfactor)
{
	ColumnVector DL,DU;
	int m=fem->GetNumberOfActiveSystemDofs();
	DL=fem->CCenter->GetCurrentDisplacement()-fem->CCenter->GetDisplacementPreviousStage();
	DL=DL.Rows(m+1,DL.Nrows());
//	printvector("DL",DL,0);
	Matrix B;B=CreateConstrainsMatrix(stage);
	DU=(fem->CCenter->GetCurrentDisplacement()-fem->CCenter->GetDisplacementPreviousStage()).Rows(1,m);
	DL= B.t()*(-DL);
//	printvector("BDL",DL,0);
	DL= DL & B*DU;
	return DL;
}


ColumnVector Constrains::GetAugmentedConstrainLoadVector(int stag,double incrfactor)
{
	return 0;
}

void Constrains::BPenaltyConstrainVector(ColumnVector& V,int stage,double incrfactor)
{
	ColumnVector DU;
	double penalty=fem->GetPenalty();
	DU=fem->CCenter->GetCurrentDisplacement()-fem->CCenter->GetDisplacementPreviousStage();
	Matrix B;B=CreateConstrainsMatrix(stage);
	ColumnVector Vc;Vc=CreateConstrainsVector(stage);
//		printvector("Vc",incrfactor*Vc,0);
//		printvector("BDU",B*DU,0);
	Vc=incrfactor*Vc-B*DU;
//		printvector("Vc",penalty*B.t()*Vc,0);
	V=V+penalty*B.t()*Vc;
}


void Constrains::BLagrangeConstrainVector(ColumnVector& V,int stage,double incrfactor)
{
	int m,n,l;
	m=fem->countfixes();n=fem->countkinimatic(stage);
	l=fem->GetNumberOfActiveSystemDofs();
	ColumnVector Vc(m+n);Vc=0;
	V= V & Vc;			//& concatenates vertically and | horizontally
}

void Constrains::BAugmentedConstrainVector(ColumnVector& V,int stage,double incrfactor)
{
	ColumnVector Vc;Vc=CreateConstrainsVector(stage);
	ColumnVector one(1); one=0;
	int countn=0;
	for(int i=1;i<=Vc.Nrows();i++)
	{
		if(Vc(i)==0)
		{
			V = V & one; //& concatenates vertically and | horizontally
		}
	}
}
