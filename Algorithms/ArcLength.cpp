// ArcLength.cpp: implementation of the ArcLength class.
//
//////////////////////////////////////////////////////////////////////

#include "ArcLength.h"
#include "FEM.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ArcLength::ArcLength()
{

}

ArcLength::ArcLength(FEM* pfem)
{
	fem = pfem ;
	cc=pfem->CCenter;
	pc=pfem->PCenter;
	ActiveSystemDofs=fem->GetNumberOfActiveSystemDofs();
}


ArcLength::~ArcLength()
{

}

void ArcLength::Solve()
{
	ArcLengthMethod();
}


void ArcLength::ArcLengthMethod()
{
/*----------------------------------------------------------\
|															|
|					Variable Declaration					|
|															|
\----------------------------------------------------------*/
	ColumnVector V1,V2,V3;
	ColumnVector Zero(ActiveSystemDofs);Zero=0;
	LowerTriangularMatrix Lower(ActiveSystemDofs) ; // for ldl decom
	DiagonalMatrix Diagonal(ActiveSystemDofs);		// for ldl decom
	StrainEnergy=0;
/*----------------------------------------------------------\
|															|
|			 Are we STARTING OR RESTARTING????				|
|															|
\----------------------------------------------------------*/
	if (fem->Restart== false)
	{
		totalloadfactor=0;
		cc->SetCurrentDisplacementZero();
		fem->UpdateNodes(cc->GetCurrentDisplacement());// update structure with node displacements
		fem->CurrentStage=1;
	}
	else
	{
//		ReadRestartParameters(Stifi,StrainEnergy,TotalLoadFactor,CurrentStage,Restart,VectorLength,U,DisplacementAtPreviousStage);//These are Current StiffnessParameter,TotalLoadFactor,Displacements
//		CurrentLoadStep=Restart;
	}

/*----------------------------------------------------------\
|															|
|			 Write Headers To Stiffness File				|
|															|
\----------------------------------------------------------*/
	fstream Stiffness;
	Stiffness.open(fem->CurrentStiffnessFile,ios::app);//Εξωτερικά φορτία
	Stiffness<<"CurrStif\t\t"<<"NegEntries\t"<<"LFactor\t\t\t"<<"StrEner\t\t\t"<<"LStep\t\tArcUsed\t\tArcHint\n";
	flush(Stiffness);
	Stiffness.close();

	for(int l=fem->CurrentStage;l<=fem->NumberOfStages;l++)
	{
		fem->CurrentStage=l;
		if(fem->CurrentStage==1)
		{
			LPS=Zero;
		}
		else
		{
			LPS = fem->AssembleStageLoads(1,fem->CurrentStage-1)+ fem->AssembleLoads();
		}
		
		if(fem->CurrentStage==1)
		{
			LCS = fem->AssembleLoads()+fem->AssembleStageLoads(1);
		}
		else
		{
			LCS = fem->AssembleStageLoads(fem->CurrentStage);
		}


		if(fem->Restart==0)
		{
			totalloadfactor=0;
			oldtotalloadfactor=0;
			cc->SetDisplacementPreviousStage(cc->GetCurrentDisplacement());
		}
		else
		{
			fem->Restart=0;//Next time it will be a normal iteration.....
		}

		cc->SetExternalLoadFromPreviousStage(LPS);
		fem->ConstrainCenter->ConstrainVector(LCS,fem->CurrentStage);
		cc->SetExternalLoad(LCS);
		ColumnVector ZeroFull(LCS.Nrows());ZeroFull=0;				// All system dofs+Lagrange
		cc->SetCurrentExternalLoad(ZeroFull);						//this will count K*DU
		DT=LCS;
		Residual=DT;												//All the external force is a residual
/*----------------------------------------------------------\
|															|
|				MAIN ROUTINE - Stage Iteration				|
|															|
\----------------------------------------------------------*/
		for(int incr=1;incr<=fem->IterationsDesired;incr++)
		{
			oldtotalloadfactor=totalloadfactor;						//Saving old Force for cut increments and old total loadlevel

			K=fem->AssembleStiffnessMatrix();						//Assemble Stiffeness Matrix
			cc->SetStiffnessMatrix(K);								//Storing the stiffness Matrix is not needed for Elastic Solve
			fem->ConstrainCenter->ConstrainMatrix(K,fem->CurrentStage);			//Impose Boundary Conditions
			if (fem->SolverType==1) { Kinv=K.i();}
			switch (fem->SolverType)
			{
			case 1 : //Gauss Elimination
				DT=Kinv*LCS;				
//		printvector("DT",DT,1);
				break;
			case 2:
				//κωδικας για αλλους solvers
				break;
			}//end switch
			ldl(K,Lower,Diagonal,neg,ActiveSystemDofs);				//ldl decomposition to get the negative pivots

/*----------------------------------------------------------\
|															|
|			 Calculate Stiffness Parameter					|
|															|
\----------------------------------------------------------*/		
			fem->CalibrateVector(Residual,true,false);				//Zero the fixed dof coords
//			fem->CalibrateVector(DT,true,false);
			double Stift;Stift = dot(DT,LCS,ActiveSystemDofs);
			dl=dot(DT,DT,ActiveSystemDofs);
			double Stif;Stif=Stift/dl;
			if(fem->CurrentLoadStep==0 && fem->Restart==false)
			{
				Stifi=Stif;
			}
			CStif=Stifi/Stif;
			dl=sqrt(dl);
/*----------------------------------------------------------\
|															|
|			Handle More than one negative points			|
|															|
\----------------------------------------------------------*/		
			if (neg>1)
			{
				/*cout<<"negative pivots "<<neg;flush(cout);
				DiagonalMatrix EM(GetNumberOfActiveSystemDofs());
				EigenValues(K.SubMatrix(1,1,GetNumberOfActiveSystemDofs(),GetNumberOfActiveSystemDofs()),EM);
				fstream eigen;
				eigen.open("eigenvalues.txt",ios::app);//Εξωτερικά φορτία
				eigen.setf(ios::showpoint);eigen.precision(7);
				printvector(eigen,EM,CurrentLoadStep);
				eigen.flush();
				eigen.close();
				exit(0);*/
			}
/*----------------------------------------------------------\
|	Now we call the scaleup routine witch calculates		|
|	the incrloadfactor and scales the Displacement vector	|
|	Then we start equilibrium iterations					|
\----------------------------------------------------------*/		
			ScaleUp();
			iterations=IterateToEquilibrium();
			fem->AssembleInternalLoads(true,fem->Method);//Update stress condition in structure
			//Γράφω τα στοιχεία στο αρχείο αν έχει επιτευχθεί σύγκλιση//
			//Γράφω τα στοιχεία στο αρχείο αν έχει επιτευχθεί σύγκλιση//
/*----------------------------------------------------------\
|															|
|			Writing to files if iterations converged		|
|															|
\----------------------------------------------------------*/		
			if(bet<fem->ErrorTolerance)
			{
				fem->CurrentLoadStep++;
				V1 = (cc->GetDU()).Rows(1,ActiveSystemDofs);
				K=cc->GetStiffnessMatrix();
//				fem->ConstrainCenter->BConstrainMatrix(K,fem->CurrentStage);
				V2 = K*V1;
				StrainEnergy += dot(V1,V2,ActiveSystemDofs );

				Stiffness.open(fem->CurrentStiffnessFile,ios::app);//Εξωτερικά φορτία
				Stiffness.setf(ios::showpoint);Stiffness.precision(7);
				Stiffness<<CStif<<"\t\t"<<neg<<"\t\t"<<totalloadfactor<<"\t\t"<<StrainEnergy<<"\t\t"<<fem->CurrentLoadStep <<"\t\t"<<fem->DesiredLengthIncrement<<"\t\t";
				flush(Stiffness);Stiffness.close();
				ColumnVector TL;//TotalLoad
//				TL=(totalloadfactor*cc->GetExternalLoad()).Rows(1,ActiveSystemDofs)+cc->GetExternalLoadFromPreviousStage();	
//				TL=TL+fem->ConstrainCenter->GetConstrainLoadVector(fem->CurrentStage,incrloadfactor).Rows(1,ActiveSystemDofs);
				TL =fem->AssembleInternalLoads(true,fem->Method,true);
//				loads<<(norm(tt,2,K))/norm(TL,2,K)<<"\n";
				fstream loadfile;loadfile.open(fem->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
				printvector(loadfile,TL,fem->CurrentLoadStep,ActiveSystemDofs);
				loadfile.close();
//				ColumnVector tt=AssembleInternalLoads(true,Method,true)-TL;
//				loads<<(norm(tt,2,K))/norm(TL,2,K)<<"\n";
				loadfile.close();
				fstream displ;displ.open(fem->DisplacementsFileOut,ios::app);//Μετακινήσεις
				printvector(displ,cc->GetCurrentDisplacement(),fem->CurrentLoadStep,ActiveSystemDofs);
				displ.close();
				fem->ReactionToFile(fem->ReactionFileOut);//Αντιδράσεις
				fem->InternalForceToFile(fem->ForceFileOut);//Εσωτερικές Δυνάμεις
				fem->StressAndStainToFile(fem->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
				fem->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
//				WriteToRestartFile(Stifi,StrainEnergy,TotalLoadFactor,CurrentStage,CurrentLoadStep,VectorLength,U);
				//end of writing to file
				cc->SetLastDisplacement(cc->GetCurrentDisplacement());
			}
			if (fem->AutoIncrement==true)
			{
				PrepareNextIncr();
			}
			if(totalloadfactor>=0.99999999999)
			{break;}


		}//end for stage



	}// end for	all stages



}

void ArcLength::PrepareNextIncr()
{
	fstream info;
	double ReductionFactor,RIterations,RIterationsDesired;
	if (iterations >= fem->MaxIterations )
	{
		info.open(fem->InfoFile,ios::app);
		info<<" Δεν Υπήρξε Σύγκλιση. Νεοι συντελεστές Φορτίου και επανεκκίνηση απο τον προηγούμενο σημείο ισορροπίας\n";
		info.close();
		// Did not Converge in Last Increment
		//Calculating Reduction Factor
		ReductionFactor = fem->ErrorTolerance /bet;
		ReductionFactor = Min(ReductionFactor,0.5);
		ReductionFactor = Max(ReductionFactor,0.1);
		cc->SetCurrentDisplacement(cc->GetLastDisplacement()) ;//U=Uprevious
		fem->UpdateNodes(cc->GetCurrentDisplacement().Rows(1,ActiveSystemDofs) );
		fem->AssembleInternalLoads(true,fem->Method);
		incrloadfactor = totalloadfactor - oldtotalloadfactor;
		totalloadfactor = oldtotalloadfactor;
	}
	else //Did  converge.Get Change Factor and make large if no Real iterations on last increament
	{
		ReductionFactor = 1000;
		if (iterations >1)
		{
			RIterations = iterations -1;
			RIterationsDesired = fem->IterationsDesired ;
			ReductionFactor = sqrt(RIterationsDesired / RIterations);
		}
	}
	if (fem->EnableArcLengthMethod  == false)
	{
		incrloadfactor = ReductionFactor *incrloadfactor;
		incrloadfactor = Min(incrloadfactor,fem->MaxLoadIncrement );
		incrloadfactor = Max(incrloadfactor,fem->MinLoadIncrement );
		fstream info;
		info.open(fem->InfoFile,ios::app);
		info<<"\tΝεό ποσοστό αύξησης φορτίου "<<incrloadfactor<<"\n";
		info.close();
	}
	else
	{
		incrloadfactor = ReductionFactor *incrloadfactor;
		incrloadfactor = Min(incrloadfactor,fem->MaxLoadIncrement );
		incrloadfactor = Max(incrloadfactor,fem->MinLoadIncrement );
		fem->DesiredLengthIncrement= ReductionFactor*dl;
		fem->DesiredLengthIncrement= Min(fem->DesiredLengthIncrement,fem->MaxLengthIncrement);
		fem->DesiredLengthIncrement= Max(fem->DesiredLengthIncrement,fem->MinLengthIncrement);
	}
	return;
}



void ArcLength::ScaleUp()
{
	int assign=1;
	incrloadfactor=fem->IncrLoadFactor;
	if (fem->AutoIncrement==true)
	{
		//Checking for switch to Arc Length Method//
		if (fem->SwitchToArcLength==true)
		{
			if( (fem->EnableArcLengthMethod==false) && (CStif<=fem->LowestStiffnessSwitch) )
			{
				fem->EnableArcLengthMethod=true;
				if(fem->DesiredLengthIncrement==0)
					{fem->DesiredLengthIncrement=incrloadfactor*dl;}
				if(fem->MaxLengthIncrement==0)
					{fem->MaxLengthIncrement=5*dl;}
				if(fem->MinLengthIncrement==0)
					{fem->MinLengthIncrement=0.01*dl;}
				fem->SwitchToArcLength=false;
			}
		}
		if (fem->EnableArcLengthMethod==false)
		{
			//Load Control//
			incrloadfactor=fem->IncrLoadFactor;
			dl=dl*incrloadfactor;
		}
		else
		{
			//Arc Length Control//
			//Set up lengths if none were given//
			if(fem->DesiredLengthIncrement==0.0)
			{
				fem->DesiredLengthIncrement=incrloadfactor*dl;
				fem->MaxLengthIncrement=5*dl;
				fem->MinLengthIncrement=0.01*dl;
			}
			//Computing incr load factor FACI
			incrloadfactor=  fem->DesiredLengthIncrement/dl;
			incrloadfactor = Min(incrloadfactor,fem->MaxLoadIncrement );
			incrloadfactor = Max(incrloadfactor,fem->MinLoadIncrement );
		    if(neg==1)
			{
				 assign = -1;
			}
			incrloadfactor=assign*incrloadfactor;
			dl=fem->DesiredLengthIncrement;
		}
	}
	totalloadfactor=totalloadfactor+incrloadfactor;

	ColumnVector U;U=cc->GetLastDisplacement()+fabs(assign*incrloadfactor)*DT;
//	fem->CalibrateVector(U,true,false);
	fem->UpdateNodes(U);
	fem->AssembleInternalLoads(true,fem->Method);
	cc->SetCurrentDisplacement(U);
	cc->SetDisplacementPreviousIteration(cc->GetLastDisplacement());
	cc->IncreaseCurrentExternalLoad(incrloadfactor*cc->GetExternalLoad());
//	printmatrix("CK",cc->GetStiffnessMatrix(),fem->CurrentLoadStep);
//	printvector("CU",U,fem->CurrentLoadStep,Residual.Nrows());
//	printvector("CLoad",cc->GetCurrentExternalLoad(),fem->CurrentLoadStep,Residual.Nrows());
//	printvector("CELoad",cc->GetExternalLoad(),fem->CurrentLoadStep,Residual.Nrows());
/*----------------------------------------------------------\
|															|
|New Updated displacments according to the calculated load  |
|factor. Write Infos To file.								|
|															|
\----------------------------------------------------------*/		
	fstream info;
	info.open(fem->InfoFile,ios::app);
	info<<" O συντελεστής φορτίου είναι "<<totalloadfactor<<"%\n";
	info.close();


}



int ArcLength::IterateToEquilibrium()
{
/*----------------------------------------------------------\
|	Iterate to equilibrium									|
|	Input: Predicted Displacements (U),Force Vector			|
|	Output: New Displacement Vector							|
\----------------------------------------------------------*/		
	int illfail=0;
	ColumnVector V1,DTbar,InternalForceVector,InitialResidual,ReactionVector,Temp,ConstrainVector;																	// Failing Contition
	ColumnVector Zero(ActiveSystemDofs);Zero=0;
	double dispbet=0;double dispbet1=0;
	double bas;double normu=0;double small = 0.001; double normuprev=0;				//norm(Uprevious , fem->NormType , K , Uprevious.Nrows() );
	cc->SetDisplacementPreviousIteration(cc->GetLastDisplacement());				//UPreviousIteration = Uprevious;

	iterations=1;																	// Starting iterations
	while(iterations<=fem->MaxIterations)
	{ 
		fem->UpdateNodes(cc->GetCurrentDisplacement().Rows(1,ActiveSystemDofs) );
		if(iterations ==1)															//Update ArcHint
		{
			fstream info;
			fem->ArcHint = sqrt(dot(cc->GetDU(),cc->GetDU(),ActiveSystemDofs) );
			info.open(fem->CurrentStiffnessFile,ios::app);
			info<<fem->ArcHint<<"\n";flush(info);
			info.close();
	
		}
		if( (iterations==1) || (fem->MaxNoLineSearches==0) || (illfail==1))
		{
			if(fem->FullNR==true)													// Update K if FUll NR
			{
				K=fem->AssembleStiffnessMatrix();									//Assemble Stiffeness Matrix
				cc->SetStiffnessMatrix(K);											//Storing the stiffness Matrix is not needed for Elastic Solve
				fem->ConstrainCenter->ConstrainMatrix(K,fem->CurrentStage);			//Impose Boundary Conditions
				Kinv = K.i();
			}
		}

		InternalForceVector = fem->AssembleInternalLoads(false,fem->Method);		// This are internal forces only due to loads
//		printvector("interV",InternalForceVector,fem->CurrentLoadStep);
		ConstrainVector = fem->ConstrainCenter->GetConstrainLoadVector(fem->CurrentStage,totalloadfactor);
		InternalForceVector-=ConstrainVector.Rows(1,ActiveSystemDofs);
		if(fem->ConstrainCenter->GetConstrainMethod()!=PENALTY_METHOD)
		{
			InternalForceVector= InternalForceVector & ConstrainVector.Rows(ActiveSystemDofs+1,LCS.Nrows());
		}
//		fem->CalibrateVector(InternalForceVector ,true,false);						// Zero at reactions to easily compare with exteral forces
//		printvector("ConLoV",ConstrainVector,fem->CurrentLoadStep);
//		printvector("PLoad",cc->GetExternalLoadFromPreviousStage(),fem->CurrentLoadStep);
		if(fem->CurrentStage==1)
		{
			Residual = cc->GetCurrentExternalLoad()-InternalForceVector;
		}
		else
		{
			Residual = cc->GetCurrentExternalLoad()+cc->GetExternalLoadFromPreviousStage()-InternalForceVector;
		}
		fem->CalibrateVector(Residual,true,false);
//*			DEBUGING STUF
//		printmatrix("K",fem->CCenter->GetStiffnessMatrix(),fem->CurrentLoadStep);
//		printvector("U",cc->GetCurrentDisplacement(),fem->CurrentLoadStep,(cc->GetCurrentDisplacement()).Nrows());
//		printvector("Residual",Residual,fem->CurrentLoadStep);
//		printvector("Load",cc->GetCurrentExternalLoad(),fem->CurrentLoadStep);
//		printvector("PLoad",cc->GetExternalLoadFromPreviousStage(),fem->CurrentLoadStep);
//		printvector("InternalForceVector",InternalForceVector,fem->CurrentLoadStep);
//		printvector("DUP",cc->GetDUPreviousIteration(),fem->CurrentLoadStep);
//*/
		if(iterations==1) 
		{
			InitialResidual = -Residual;
		}
		//Check for convergence
		if (fem->EnableArcLengthMethod == true)
		{
//			DT = cc->GetExternalLoad();
		}
		if (fem->DisplacementControl == false ) //Load control
		{
//			V1 = cc->GetCurrentExternalLoad().Rows(1,ActiveSystemDofs)+ConstrainVector;
			V1 = InternalForceVector;
//		printvector("V1",V1,fem->CurrentLoadStep);
			bas = Max( norm(V1 , fem->NormType , K , ActiveSystemDofs ),small);
		}
		else
		{
			ReactionVector = fem->ReactionVector(false,fem->Method);
			bas = Max( norm(ReactionVector , fem->NormType , K, ActiveSystemDofs) , small );
		}

		bet = norm(Residual,fem->NormType,K,ActiveSystemDofs)/bas;
//		normu=norm(U , fem->NormType , K , U.Nrows() );
//		betu= fabs(normu-normuprev)/normu;normuprev=normu;
//		cout<<" InternForceVector "<<norm(InternalForceVector , fem->NormType ,cc->GetStiffnessMatrix(),ActiveSystemDofs)<<" Load "<<norm(V1 , fem->NormType , cc->GetStiffnessMatrix(), ActiveSystemDofs)<<"\n";flush(cout);
//		cout<<" Residual "<<norm(Residual,fem->NormType,K,ActiveSystemDofs)<<"\n";flush(cout);
		if(iterations==1){dispbet1=norm(cc->GetDUPreviousIteration(),fem->NormType,K,ActiveSystemDofs);}
		dispbet = norm(cc->GetDUPreviousIteration(),fem->NormType,K,ActiveSystemDofs)/dispbet1;
		if(fem->EnableArcLengthMethod==true){dispbet=fem->ErrorTolerance;}
		cout<<" Error "<<bet<<" ErrorDisp "<<dispbet<<"\n";flush(cout);

		if (bet<=fem->ErrorTolerance && dispbet<=fem->ErrorDispTolerance)
		{
			bet=Min(bet,dispbet);
			fstream info;
			info.open(fem->InfoFile,ios::app);
			info<<"\tEπιτεύχθηκε Συγκλιση. Το αριθμητικό σφάλμα είναι ίσο με "<<bet<<"\n";
			info.close();
			return iterations;
		}
		else
		{
			fstream info;
			info.open(fem->InfoFile,ios::app);
			info<<"\tΤο αριθμητικό σφάλμα είναι ίσο με "<<bet<<"\n";
			info.close();
		}
		if (fem->FullNR ==true && iterations!=1)
		{
			K= fem->AssembleStiffnessMatrix();
			cc->SetStiffnessMatrix(K);											//Storing the stiffness Matrix is not needed for Elastic Solve
			fem->ConstrainCenter->ConstrainMatrix(K,fem->CurrentStage);			//Impose Boundary Conditions
			Kinv = K.i();
		}
		if ( (fem->EnableArcLengthMethod == true) && (fem->FullNR == true) )
		{
			DT = Kinv*LCS.Rows(1,Kinv.Nrows());// NA TO DW
//			ColumnVector tt = fem->AssembleLoads()+fem->AssembleStageLoads(1);
//			fem->ConstrainCenter->BPenaltyConstrainVector(tt,fem->CurrentStage,totalloadfactor);
//			DT = Kinv*tt;// NA TO DW
//			fem->CalibrateVector(DT,true,true);
		}
//		printvector("Residual",Residual,0);
//		fem->ConstrainCenter->BConstrainVector(Residual,fem->CurrentStage,totalloadfactor);
//		printvector("Residual",Residual,0);
		DTbar = Kinv*Residual.Rows(1,Kinv.Nrows());
		fem->CalibrateVector(DTbar,true,false);

//*			DEBUGING STUF
//		printmatrix("Kinv",Kinv,fem->CurrentLoadStep);
//		printmatrix("K",K,fem->CurrentLoadStep);
//		printvector("DTbar",DTbar,fem->CurrentLoadStep);
//		printvector("Load",cc->GetCurrentExternalLoad(),fem->CurrentLoadStep);
//		printvector("InternalForceVector",InternalForceVector,fem->CurrentLoadStep);
//		printvector("DTbar",DTbar,fem->CurrentLoadStep);
//		printvector("Uprin",cc->GetCurrentDisplacement(),fem->CurrentLoadStep);
//		printvector("DT",DT,1);
//*/
		if (fem->EnableArcLengthMethod ==true)
		{
			ArcLengthRoutine(DT,DTbar,illfail,iterations);
//		printvector("U",cc->GetCurrentDisplacement(),fem->CurrentLoadStep);
//		printvector("DUP",cc->GetDUPreviousIteration(),fem->CurrentLoadStep);
			if (illfail==2)
			{
				iterations= fem->MaxIterations ;
				//WRITE WHY
			}
		}
		else
		{
			if (fem->Accelerate ==true)
			{
				Accelerate();
			}
			ColumnVector U;U=cc->GetCurrentDisplacement();
			cc->SetDisplacementPreviousIteration(U);								//Update Previous Iteration
			U = U + DTbar;															//Update Displacements
//			fem->CalibrateVector(U,true,false);
			fem->UpdateNodes(U.Rows(1,ActiveSystemDofs) );
			fem->AssembleInternalLoads(true,fem->Method);
			cc->SetCurrentDisplacement(U);
//		printvector("Umeta",cc->GetCurrentDisplacement(),fem->CurrentLoadStep);
			if (fem->MaxNoLineSearches > 0)
			{
//				LineSearchLoop(fem,U,DTbar,Residual,totalloadfactor,ExternalLoad,K,illfail,Slol);
			}
			if (illfail==2)
			{
				//WRITE WHY
				fstream info;
				info.open(fem->InfoFile,ios::app);
				info<<" Η Line Search δεν ήτανε επιτυχής. Ξεπεράστηκε ο επιτρεπτός αριθμός βημάτων\n";
				info.close();
				iterations= fem->MaxIterations ;
			}
		}
		//Write results if needed to check every step of displacements
		iterations++;
	}//End while
	//Max iterations reached . No Converge Write Info
	if (fem->AutoIncrement ==true)
	{
		return iterations;
	}
	else
	{
		fstream info;
		info.open(fem->InfoFile,ios::app);
		info<<"\t\tΔεν επιτεύχθηκε Συγκλιση\n";
		info.close();
		exit(1);
	}
}


void ArcLength::ldl (Matrix A,LowerTriangularMatrix &L,DiagonalMatrix &D,int &neg, int n) 
{
    int i, j, k;
	RowVector v(n);
	neg=0;
// Indices are zero-based in C, but we start with row 1
// because the 0th element is correct already 
	D(1)=A(1,1);
    for (i = 1;i <= n;i++) 
	{
		L(i,i)=1;
        for (j = 1;j<= i - 1;j++)
		{
            if (A(j,j) == 0.0)
                D(i) = 0.0;
            else 
			{
                v(j)=L(i,j)*D(j);
			}
		}
		D(i)=A(i,i);
        for (j = 1;j<= i - 1;j++)
		{
			 D(i)-=L(i,j)*v(j);
		}
        for (j = i+1 ;j<= n;j++)
		{
			 L(j,i)=A(j,i);
             for(k=1;k<=i-1;k++)
			 {
				 L(j,i)-=L(j,k)*v(k);
			 }
			 L(j,i)=L(j,i)/D(i);
		}
     }
	for(i=1;i<=n;i++)
	{
		if(D(i)<0)
		{
			neg++;
		}
	}
 }


void ArcLength::SolveQuadraticEquation(double a,double b,double c, double &R1,double &R2,double &Rlinear,int &failcondition)
{
	//failcondition : 0 : 2 Real Roots
	//              : 1 : 1 Real Root
	//              : 2 : No Real Roots
	double diakrinousa;
	if (a == 0) //Linear Equation
	{
		if (b==0) //No equation.No Roots
		{
			failcondition=2;
			return;
		}
		Rlinear = -c/b;
	}
	else 
	{
		diakrinousa = b*b - 4*a*c;
	}
	if (diakrinousa ==0)
	{
		failcondition = 1;
		Rlinear = - 0.5*b/a;
	}
	else if(diakrinousa<0)
	{
		failcondition=2;
	}
	else
	{
		failcondition = 0;
		R1=(-b + sqrt(diakrinousa))/(2*a);
		R2=(-b - sqrt(diakrinousa))/(2*a);
	}
	return;
}


void ArcLength::ArcLengthRoutine(ColumnVector &DT,ColumnVector &DELBAR,int &illfail,int iterations)
{
/*--------------------------------------------------------------\
|				ArcLengthMethod Routine							|
|	ArcLengthAlgorithm = 1 : Crisfield' s ArcLength				|
|	ArcLengthAlgorithm = 2 : Rams's ArcLength					|
|	ArcLengthAlgorithm = 3 : Frieds's ArcLength					|
|	ArcLengthAlgorithm = 4 : Reinbolt's ArcLength				|
|	ArcLengthAlgorithm = 5 : Powells's Work Control ArcLength	|
\--------------------------------------------------------------*/		

	
	ColumnVector W1,UChange,DPbar;
	
	fstream info;
	info.open(fem->InfoFile,ios::app);
	info<<"\t\tEntering ArcLength Method\n";flush(info);
	info.close();
	UChange = cc->GetDUPreviousIteration();
	double t0,t1,t2,t3,Solution;
//	int k;//Fastest Changing component in Dangential Predictor
	if (fem->ArcLengthAlgorithm==1)
	{
		double Root1,Root2,RLin;
		DPbar = cc->GetDUPreviousIteration()+DELBAR;
		double a,b,c,c1,a4,a5,l2,d;
		a=0;b=0;a4=0;a5=0;d=incrloadfactor/10.;
		a=dot(DT,DT,ActiveSystemDofs );
		b=2*dot(DT,DPbar,ActiveSystemDofs );
		c1=dot(DPbar,DPbar,ActiveSystemDofs );
/*		if(iterations)
		{
			l2= -dot(UChange,UChange,ActiveSystemDofs );arc=l2;
		}
		else
		{
			l2=arc;
		} //*/
		l2=fem->DesiredLengthIncrement*fem->DesiredLengthIncrement*(-1);
		c=c1+l2;
		a4=dot(DPbar,UChange,ActiveSystemDofs);
		a5=dot(DT,UChange,ActiveSystemDofs );
//		printvector("DT",DT,1,ActiveNumberOfDofs);
//		printvector("U",U,1,ActiveNumberOfDofs);
//		printvector("Uprevious",UPreviousIteration,1,ActiveNumberOfDofs);
//		printvector("DELBAR",DELBAR,1,ActiveNumberOfDofs);
//		printvector("UChange",UChange,1,ActiveNumberOfDofs);
//		printvector("DPbar",DPbar,1,ActiveNumberOfDofs);
//		info<<"b = " <<b;
		SolveQuadraticEquation(a,b,c,Root1,Root2,RLin,illfail);
		if(illfail==2)
		{
			//write to file that we are anable to find a solution
			info.open(fem->InfoFile,ios::app);
			info<<"\t\tΔεν βρέθηκε ρίζα της εξίσωσης 2ου Βαθμού\n";flush(info);
			info.close();
			return;
		}
		else if (illfail ==1)
		{
			Solution = RLin;
		}
		else
		{
			double cost1,cost2;
			cost1 = a4+a5*Root1;
			cost2 = a4+a5*Root2;
			Solution = Root1;
			if( (cost2>0) && (cost1<0) )
			{
				Solution = Root2;
			}
			if( ( (cost1>0) && (cost2>0) ) || ( (cost1<0) && (cost2<0) ) )
			{
				RLin= -c/b;
				if ( fabs(RLin-Root1)>fabs(RLin-Root2))
				{
					Solution = Root2;
				}
				else
				{
					Solution = Root1;
				}
			}
		}
	}//End of Crisfields Arc Length Rootine
	else if(fem->ArcLengthAlgorithm == 2)
	{
		W1 = cc->GetDUPreviousIteration();;
		W1 = (1/dot(W1,W1,ActiveSystemDofs ))*W1;
		t0=dot(W1,DT,ActiveSystemDofs );
		t3=dot(W1,DELBAR,ActiveSystemDofs );
		Solution = -t3/t0;
	}//End of Rams's Arc Length
	else if(fem->ArcLengthAlgorithm == 3)
	{
		W1 = DT;
		W1 = (1/dot(W1,W1,ActiveSystemDofs))*W1;
		t0=dot(W1,DT,ActiveSystemDofs);
		t3=dot(W1,DELBAR,ActiveSystemDofs );
		Solution = -t3/t0;
	}//End of Frieds's Arc Length
	else if(fem->ArcLengthAlgorithm == 4)
	{
	}//End of ReinBolds's Arc Length
	else if(fem->ArcLengthAlgorithm == 5)
	{
//		t1=dot(Loads,DT,ActiveSystemDofs );
//		t2=dot(Loads,DELBAR,ActiveSystemDofs );
//		Solution = -t2/t1;
	}//End of Powells's Arc Length
	totalloadfactor += Solution;
	incrloadfactor+=Solution;

	ColumnVector U;U=cc->GetCurrentDisplacement();
	cc->SetDisplacementPreviousIteration(U);								//Update Previous Iteration
	U = U + DELBAR + Solution*DT;											//Update Displacements
//	fem->CalibrateVector(U,true,false);
	fem->UpdateNodes(U.Rows(1,ActiveSystemDofs) );
	fem->AssembleInternalLoads(true,fem->Method);
	cc->SetCurrentDisplacement(U);
	cc->IncreaseCurrentExternalLoad(Solution*cc->GetExternalLoad());

	//Writing Infos
	info.open(fem->InfoFile,ios::app);
	info<<"\t\tΗ τιμή της Δλ είναι "<<Solution<<" και ο συντελεστής φορτίου είναι "<<totalloadfactor<<"%\n";flush(info);
	info.close();

	return;
}
