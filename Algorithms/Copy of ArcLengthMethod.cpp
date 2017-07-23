#include "FEM.h"
//#include "FemIncludes.h"
//#include <fstream.h>
//#include "Newmatio.h" // for testing


void ldl (Matrix A,LowerTriangularMatrix &L,DiagonalMatrix &D,int &neg, int n) 
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



void ReadRestartParameters(double& InitialStiffness,double& StrainEnergy,double& TotalLoadFactor,int& CurrentStage,int RestartingLoadStep,int VectorLength,ColumnVector& U,ColumnVector& DispPrev)
{
	int CurrentLoadStep=0;
	U=0;
	fstream restart;
	restart.open("restart.txt",ios::in);//Μετακινήσεις
	double previousreadstage=0;
	while(!restart.eof())
	{
		restart>>InitialStiffness>>StrainEnergy>>TotalLoadFactor>>CurrentStage>>CurrentLoadStep>>VectorLength;
		//check if i am entering changing stage to note the last disp vector of the previous stage.
		if(previousreadstage<CurrentStage && CurrentLoadStep==1)
		{
			DispPrev=U;
			previousreadstage=CurrentStage;
		}
		for(int i = 1;i<=VectorLength;i++)
		{
			int temp;
			restart>>temp>>temp>>U(i);
		}
		if(RestartingLoadStep==CurrentLoadStep)
		{
			break;
		}
	}//end while

	restart.close();
	return;
}

void WriteToRestartFile(double InitialStiffness,double StrainEnergy,double TotalLoadFactor,int CurrentStage,int CurrentLoadStep,int VectorLength,ColumnVector U)
{
	fstream restart;
	restart.open("restart.txt",ios::app);//Μετακινήσεις
	restart.width(25);restart.precision(6);restart<<InitialStiffness;
	restart.width(25);restart.precision(6);restart<<StrainEnergy;
	restart.width(10);restart.precision(6);restart<<TotalLoadFactor;
	restart.width(10);restart.precision(6);restart<<CurrentStage;
	restart.width(10);restart.precision(6);restart<<CurrentLoadStep;
	restart.width(10);restart.precision(6);restart<<VectorLength<<"\n";
	printvector(restart,U,CurrentLoadStep,VectorLength);
	restart.close();

}

void SolveQuadraticEquation(double a,double b,double c, double &R1,double &R2,double &Rlinear,int &failcondition)
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


void ArcLengthRoutine(FEM *fem,ColumnVector &DT,ColumnVector &DELBAR,double &TotalLoadFactor,ColumnVector &U,ColumnVector &UPreviousIteration,ColumnVector Loads,int &illfail,int iterations)
{
	// ArcLengthMethod Routine
	// ArcLengthAlgorithm = 1 : Crisfield' s ArcLength
	// ArcLengthAlgorithm = 2 : Rams's ArcLength
	// ArcLengthAlgorithm = 3 : Frieds's ArcLength
	// ArcLengthAlgorithm = 4 : Reinbolt's ArcLength
	// ArcLengthAlgorithm = 5 : Powells's Work Control ArcLength
	int ActiveNumberOfDofs = fem->GetNumberOfActiveSystemDofs();
	fstream info;
	info.open(fem->InfoFile,ios::app);
	info<<"\t\tEntering ArcLength Method\n";flush(info);
	info.close();
	ColumnVector W1;
	ColumnVector UChange;
	UChange = U- UPreviousIteration ;
	double t0,t1,t2,t3,Solution;
//	int k;//Fastest Changing component in Dangential Predictor
	if (fem->ArcLengthAlgorithm ==1)
	{
		double Root1,Root2,RLin;
		ColumnVector DPbar;
		DPbar = U-UPreviousIteration+DELBAR;
		double a,b,c,a4,a5;
		c=fem->DesiredLengthIncrement*fem->DesiredLengthIncrement*(-1);
//		c= -dot(UChange,UChange,ActiveNumberOfDofs );
		a=0;b=0;a4=0;a5=0;
		a=dot(DT,DT,ActiveNumberOfDofs );
		b=2*dot(DT,DPbar,ActiveNumberOfDofs );
		c+=dot(DPbar,DPbar,ActiveNumberOfDofs );
		a4=dot(DPbar,UChange,ActiveNumberOfDofs );
		a5=dot(DT,UChange,ActiveNumberOfDofs );
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
		W1 = U-UPreviousIteration;
		W1 = (1/dot(W1,W1,ActiveNumberOfDofs ))*W1;
		t0=dot(W1,DT,ActiveNumberOfDofs );
		t3=dot(W1,DELBAR,ActiveNumberOfDofs );
		Solution = -t3/t0;
	}//End of Rams's Arc Length
	else if(fem->ArcLengthAlgorithm == 3)
	{
		W1 = DT;
		W1 = (1/dot(W1,W1,ActiveNumberOfDofs))*W1;
		t0=dot(W1,DT,ActiveNumberOfDofs);
		t3=dot(W1,DELBAR,ActiveNumberOfDofs );
		Solution = -t3/t0;
	}//End of Frieds's Arc Length
	else if(fem->ArcLengthAlgorithm == 4)
	{
	}//End of ReinBolds's Arc Length
	else if(fem->ArcLengthAlgorithm == 5)
	{
		t1=dot(Loads,DT,ActiveNumberOfDofs );
		t2=dot(Loads,DELBAR,ActiveNumberOfDofs );
		Solution = -t2/t1;
	}//End of Powells's Arc Length
	TotalLoadFactor += Solution;
	U = U + DELBAR + Solution*DT;
	fem->CalibrateVector(U,true,false);
	fem->UpdateNodes(U.Rows(1,ActiveNumberOfDofs) );
	fem->AssembleInternalLoads(true,fem->Method);
	//Writing Infos
	info.open(fem->InfoFile,ios::app);
	info<<"\t\tΗ τιμή της Δλ είναι "<<Solution<<" και ο συντελεστής φορτίου είναι "<<TotalLoadFactor<<"%\n";flush(info);
	info.close();

	return;
}

void Accelerate()
{
}


void scaleup(FEM *fem,double &incrloadfactor,double &totalloadfactor,double &dl,double &cstif,int &neg,ColumnVector &U,ColumnVector &Uprevious,ColumnVector &DT)
{
	int assign=1;
	///////////////////////////////////////////////////////////
	///Updates LoadFactor FACT via increase FACI //////////////
	///Updates total displacement U=Pt/////////////////////////
	///////////////////////////////////////////////////////////
	if (fem->AutoIncrement==true)
	{
		//Checking for switch to Arc Length Method//
		if (fem->SwitchToArcLength==true)
		{
			if( (fem->EnableArcLengthMethod==false) && (cstif<=fem->LowestStiffnessSwitch) )
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
//				 assign = -1;
			}
			incrloadfactor=assign*incrloadfactor;
			dl=fem->DesiredLengthIncrement;
		}
	}
	totalloadfactor=totalloadfactor+incrloadfactor;
//	if(totalloadfactor>1)
//	{
//		incrloadfactor=incrloadfactor - (totalloadfactor-1);
//		totalloadfactor= 1;
//	}//dont go over Factor 1,0
	U=Uprevious+incrloadfactor*DT;
	fem->CalibrateVector(U,true,false);
	fem->AssembleInternalLoads(true,fem->Method);
	//New Updated displacments according to the calculated load factor  
    //Write Infos To file
		//Writing Infos
	fstream info;
	info.open(fem->InfoFile,ios::app);
	info<<" O συντελεστής φορτίου είναι "<<totalloadfactor<<"%\n";
	info.close();

}//end of scale up routine//


void Search(FEM *fem,double *Prodr,double *ETA, int LineSearchStep, int &ico)
{
	// Obtain No of Previous S-L with negative ratio nearest the origin
	// as well as MAX Previous Step Length
	//If No negative Products INEG=999
	int ineg = 999;
	int ils = LineSearchStep ; // for abbreviation
	int IPOS;
	double ETANEG = 1.0E5;
	double ETAMAXPOSITIVE = 0.0;
	for(int i = 1;i<=ils+1 ;i++)
	{
		if(ETA[i]>ETAMAXPOSITIVE)
		{
			ETAMAXPOSITIVE = ETA[i];
		}
		if( (Prodr[i]<0) && (ETA[i]<= ETANEG) )
		{
			ETANEG = ETA[i]; //Negative Product Found
			ineg = i;
		}
	}//END FOR
	//Begin Interpolation Step
	if( ineg < 999)
	{
		// find IPOS = No of Previous S-L with Pos Ratio that is Closest to ineg
		//But with smaller S-L
		IPOS = 1;
		for(i = 1;i<=ils;i++)
		{
			if ( (Prodr[i]>=0) && (ETA[i]<=ETA[ineg]) && (ETA[i]>=ETA[IPOS]) )
			{
				IPOS = i;
			}
		}//Next i
		//Now we are interpolationg to find  ETAINT (erpolation)
		double EtaInterpolation;
		EtaInterpolation = Prodr[ineg]*ETA[IPOS] - Prodr[IPOS]*ETA[ineg];
		EtaInterpolation = EtaInterpolation / (Prodr[ineg] - Prodr[IPOS] );
		// Computing ETAALT (ernative) to ensure a resonable change
		double EtaAlternative;
		EtaAlternative = ETA[IPOS] +0.2 * (ETA[ineg]-ETA[IPOS]);
		//Take Max of Etainterpolation and Etaalternative
		EtaInterpolation = Max(EtaInterpolation ,EtaAlternative );
		//Check if Lower than Min Step Length
		if (EtaInterpolation <= fem->MinTotalStepLength)
		{
			EtaInterpolation = Max(EtaInterpolation ,fem->MinTotalStepLength );
			if (ico == 1)
			{
				ico++;
				fstream info;
				info.open(fem->InfoFile,ios::app);
				info<<"\t\t\tMin step Length reached twice.";
				info.close();
				//Write that Min step Length Reached Twice
			}
			else if ( ico == 0 )
			{
				ico++;
			}
		}
		ETA[ils+2]=EtaInterpolation;
		// write infos to file
		return;
	}//end if ineg<>99
	// if out here then ineg = 999
	//Proceed with extrapolation
	double EtaExtrapolation,Etamxt;
	Etamxt  = fem->MaxAmpAtAnyStemp *ETAMAXPOSITIVE ;
	Etamxt = Min ( Etamxt,fem->MaxTotalStepLength );
	//Extrapolationg between current and previous step
	EtaExtrapolation = Prodr[ils+1]*ETA[ils] - Prodr[ils]*ETA[ils+1];
	EtaExtrapolation = EtaExtrapolation / (Prodr[ils+1]-Prodr[ils]);
	ETA[ils+2] = EtaExtrapolation ;
	//Check if Value is accepted
	if ( ( EtaExtrapolation <= 0 ) || (EtaExtrapolation > Etamxt) )
	{
		ETA[ils+2]=Etamxt;
	}
	if( (ETA[ils+2]==fem->MaxTotalStepLength ) && (ico == 1) )
	{
		// write to file that we reached Max step Length twice
		fstream info;
		info.open(fem->InfoFile,ios::app);
		info<<"\t\t\tMax step Length reached twice.";
		info.close();
		ico = 2;
		return;
	}
	if( (ETA[ils+2]==fem->MaxTotalStepLength ) && (ico == 0) )
	{
		// write to file that we reached Max step Length once
		ico = 1;
		return;
	}

}// end of search routine




void LineSearchLoop(FEM * fem,ColumnVector &U,ColumnVector &PBar,ColumnVector &Residual,double &totalloadfactor,ColumnVector &Load,StiffnessMatrix &K,int &illfail,double &Slol)
{
	fstream info;
	ColumnVector TempVector;
	double Prodr[20],ETA[20];
	double So;
	illfail = 0;
	So = -dot(PBar,Residual,fem->GetNumberOfActiveSystemDofs() );
	if( So >= 0)
	{
		illfail = 1; // We are Uphill.No need For Line Search
		info.open(fem->InfoFile,ios::app);
		info<<"\t\tWe are Uphill. Exiting LS\n";
		info.close();
		return;
	}
	//Preparing Produact Ratios and step Lengths
	Prodr[1]= 1.0;ETA[1]=0.0;ETA[2]=1.0;
	int ico;// A Counter that will become 
	        // : 1 when max or min S-Length reached or
	        // : 2 When twice reached
	ico=0;
	ColumnVector PTO,InternalForceVector;
	double Seta;
	PTO = U - PBar;
	for(int ils = 1 ; ils <= fem->MaxNoLineSearches ;ils++)
	{
		//Write to info file that we are in LS number ils
		info.open(fem->InfoFile,ios::app);
		info<<"\t\tLine Search No "<<ils;
		///////////////////////////////////////////////////////
		Seta = 0.0;
		InternalForceVector = fem->AssembleInternalLoads(true,fem->Method);
		TempVector = InternalForceVector-totalloadfactor*Load;
		fem->CalibrateVector(TempVector,true,true);
		Seta+=dot(PBar,TempVector,fem->GetNumberOfActiveSystemDofs() );
		Seta = Seta/So;
		Prodr[ils+1] = Seta;
		///////////////////////////////////////////////////
		////////////Write to file The Product Ratio////////
		///////////////////////////////////////////////////
		info<<"\tRatio "<<Seta<<"\n";
		info.close();
		if(fabs(Seta)<=fem->ToleranceOnRatio)
		{
			Slol = ETA[ils+1];
			if(fem->FullNR==true)
			{
				K= fem->AssembleStiffnessMatrix();
//				fem->Penalty(K);
				fem->LagrangeStiffnessMatrix(K);
			}
			return;
		}//endif
		Search(fem,Prodr,ETA,ils,ico);
		if(ico==2)
		{
			if(fem->AutoIncrement == false)
			{
				//Write why we are exiting
				exit(1);
			}
			else
			{
				illfail = 2;
				return;
			}
		}
		//So Ico<>2
		U = PTO + ETA[ils+2]*PBar;
		fem->CalibrateVector(U,true,false);
		fem->UpdateNodes(U.Rows(1,fem->GetNumberOfActiveSystemDofs()) );
			info.open(fem->InfoFile,ios::app);
			info<<"\t\tΕίναι U = U - ( "<<ETA[ils+2]<<" - 1) * PBAR\n";
			info.close();
	}//next ils
	fem->AssembleInternalLoads(true,fem->Method);
	//If we manage to reach here then Max Number of L-S was reached
	if(fem->AutoIncrement==false)
	{
		//Write why we are stoping
		exit(1);
	}
	else
	{
		illfail = 2;
		return;
	}
	//end of story
}





int IterateToEquilibrium(FEM *fem,StiffnessMatrix &K,StiffnessMatrix &Kinv,double &incrloadfactor,double &totalloadfactor,double &dl,ColumnVector &U,ColumnVector &Uprevious,ColumnVector &DT,ColumnVector &ExternalLoad,ColumnVector &LPS,double &bet)
{
	/// Iterate To Equilibrium//
	//Input: Predicted Displacements (U),Force Vector
	//Output: New Displacement Vector
	int illfail=0; // Failing Contition
	int ActiveDofs = fem->GetNumberOfActiveSystemDofs();int iterations =1;
	double bas,Slol;double normu=0;double small = 0.001; double normuprev=0;//norm(Uprevious , fem->NormType , K , Uprevious.Nrows() );
	ColumnVector Residual;
	ColumnVector InitialResidual;
	ColumnVector InternalForceVector;
	ColumnVector ReactionVector;
	ColumnVector DTbar;
	ColumnVector Temp;
	ColumnVector DU;
	ColumnVector UPreviousIteration;
	ColumnVector UNonKinimatic(U.Nrows());
	ColumnVector LagrangeResidual(U.Nrows());
	UPreviousIteration = Uprevious;
	while(iterations<=fem->MaxIterations)
	{ 
		DU=U-UPreviousIteration;
		LagrangeResidual.Rows(1,ActiveDofs) = K.SubMatrix(1,ActiveDofs,1,ActiveDofs)*DU.Rows(1,ActiveDofs);
		fem->CalibrateVector(LagrangeResidual,true,true);
		printmatrix("K",K,1);
		fem->UpdateNodes(U.Rows(1,ActiveDofs) );
		if(iterations ==1)//Update ArcHint
		{
			fstream info;
			fem->ArcHint = sqrt(dot(DU,DU,ActiveDofs) );
			info.open(fem->CurrentStiffnessFile,ios::app);
			info<<fem->ArcHint<<"\n";flush(info);
			info.close();
	
		}
		if( (iterations==1) || (fem->MaxNoLineSearches==0) || (illfail==1))
		{
			if(fem->FullNR==true)// Update K if FUll NR
			{
				K= fem->AssembleStiffnessMatrix();
				fem->LagrangeStiffnessMatrix(K);
			}
		}
		InternalForceVector = fem->AssembleInternalLoads(false,fem->Method);// this are internal forces only due to loads
		fem->CalibrateVector(InternalForceVector ,true,true);// Zero at reactions to easily compare with exteral forces
		Residual = totalloadfactor*ExternalLoad+LPS+LagrangeResidual-InternalForceVector;
		fem->CalibrateVector(Residual,true,false);
//*			DEBUGING STUF
		printvector("U",U,1,U.Nrows());
		printvector("K**U",LagrangeResidual,1,ActiveDofs);
//		printvector("UNonKinimatic",UNonKinimatic,1,UNonKinimatic.Nrows());
		printvector("Residual",Residual,1,Residual.Nrows());
		printvector("Load",totalloadfactor*ExternalLoad+LPS,1,Residual.Nrows());
		printvector("InternalForceVector",InternalForceVector,1,InternalForceVector.Nrows());
//*/
		if(iterations==1) 
		{
			InitialResidual = -Residual;
		}
		//Check for convergence
		if (fem->EnableArcLengthMethod == true)
		{
			DT = ExternalLoad;
		}
		if (fem->DisplacementControl == false ) //Load control
		{
			Temp = totalloadfactor*ExternalLoad+LPS+LagrangeResidual;//CurrentLoad
			bas = Max( norm(Temp , fem->NormType , K , ActiveDofs ),small);
		}
		else
		{
			ReactionVector = fem->ReactionVector(false,fem->Method);
			bas = Max( norm(ReactionVector , fem->NormType , K, ActiveDofs) , small );
		}

		bet = norm(Residual,fem->NormType,K,ActiveDofs)/bas;
//		normu=norm(U , fem->NormType , K , U.Nrows() );
//		betu= fabs(normu-normuprev)/normu;normuprev=normu;
		cout<<" InternForceVector "<<norm(InternalForceVector , fem->NormType , K,ActiveDofs)<<" Load "<<norm(Temp , fem->NormType , K, ActiveDofs)<<"\n";flush(cout);
		if (bet<=fem->ErrorTolerance )
		{
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
		if (fem->FullNR == true)
		{
			K= fem->AssembleStiffnessMatrix();
			fem->LagrangeStiffnessMatrix(K);
			Kinv = K.i();
		}
		if ( (fem->EnableArcLengthMethod == true) && (fem->FullNR == true) )
		{
			DT = Kinv*ExternalLoad;
			fem->CalibrateVector(DT,true,false);
		}
		DTbar = Kinv*Residual;

/*			DEBUGING STUF
		printmatrix("kinv",Kinv,1);
		printvector("Residual",Residual,1,ActiveDofs);
		printvector("K*U",K*U,1,ActiveDofs);
		printvector("Load",Temp,1,ActiveDofs);
		printvector("InternalForceVector",InternalForceVector,1,ActiveDofs);
		printvector("DTbar",DTbar,1,DTbar.Nrows());
		printvector("DT",DT,1,ActiveDofs);
		printvector("Uprin",U,1,ActiveDofs);
*/
		fem->CalibrateVector(DTbar,true,false);
		if (fem->EnableArcLengthMethod ==true)
		{
			ArcLengthRoutine(fem,DT,DTbar,totalloadfactor,U,UPreviousIteration,ExternalLoad,illfail,iterations);
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
			U = U + DTbar;
			fem->CalibrateVector(U,true,false);
			fem->UpdateNodes(U.Rows(1,ActiveDofs) );
			fem->AssembleInternalLoads(true,fem->Method);
			if (fem->MaxNoLineSearches > 0)
			{
				LineSearchLoop(fem,U,DTbar,Residual,totalloadfactor,ExternalLoad,K,illfail,Slol);
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


void PrepareNextIncr(FEM *fem,double &totalloadfactor,double &incrloadfactor,int iterations,double &oldtotalloadfactor,double bet,ColumnVector &U,ColumnVector &UPrevious,double DL)
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
		U=UPrevious ;
		fem->UpdateNodes(U.Rows(1,fem->GetNumberOfActiveSystemDofs()) );
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
		fem->DesiredLengthIncrement= ReductionFactor*DL;
		fem->DesiredLengthIncrement= Min(fem->DesiredLengthIncrement,fem->MaxLengthIncrement);
		fem->DesiredLengthIncrement= Max(fem->DesiredLengthIncrement,fem->MinLengthIncrement);
	}
	return;
}




void FEM::ArcLengthMethod()
{
	FEM *fem;
	fem=this;//Not good programming but i need this to enable intellisense whtich does not work else
	StiffnessMatrix K;
	StiffnessMatrix Kinv; // Stiffness Matrix Vector
	CurrentLoadStep=0;
	double bet;
	double TotalLoadFactor=0;
	double OldTotalLoadFactor=0;
	double CatholicLoadFactor=0;
	double Stift;//This is <DT,OutOfBalanceForcesVector>
	double Stif;//Stiffness Parameter
	double Stifi;//Initial Stiffness Parameter
	double CStif;//Current Stiffness Parameter
	double StrainEnergy=0;
	int neg=0;//number of negative pivot points of K in an LDL analysis
	int ActiveSystemDofs = fem->GetNumberOfActiveSystemDofs();
	double DL;//this is <DT,DT>
	ColumnVector L; //Load Vector (Permenant)
	int VectorLength = GetNumberOfActiveSystemDofs();//μέγεθος του διανυσματος L
	ColumnVector LPS(VectorLength); //Load Vector (Previous Stages)
	ColumnVector LCS(VectorLength); //Load Vector (Current Stages)
	ColumnVector U(VectorLength); //Displacement Vector
	ColumnVector UPrevious(VectorLength); //Previous Displacement Vector
	ColumnVector LNext(VectorLength);//Φορτίο που προσπαθώ να φτάσω στο επόμενο βήμα ισορροπίας
	ColumnVector LPrevious(VectorLength);//Φορτίο στο προηγούμενο βήμα ισορροπίας ή στο τρέχον φορτίο
	ColumnVector DT(VectorLength); //Displacement Vector from Total Load Vector = Kinv*LoadVector
	//that is DT is Tangential Displacements due to external applied loads
	ColumnVector Residual(VectorLength); //out of balance forces
	LowerTriangularMatrix Lower(VectorLength) ;
	DiagonalMatrix Diagonal(VectorLength);	
	ColumnVector Temp=0;
	ColumnVector DU=0;

	Matrix Tempp;
	int iterations=0; //iterations made by iterateto equilibrium function. Not Load Steps.Load Steps = incr
	// Are we STARTING OR RESTARTING????//
	if (Restart== false)
	{
		TotalLoadFactor=0;
		U=0;
		UpdateNodes(U);
		CurrentStage=1;
	}
	else
	{
		ReadRestartParameters(Stifi,StrainEnergy,TotalLoadFactor,CurrentStage,Restart,VectorLength,U,DisplacementAtPreviousStage);//These are Current StiffnessParameter,TotalLoadFactor,Displacements
		CurrentLoadStep=Restart;
	}
	////////////////Write Headers To Stiffness File/////////////
	fstream Stiffness;
	Stiffness.open(this->CurrentStiffnessFile,ios::app);//Εξωτερικά φορτία
	Stiffness<<"CurrStif\t\t"<<"NegEntries\t"<<"LFactor\t\t\t"<<"StrEner\t\t\t"<<"LStep\t\tArcUsed\t\tArcHint\n";
	flush(Stiffness);
	Stiffness.close();

	for(int l=CurrentStage;l<=NumberOfStages;l++)
	{
		CurrentStage=l;
		if(CurrentStage==1)
		{
			LPS=0;
		}
		else
		{
			LPS = this->AssembleStageLoads(1,CurrentStage-1)+ this->AssembleLoads();
		}
		
		if(CurrentStage==1)
		{
			LCS = this->AssembleLoads()+this->AssembleStageLoads(1);
		}
		else
		{
			LCS = this->AssembleStageLoads(CurrentStage);
		}


		LagrangeLoadVector(LCS);
		ColumnVector Zeros(countfixes()+countkinimatic(CurrentStage));Zeros=0;
		LPS &= Zeros; 
		U &=Zeros;
		if(Restart==0)
		{
			TotalLoadFactor=0;
			OldTotalLoadFactor=0;
			DisplacementAtPreviousStage=U;
		}
		else
		{
			Restart=0;//Next time it will be a normal iteration.....
		}
		


		// MAIN ROUTINE//
		for(int incr=1;incr<=fem->IterationsDesired;incr++)
		{
			if( dot(UPrevious,U,ActiveSystemDofs)<0 ){neg=1;}
			UPrevious = U; //Saving Previous displacement Vector
			DT=LCS;
			OldTotalLoadFactor=TotalLoadFactor;//Saving old Force for cut increments and old total loadlevel

			//Calculate Stiffness Matrix And Solve for U
			fem->UpdateNodes(U.Rows(1,ActiveSystemDofs));
			fem->AssembleInternalLoads(true,fem->Method);
			K= this->AssembleStiffnessMatrix();
			//Impose boundary conditions//
			//Penalty(K);
			//ldl(K,Lower,Diagonal,neg,VectorLength);
			LagrangeStiffnessMatrix(K);
			Residual=DT;
			if (SolverType==1) { Kinv=K.i();}
			switch (SolverType)
			{
			case 1 : //Gauss Elimination
				DT=Kinv*LCS;
				break;
			case 2:
				//κωδικας για αλλους solvers
				break;
			}//end switch
			
			//Now compute the number of negative pivots and parameters 
			//for current stiffness parameter and length Increment
			Stift = 0;
			DL = 0;
			CalibrateVector(Residual,true,false);//Zero the fixed dof coords
			CalibrateVector(DT,true,false);
			Stift = dot(DT,Residual,ActiveSystemDofs);
			DL=dot(DT,DT,ActiveSystemDofs);
			//Now Crisfield Suggest Analysing K in LDL format
			//To determin if K is Positive Definite
			//By Counting the non negative entries in D
			//I avoid doing it here
			//Because i want to use iterative solvers too.
			//So an LDL analysis is too much work.
			//Anyway i ll do the Ldl analysis for now
			Stif=Stift/DL;
			if(CurrentLoadStep==0 && Restart==false)
			{
				Stifi=Stif;
			}
			CStif=Stif/Stifi;
			DL=sqrt(DL);
			if (neg>1)
			{
/*				cout<<"negative pivots "<<neg;flush(cout);
				DiagonalMatrix EM(GetNumberOfActiveSystemDofs());
				EigenValues(K.SubMatrix(1,1,GetNumberOfActiveSystemDofs(),GetNumberOfActiveSystemDofs()),EM);
				fstream eigen;
				eigen.open("eigenvalues.txt",ios::app);//Εξωτερικά φορτία
				eigen.setf(ios::showpoint);eigen.precision(7);
				printvector(eigen,EM,CurrentLoadStep);
				eigen.flush();
				eigen.close();
//				exit(0);*/
			}
			//Now We will call the scaleup routine
			//Its main function is to compute Δλ and Δl
			//and update λ and displacements p calculated in the predictor solution
			scaleup(fem,IncrLoadFactor,TotalLoadFactor,DL,CStif,neg,U,UPrevious,DT);
			//Now Its time to iterate utill we reach equilibrium
			//using the iteratetoequilibrium function
			iterations=IterateToEquilibrium(fem,K,Kinv,IncrLoadFactor,TotalLoadFactor,DL,U,UPrevious,DT,LCS,LPS,bet);
			//Write results to file
			AssembleInternalLoads(true,Method);//Update stress condition in structure
			//Γράφω τα στοιχεία στο αρχείο αν έχει επιτευχθεί σύγκλιση//
			if(bet<fem->ErrorTolerance)
			{
				CurrentLoadStep++;
				DU = U-UPrevious;
				Temp = K*DU;
				StrainEnergy += dot(DU,Temp,ActiveSystemDofs );
				//Write Info To File
				Stiffness.open(this->CurrentStiffnessFile,ios::app);//Εξωτερικά φορτία
				Stiffness.setf(ios::showpoint);Stiffness.precision(7);
				Stiffness<<CStif<<"\t\t"<<neg<<"\t\t"<<TotalLoadFactor<<"\t\t"<<StrainEnergy<<"\t\t"<<CurrentLoadStep <<"\t\t"<<fem->DesiredLengthIncrement<<"\t\t";
				flush(Stiffness);
				Stiffness.close();
				//Written Info To file
				//Write only when a converge was reached
				ColumnVector TL;//TotalLoad
				TL=(TotalLoadFactor*LCS)+LPS;				
//				ColumnVector tt=AssembleInternalLoads(true,Method,true);
//				loads<<(norm(tt,2,K))/norm(TL,2,K)<<"\n";
				fstream loadfile;loadfile.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
				printvector(loadfile,TL,CurrentLoadStep,ActiveSystemDofs);
				loadfile.close();
//				ColumnVector tt=AssembleInternalLoads(true,Method,true)-TL;
//				loads<<(norm(tt,2,K))/norm(TL,2,K)<<"\n";
				loadfile.close();
				fstream displ;displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
				printvector(displ,U,CurrentLoadStep,GetNumberOfActiveSystemDofs());
				displ.close();
				this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
				this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
				this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
				this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
				WriteToRestartFile(Stifi,StrainEnergy,TotalLoadFactor,CurrentStage,CurrentLoadStep,VectorLength,U);
				//end of writing to file
			}
			if (AutoIncrement==true)
			{
				PrepareNextIncr(fem,TotalLoadFactor,IncrLoadFactor,iterations,OldTotalLoadFactor,bet,U,UPrevious,DL);
			}
			if(TotalLoadFactor>=0.99999999999)
			{break;}
		}//end for
	}//end for stages
	//We are over and we write info to restart file
	//so we can start over from a previous step	
}










