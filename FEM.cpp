// FEM.cpp: implementation of the FEM class.
//
//////////////////////////////////////////////////////////////////////

#include "FEM.h"
#include "Statheres.h"
#include <fstream.h>
#include "iostream.h" // for testing
#include "Newmatio.h" // for testing

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


FEM::FEM()
{
	CCenter = new ControlCenter(this);// Initialise the control Center
	ConstrainCenter = new Constrains(this); 
		ConstrainCenter->SetConstrainMethod(PENALTY_METHOD);
//		ConstrainCenter->SetConstrainMethod(LAGRANGE_MULTIPLIERS);
	PCenter = new ParametersCenter(this);
	geometry = 2;
	variables = 3;
	PenaltyNumber = 0;
	FixCommandsNumber=0;
	NumberofItegrationLoops=1;
	MaxIterations=10;
	ElementFileOut="Results.txt";
	GeometryFileOut="Stoixeia.txt";
	DisplacementsFileOut="Displacements.txt";
	ReactionFileOut="Reactions.txt";
	ForceFileOut="ForceFile.txt";
	LoadVectorFileOut="Loads.txt";//Εδω εχω τα φορτία σε κάθε επαναληπτικό βήμα
	CurrentStiffnessFile="CurrentStiffnessFile.txt";
	GeneralDataFileName="GeneralData.txt";
	InfoFile="InfoFile.txt";
	IterationsDesired=10;
	ValidateStress=true;
	SolverType=1;
	SolvingMethod=0;//Elastic Solver is default
	Method=2;
	LoadIncreament=10;
	ErrorTolerance=0.01;
	ErrorDispTolerance=0.01;
	InitialDeltaL=0.10;// Initial Value of Δλ
	WriteExtraInfo=true;//Write extended information in file
	AutoIncrement=true;//Enable autoincrementation
	EnableArcLengthMethod=false;//Enable ArcLength
	Accelerate=false; //Enable secant type accelaration methods (such as Broyden method)
	Restart=0;//are we restarting??
	MaxNoLineSearches=0;
	MaxResidualCalculations=5;
	FullNR=true;
	DisplacementControl = false;
	NormType = 1;
	ArcLengthAlgorithm = 1;
	DesiredNumberOfIterations= 10;
	DesiredLengthIncrement= 0.0;
	IncrLoadFactor= 0.1;
	InitialDeltaL =1.0;
	LowestStiffnessSwitch = 0.3;
	MaxAmpAtAnyStemp =5.0;
	MaxResidualCalculations  = 10;
	MaxLengthIncrement =0.5;
	MinLengthIncrement =0.5;
	ToleranceOnRatio =0.8;
	Restart= false;
	MaxLoadIncrement= 0.01;
	MinLoadIncrement= 0.0001;
	MaxTotalStepLength=25;
	MinTotalStepLength=0.01;
	ArcHint=0;
	CurrentStage = 1;
	NumberOfStages = 1;
	NumberOfStageCommands=0;
	CurrentLoadStep= 0;
	PenaltyScaleFactor=8;
	IterationsDesired=10;
	//Line Search stuff//



}

FEM::~FEM()
{

}

void FEM::WriteFileheaders()
{
	fstream kos;
	kos.open(this->DisplacementsFileOut,ios::app);
	kos<<"   LS    Dof         Disp\n";
	kos.open(this->ReactionFileOut,ios::app);
	kos<<"   LS    Dof         Load\n";
	kos.close();
	kos.open(this->ForceFileOut,ios::app);
	kos<<"   LS    Elem   Dof       Load\n";
	kos.close();
	kos.open(this->ElementFileOut,ios::app);
	kos<<"   Load  Elem    e/s      gp    comp       Value\n";
	kos.close();
	kos.open(this->LoadVectorFileOut,ios::app);
	kos<<"   LS    Dof         Load\n";

}



StiffnessMatrix FEM::AssembleStiffnessMatrix()
{
	int matrixlength;
	matrixlength=GetNumberOfActiveSystemDofs();
//	StiffnessMatrix  Stiffness(matrixlength,band,band);
	StiffnessMatrix  Stiffness(matrixlength);
	Stiffness=0;
	for(int i=1;i<=this->numberofelements;i++)
	{
		Matrix A=elements[i-1]->CreateStiffnessMatrix();
		for(int j=1;j<=A.Nrows();j++)
		{
			for(int k=1;k<=A.Nrows();k++)
			{
				//βρισκω την θέση που θα βαλω το στοιχείο A(j,k)
				Stiffness(this->elements[i-1]->ConnectivityVector()(j),this->elements[i-1]->ConnectivityVector()(k))+=A(j,k);
			}//enf for(k)
		}//end for(j)
	}//end for(i)
	return Stiffness;
}

int FEM::systemband()
{
	int bandmax=0;
	int temp;
	for(int i=1;i<=this->numberofelements;i++)
	{
		temp=this->elements[i-1]->band();
		if(bandmax<temp)
		{
			bandmax=temp;
		}
	}
	return bandmax;
}

FEM::Initialize()
{
	WriteFileheaders();
	for(int i=1;i<=numberofelements;i++)
	{
		elements[i-1]->ActivateExistingDofs();
	}
//	band=systemband();
	LostDofsVector=LostDofs();
	this->InitialiseGaussLibrary();
	CalculatePenalty();
	int VectorLength = this->numberofNodes*this->GetVariables()-LostDofsVector(this->numberofNodes+1);//μέγεθος του διανυσματος L
	DisplacementAtPreviousStage.ReSize(VectorLength);
	DisplacementAtPreviousStage=0;

}

ColumnVector FEM::LostDofs()
{
	ColumnVector Ldofs(this->numberofNodes+1);
	Ldofs(1)=0;
	int sum=0;
	//Για καθε κόμβο της κατασκευής βλέπω πόσοι απο τους συνολικούς βαθμούς
	//ελευθεριάς που ειναι (DIM-1)*3 δεν υπάρχουνε απο τον κομβο 1 εως και τον n-1
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			//Ενδιαφερομαι για τον κόμβο n.Αρα πρέπει να ψαξω τους χαμένους βαθμούς ελευθεριας
			//Μέχρι τον κομβο n-1.Αρα αν ειναι n = i+1 τοτε n-1= i και λόγω τις
			//C++ αυτος φιλασεται στιν θέση nodes[i-1]
			if(this->nodes[i-1]->GetExist(j)==0)
			{
				sum+=1;
			}
		}
		Ldofs(i+1)=sum;
	}
	//Το τελευταίο στοιχείο του Ldofs αν αφαιρεθεί από 
	//το (DIM-1)*3 ειναι ταυτόχρονα και μέγεθος του πίνακα ακαμψίας
	return Ldofs;

}

ColumnVector FEM::AssembleLoads()
{
	int Vectorlength;
	Vectorlength=GetNumberOfActiveSystemDofs();
	ColumnVector loads(Vectorlength);
	int tempindex=0;
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			if(nodes[i-1]->GetExist(j)==1)
			{
				tempindex+=1;
				loads(tempindex)=nodes[i-1]->GetLoad(j);
			}//endif
		}//end for(j)
	}//end for(i)
	//////////Now I impose Penalty vector Loads.....///////////// (Only when using Penalty method)
/*	double BIG = GetPenalty();
	for(i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= this->Fixes[i-1];
		int vecsize = TempCommand->coeff.size();
		if(CurrentStage==1 && CurrentStage>=TempCommand->startingstage && CurrentStage<=TempCommand->endingstage)//Process Fix Command
		{
			double pivot=0;
			for(int l=1;l<=TempCommand->shiftingtimes;l++)
			{
				int svalue= TempCommand->shiftvalue;
				for(int j=1;j<=vecsize;j++)
				{
					pivot = TempCommand->coeff[j-1]*TempCommand->result*BIG;//calculate constant to add to Load...
					int dofi = GetGlobalDof(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);

					loads(dofi)+=pivot;
				}
			}
		}
	}//Next Command
*/
	return loads;
}

ColumnVector FEM::AssembleStageLoads(int fromstage,int tostage)
{

	int Vectorlength;
	Vectorlength=GetNumberOfActiveSystemDofs();
	ColumnVector loads(Vectorlength);
	loads=0;
	for(int i=1;i<=this->NumberOfStageCommands;i++)
	{
		if(Stages[i-1]->tostage<=tostage && Stages[i-1]->keep==0 && fromstage<=Stages[i-1]->tostage)
		{
		}
		else if(Stages[i-1]->tostage<fromstage)
		{
		}
		else if (tostage<Stages[i-1]->fromstage)
		{
		}
		else
		{
			int snode = Stages[i-1]->fromnode;
			int enode = Stages[i-1]->tonode;
			int step =  Stages[i-1]->step;
			int dof = Stages[i-1]->dof;
			for(int j=snode;j<=enode;j+=step)
			{
				int vecposition = GetGlobalDof(j,dof);
				loads(vecposition)+=Stages[i-1]->load;
			}
		}
	}// end for
	return loads;
}

ColumnVector FEM::AssembleStageLoads(int stage)
{

	int Vectorlength;
	Vectorlength= GetNumberOfActiveSystemDofs();
	ColumnVector loads(Vectorlength);
	loads=0;
	for(int i=1;i<=this->NumberOfStageCommands;i++)
	{
		if(Stages[i-1]->fromstage==stage)
		{
			int snode = Stages[i-1]->fromnode;
			int enode = Stages[i-1]->tonode;
			int step =  Stages[i-1]->step;
			int dof = Stages[i-1]->dof;
			for(int j=snode;j<=enode;j+=step)
			{
				int vecposition = GetGlobalDof(j,dof);
				loads(vecposition)+=Stages[i-1]->load;
			}
		}
	}// end for
	return loads;
}

ColumnVector FEM::AssemblePenaltyStageLoads(int stage)
{
	int Vectorlength;
	Vectorlength=GetNumberOfActiveSystemDofs();
	ColumnVector loads(Vectorlength);
	int tempindex=0;
	loads=0;
	//////////Now I impose Penalty vector Loads...../////////////
	double BIG = GetPenalty();
	for(int i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= this->Fixes[i-1];
		int vecsize = TempCommand->coeff.size();
		if(CurrentStage>=TempCommand->startingstage && CurrentStage<=TempCommand->endingstage)//Process Fix Command
		{
			double pivot=0;
			for(int l=1;l<=TempCommand->shiftingtimes;l++)
			{
				int svalue= TempCommand->shiftvalue;
				for(int j=1;j<=vecsize;j++)
				{
					pivot = TempCommand->coeff[j-1]*TempCommand->result*BIG;//calculate constant to add to Load...
					int dofi = GetGlobalDof(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);

					loads(dofi)+=pivot;
				}
			}
		}
	}//Next Command

	return loads;
}



FEM::ElasticSolve(int solvertype)
{
	ElasticSolver el(this);
	el.Solve();
/*	StiffnessMatrix K;  // Stiffness Matrix Vector
	ColumnVector L; //Load Vector
	ColumnVector U; //Displacement Vector
	K=this->AssembleStiffnessMatrix();
	L=this->AssembleLoads()+ AssembleStageLoads(1);
	int n=L.Nrows();// this should be the same as ""GetNumberOfActiveSystemDofs()""
//	Penalty(K);// THis imposes boundary conditions by using penalty method
	LagrangeStiffnessMatrix(K);
	LagrangeLoadVector(L);
	if (solvertype == 1)
	{
		U=GaussSolver(K,L);
	}
	UpdateNodes(U.Rows(1,n));// Update Only the active dofs Ignoring Lagrange Multipliers
	//////////////Αποθήκευση των αποτελεσμάτων σε αρχεία/////////////
	fstream loads;
	loads.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
	printvector(loads,L,1,GetNumberOfActiveSystemDofs());
	loads.close();
	fstream displ;
	displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
	printvector(displ,U,1,GetNumberOfActiveSystemDofs());
	displ.close();
	this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
	this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
	this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
	this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις*/

}

void FEM::CalculatePenalty()
{
	StiffnessMatrix K;  // Stiffness Matrix Vector
	K=this->AssembleStiffnessMatrix();
	SetPenalty(K.MaximumAbsoluteValue()*pow(10,GetPenaltyScaleFactor() ));
	if(GetPenaltyScaleFactor()==0){SetPenalty(0.0);}
}

FEM::Penalty(StiffnessMatrix& pinakas)
{
	int tempindex=0;
	double BIG = GetPenalty();
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			if(nodes[i-1]->GetExist(j)==1)
			{
				tempindex+=1;
				if(nodes[i-1]->GetFixed(j)==1)
				{
					pinakas(tempindex,tempindex)+=BIG;
				}//endif
			}//endif
		}//end for(j)
	}//end for(i)
	for(i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= this->Fixes[i-1];
		int vecsize = TempCommand->coeff.size();
		if(CurrentStage>=TempCommand->startingstage && CurrentStage<=TempCommand->endingstage)//Process Fix Command
		{
			double pivot=0;
			for(int l=1;l<=TempCommand->shiftingtimes;l++)
			{
				int svalue= TempCommand->shiftvalue;
				for(int j=1;j<=vecsize;j++)
				{
					for(int k=1;k<=vecsize;k++)
					{
						pivot = TempCommand->coeff[j-1]*TempCommand->coeff[k-1]*BIG;//calculate constant to add to Matrix...
						int dofi = GetGlobalDof(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);
						int dofj = GetGlobalDof(TempCommand->dofs[2*(k-1)]+svalue*(l-1),TempCommand->dofs[2*(k-1)+1]);
						pinakas(dofi,dofj)+=pivot;
					}
				}
			}
		}
	}//Next Command

}

FEM::UpdateNodes(const ColumnVector& disp)
{
	int tempindex=0;
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			if(nodes[i-1]->GetExist(j)==1)
			{
				tempindex+=1;
				if(nodes[i-1]->GetFixed(j)==1)
				{
					nodes[i-1]->SetDisplacement(j,0);
				}
				else
				{
					nodes[i-1]->SetDisplacement(j,disp(tempindex));
				}//endif
			}//endif
		}//end for(j)
	}//end for(i)

}

ColumnVector FEM::ReactionAtNode(int nodenumber)
{

	ColumnVector Temp;
	Temp.ReSize(GetVariables());
	Temp=0;
	int elementnumber=0;//το στοιχείο υπο επεξεργασία
	int nodenumberinelement=0;
	for(int i=1;i<=this->nodes[nodenumber-1]->GetNumberOfElements();i++)
	{
		int tempindex=0;
		elementnumber=nodes[nodenumber-1]->ElementGet(i);//επιστρέφει τον αριθμό του στοιχείου στον κομβο
		//ColumnVector V=(elements[elementnumber-1]->CreateStiffnessMatrix())*(elements[elementnumber-1]->GetDisplacementOfNodes());
		ColumnVector V=elements[elementnumber-1]->InternalForce(false,Method);
		//Βρες ποιος κομβος του στοιχείου ειναι ο προς επεξεργασία κόμβος
		for(int k=1;k<=elements[elementnumber-1]->GetNumberOfNodes();k++)
		{
			if(elements[elementnumber-1]->GetNode(k)==nodenumber)
			{
				nodenumberinelement=k;
			}//endif
		}//endfor k
		for(int j =1;j<=GetVariables();j++)
		{
			if((this->elements[elementnumber-1]->Getdofcode(j))!=0)
			{
				tempindex++;
				Temp(j)+=V((nodenumberinelement-1)*elements[elementnumber-1]->dofs()+tempindex);
			}//endif
		}//endfor j
	}//endfor i
	return Temp;

}

void FEM::NonLinearSolve()
{
	switch (this->SolvingMethod)
	{
		case 1: 
			TangentStiffnessN_RComplete();
			break;
		case 2:
			TangentStiffnessN_RIncomplete();
			break;
		case 3:
			SecantStiffnessMethod();
			break;
		case 4:
			OrthogonalResidualMethod();
			break;
		case 5:
			SecantStiffnessStepMethod();
			break;
		case 6:
//			ArcLengthMethod();
			ArcLength arc(this);
			arc.Solve();
			break;
	}
}


ColumnVector FEM::AssembleInternalLoads(bool update,int method,bool withpenaltyloads)
{
	//Την Ρουτίνα αυτη την χρησιμοποιώ για να φτίαξω το διάνυσμα των εσωτερικών δυνάμεων 
	//που αναπτύσονται στην κατασκευή και να δώ αν βρίσκεται σε ισορρόπία με αυτό τον επιβαλλόμενων φορτίων
	int Vectorlength;
	Vectorlength=(GetVariables()*this->numberofNodes)-(this->LostDofsVector(this->numberofNodes+1));
	ColumnVector loads(Vectorlength);
	loads=0;
	for(int i=1;i<=this->numberofelements;i++)
	{
		int numberofnodes=elements[i-1]->dofs()*elements[i-1]->GetNumberOfNodes();
		ColumnVector elementforce=elements[i-1]->InternalForce(update,method);
		for(int j=1;j<=numberofnodes;j++)
		{
			loads( elements[i-1]->ConnectivityVector()(j) )+=elementforce(j);
		}//end for j
	}//end for i
	// Τωρα αφαιρώ απο το ανω διανυσμα τις αντιδράσεις των στηρίξεων
	//Ετσι ώστε να μπορώ να το συγκρινω με το διανυσμα των φορτίων
	int tempindex=0;//βοηθητικός μετρητης
	for(int k=1;k<=this->numberofNodes;k++)
	{
		for(int l=1;l<=GetVariables();l++)
		{
			if(nodes[k-1]->GetExist(l)==1)
			{
				tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
				if(nodes[k-1]->GetFixed(l)==1) //Ο βαθμός ελευθερίας ειναι δεσμευμένος
				{
					loads(tempindex)=nodes[k-1]->GetLoad(l);
				}//endif
			}//endif
		}//endfor l
	} //endfor k

/*	//////////Now Calculate the penalty loads imposed previously.../////////////
	if(withpenaltyloads==true)
	{
		double BIG = GetPenalty();
		for(i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
		{
			FixCommand* TempCommand= this->Fixes[i-1];
			if(CurrentStage>=TempCommand->startingstage && CurrentStage<=TempCommand->endingstage)//Process Fix Command
			{
				int vecsize = TempCommand->coeff.size();
				double pivot=0;
				for(int l=1;l<=TempCommand->shiftingtimes;l++)
				{
					int svalue= TempCommand->shiftvalue;
					for(int j=1;j<=vecsize;j++)
					{
						pivot = BIG*TempCommand->coeff[j-1];
						double tempsum=0;
						for(int k=1;k<=vecsize;k++)
						{
							int nodetoadd=TempCommand->dofs[2*(k-1)]+svalue*(l-1);
							int doftoadd=TempCommand->dofs[2*(k-1)+1];
							double dispfrompreviousstage=DisplacementAtPreviousStage(GetGlobalDof(nodetoadd,doftoadd));
							tempsum+= TempCommand->coeff[k-1]*nodes[nodetoadd-1]->GetDisplacement(doftoadd);
							tempsum-= TempCommand->coeff[k-1]*dispfrompreviousstage;
							//calculate constant to add to Load...
						}
						int dofi = GetGlobalDof(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);

						loads(dofi)+=pivot*(tempsum);
					}
				}
			}
		}//Next Command
	}

	int m,n;
	m=countfixes();n=countkinimatic(CurrentStage);
	ColumnVector Zeros(n+m);Zeros=0;
	loads=loads & Zeros;

	
	if(withpenaltyloads==true)
	{
		double BIG = GetPenalty();
		for(i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
		{
		FixCommand* TempCommand= this->Fixes[i-1];
		int vecsize = TempCommand->coeff.size();
		double pivot=0;
		for(int l=1;l<=TempCommand->shiftingtimes;l++)
		{
			int svalue= TempCommand->shiftvalue;
			for(int j=1;j<=vecsize;j++)
			{
//				double currdisp = nodes[TempCommand->dofs[2*(j-1)]+svalue*(l-1)-1]->GetDisplacement(TempCommand->dofs[2*(j-1)+1]);
				pivot = TempCommand->coeff[j-1]*LoadIncr*BIG*TempCommand->result;//calculate constant to add to Load...
				int dofi = dofinSystem(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);

				loads(dofi)+=pivot;
				}
			}
		}//Next Command
	}*/


	return loads;
}

void FEM::TangentStiffnessN_RComplete()
{
	//Solve Method = 1 FULL NEWTON-RAWSHON = Update Stiffness matrix at every iteration
	//Solver Type = 1 ----> GaussElimination
	CurrentLoadStep=0; //for testing
	StiffnessMatrix K;  // Stiffness Matrix Vector
	ColumnVector L; //Load Vector
	L=this->AssembleLoads();
	int VectorLength = L.Nrows();//μέγεθος του διανυσματος L
	ColumnVector LPrevious(VectorLength);//Φορτίο στο προηγούμενο βήμα ισορροπίας ή στο τρέχον φορτίο
	LPrevious=0;
	ColumnVector LNext(VectorLength);//Φορτίο που προσπαθώ να φτάσω στο επόμενο βήμα ισορροπίας
	LNext=0;
	ColumnVector U(VectorLength); //Displacement Vector
	U=0;
	ColumnVector DU(VectorLength);
	DU=0;
	ColumnVector LoadImbalance(VectorLength);
	//Θέτω τις τεταγμένες που αντιστοιχούν σε fixαρισμένους βαθμούς ελευθερίας
	//ίσες με το φορτίο που ασκείται σε αυτές γιατί δεν παίζει ρόλο η τμηματική αυξηση του φορτίου
	// σε αυτές καθώς η όποια δύναμη μεταβιβάζεται κατευθείαν στην στήριξη
	int tempindex=0;//βοηθητικός μετρητης
	double icrpercent = this->IncrLoadFactor;//Pososto auksisis fortioy
	int SolverType = this->SolverType;
	double ErrorTolerance = this->ErrorTolerance;

	for(int k=1;k<=this->numberofNodes;k++)
	{
		for(int l=1;l<=GetVariables();l++)
		{
			if(nodes[k-1]->GetExist(l)==1)
			{
				tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
				if(nodes[k-1]->GetFixed(l)==1) //Ο βαθμός ελευθερίας ειναι δεσμευμένος
				{
					LPrevious(tempindex)=nodes[k-1]->GetLoad(l);
					LNext(tempindex)=nodes[k-1]->GetLoad(l);
				}//endif
			}//endif
		}//endfor l
	} //endfor k
	while ( (L-LNext).Norm1() >=0.00001 && CurrentLoadStep<MaxIterations && icrpercent*CurrentLoadStep<=1) //Αν δεν εχω φτασει ακόμη το φορτίο που επιβάλεται
	{
		CurrentLoadStep++;
		for (int i=1;i<=VectorLength;i++)//Ευρεση του φορτίου ισορροπίας
		{
			LNext(i)+=icrpercent*L(i);
	
		}//endfor i  Βρήκα το νέο φορτίο ισορροπίας
		LoadImbalance=LNext-LPrevious;
		if(CurrentLoadStep>1)
		{
			fstream loads;
			loads.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
			loads<<LPrevious.t();
			printvector(loads,LPrevious,CurrentLoadStep);
			flush(loads);
			loads.close();
			fstream displ;
			displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
			printvector(displ,U,CurrentLoadStep);
			flush(displ);
			displ.close();
			this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
			this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
			this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
			this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
		}//Μετα την πρώτη βόλτα που τα στοιχεία εχουν αποκτήσει έντάσεις
		//Γράφω στο αρχείο τα στοιχεία κάθε στοιχείου
		while ( LoadImbalance.Norm1()>ErrorTolerance )
		{
			K= this->AssembleStiffnessMatrix();
			Penalty(K);
			switch (SolverType)
			{
			case 1 : //Gauss Elimination
				DU=GaussSolver(K,LoadImbalance);
				break;
			case 2:
				//κωδικας για αλλους solvers
				break;
			}//end switch
			U+=DU;
			this->UpdateNodes(U);
			LPrevious=AssembleInternalLoads(true,Method);
			LoadImbalance=LNext-LPrevious;
		}//end while
	}//end while

}

void FEM::TangentStiffnessN_RIncomplete()
{
	//Solve Method = 2 SEMI NEWTON-RAWSHON = Update Stiffness matrix at every equilibrum;
	//Solver Type = 1 ----> GaussElimination

	double icrpercent = this->IncrLoadFactor;//Pososto auksisis fortioy
	int SolverType = this->SolverType;
	double ErrorTolerance = this->ErrorTolerance;

	CurrentLoadStep=0; //for testing
	StiffnessMatrix K;  // Stiffness Matrix Vector
	Matrix Kinv;//O αντιστροφος του Stiffness Matrix
	ColumnVector L; //Load Vector
	L=this->AssembleLoads();
	int VectorLength = L.Nrows();//μέγεθος του διανυσματος L
	ColumnVector LPrevious(VectorLength);//Φορτίο στο προηγούμενο βήμα ισορροπίας ή στο τρέχον φορτίο
	LPrevious=0;
	ColumnVector LNext(VectorLength);//Φορτίο που προσπαθώ να φτάσω στο επόμενο βήμα ισορροπίας
	LNext=0;
	ColumnVector U(VectorLength); //Displacement Vector
	U=0;
	ColumnVector DU(VectorLength);
	DU=0;
	ColumnVector LoadImbalance(VectorLength);
	//Θέτω τις τεταγμένες που αντιστοιχούν σε fixαρισμένους βαθμούς ελευθερίας
	//ίσες με το φορτίο που ασκείται σε αυτές γιατί δεν παίζει ρόλο η τμηματική αυξηση του φορτίου
	// σε αυτές
	int tempindex=0;//βοηθητικός μετρητης
	for(int k=1;k<=this->numberofNodes;k++)
	{
		for(int l=1;l<=GetVariables();l++)
		{
			if(nodes[k-1]->GetExist(l)==1)
			{
				tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
				if(nodes[k-1]->GetFixed(l)==1) //Ο βαθμός ελευθερίας ειναι δεσμευμένος
				{
					LPrevious(tempindex)=nodes[k-1]->GetLoad(l);
					LNext(tempindex)=nodes[k-1]->GetLoad(l);
				}//endif
			}//endif
		}//endfor l
	} //endfor k
	while ( (L-LNext).Norm1() >=0.00001 && CurrentLoadStep<MaxIterations && icrpercent*CurrentLoadStep<=1) //Αν δεν εχω φτασει ακόμη το φορτίο που επιβάλεται
	{
		CurrentLoadStep++;
		for (int i=1;i<=VectorLength;i++)//Ευρεση του φορτίου ισορροπίας
		{
			LNext(i)+=icrpercent*L(i);
	
		}//endfor i  Βρήκα το νέο φορτίο ισορροπίας
		if(CurrentLoadStep>1)
		{
			fstream loads;
			loads.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
			printvector(loads,LPrevious,CurrentLoadStep);
			flush(loads);
			loads.close();
			fstream displ;
			displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
			printvector(displ,U,CurrentLoadStep);
			flush(displ);
			displ.close();
			this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
			this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
			this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
			this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
		}//Μετα την πρώτη βόλτα που τα στοιχεία εχουν αποκτήσει έντάσεις
		//Γράφω στο αρχείο τα στοιχεία κάθε στοιχείου
		LoadImbalance=LNext-LPrevious;
		K= this->AssembleStiffnessMatrix();
		Penalty(K);
		if( SolverType==1)
		{
			Kinv=K.i();
		}
		while ( LoadImbalance.Norm1()>ErrorTolerance )
		{
			switch (SolverType)
			{
			case 1 : //Gauss Elimination
				DU=Kinv*LoadImbalance;
				break;
			case 2:
				//κωδικας για αλλους solvers
				break;
			}//end switch
			U+=DU;
			this->UpdateNodes(U);
			LPrevious=AssembleInternalLoads(true,Method);
			LoadImbalance=LNext-LPrevious;
		}//end while
	}//end while

}


void FEM::SecantStiffnessStepMethod()
{
	//Solve Method = 5 SECANT STIFFNESS APROACH - COMPLETE PATH
	//Solver Type = 1 ----> GaussElimination
	double LoadFactorIncrement = this->IncrLoadFactor;//Pososto auksisis fortioy
	int SolverType = this->SolverType;
	double ErrorTolerance = this->ErrorTolerance;

	long iteration; //for testing
	CurrentLoadStep =0;//Βημα φόρτισης
	StiffnessMatrix K;  // Stiffness Matrix Vector
	ColumnVector L; //Load Vector
	ColumnVector CurrentLoad; //CurrentLoad = LoadFactorIncerement*Load
	L=this->AssembleLoads();
	int VectorLength = L.Nrows();//μέγεθος του διανυσματος L
	ColumnVector LoadImbalance(VectorLength);//Διαφορά Εξωτερικού φορτίου και εσωτερικών δυνάμεων
	ColumnVector U(VectorLength); //Displacement Vector
	ColumnVector DU(VectorLength); //Displacement Vector Previous Step
	U=0;
	UpdateNodes(U);
	K= this->AssembleStiffnessMatrix();
	Penalty(K);
	for(double i=LoadFactorIncrement;i<=1;i+=LoadFactorIncrement)
	{
		CurrentLoad = i*L;
		LoadImbalance=CurrentLoad;
		DU=CurrentLoad;
		CurrentLoadStep++;
		iteration = 0;
		while ( sqrt(LoadImbalance.SumSquare())/sqrt(CurrentLoad.SumSquare()) >=ErrorTolerance  && iteration<MaxIterations ) //Αν δεν εχω φτασει ακόμη το φορτίο που επιβάλεται
		{
			iteration++;
			cout<<"Secant Stiffness Method: LoadStep =  " <<CurrentLoadStep<<" Equilibrium Iteration = "<<iteration<<"\n";
			flush(cout);
			switch (SolverType)
			{
			case 1 : //Gauss Elimination
				DU=GaussSolver(K,LoadImbalance);
				break;
			case 2:
				//κωδικας για αλλους solvers
				break;
			}//end switch
			U=U+DU;
			this->UpdateNodes(U);//Ενημερώνω τους κόμβους
			K= this->AssembleStiffnessMatrix();//Υπολογίζω τον καινούργιο πίνακα δυσκαμψίας
			LoadImbalance = CurrentLoad - K*U;
			Penalty(K);
			AssembleInternalLoads(true,Method);
			// Τωρα αφαιρώ απο το ανω διανυσμα τις αντιδράσεις των στηρίξεων
			//Ετσι ώστε να μπορώ να το συγκρινω με το διανυσμα των φορτίων
			int tempindex=0;//βοηθητικός μετρητης
			for(int k=1;k<=this->numberofNodes;k++)
			{
				for(int l=1;l<=GetVariables();l++)
				{
					if(nodes[k-1]->GetExist(l)==1)
					{
						tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
						if(nodes[k-1]->GetFixed(l)==1) //Ο βαθμός ελευθερίας ειναι δεσμευμένος
						{
							LoadImbalance(tempindex)=LoadFactorIncrement*nodes[k-1]->GetLoad(l);
						}//endif
					}//endif
				}//endfor l
			} //endfor k
			cout<<" Error = "<<(sqrt(LoadImbalance.SumSquare())/sqrt(CurrentLoad.SumSquare()))<<"\n";
			cout<<" Work = "<<0.5*U.t()*K*U<<"\n";
		}//end while
		if (iteration>=MaxIterations)
		{
			cout<< "Exceeded MaxIterations. Please increase the MaxIteration Parameter and rerun the program"<<"\n";
			exit(1);
		}
		/////////////////////////////
		///Writing results to files//
		/////////////////////////////
		fstream loads;
		loads.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
		printvector(loads,CurrentLoad,CurrentLoadStep);
		flush(loads);
		loads.close();
		fstream displ;
		displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
		printvector(displ,U,CurrentLoadStep);
		flush(displ);
		displ.close();
		this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
		this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
		this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
		this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
	}//end for i


}


void FEM::SecantStiffnessMethod()
{
	//Solve Method = 3 SECANT STIFFNESS APROACH 
	//Solver Type = 1 ----> GaussElimination
	double LoadFactorIncrement = this->IncrLoadFactor;//Pososto auksisis fortioy
	int SolverType = this->SolverType;
	double ErrorTolerance = this->ErrorTolerance;

	CurrentLoadStep=0; //for testing
	StiffnessMatrix K;  // Stiffness Matrix Vector
	ColumnVector L; //Load Vector
	L=this->AssembleLoads();
	int VectorLength = L.Nrows();//μέγεθος του διανυσματος L
	ColumnVector LNext(VectorLength);//Φορτίο που προσπαθώ να φτάσω στο επόμενο βήμα ισορροπίας
	LNext=0;
	ColumnVector U(VectorLength); //Displacement Vector
	U=0;
	K= this->AssembleStiffnessMatrix();
	Penalty(K);
	while ( sqrt((L-LNext).SumSquare())/sqrt(L.SumSquare()) >=ErrorTolerance && CurrentLoadStep<MaxIterations ) //Αν δεν εχω φτασει ακόμη το φορτίο που επιβάλεται
	{
		CurrentLoadStep++;
		if(CurrentLoadStep>1)
		{
			cout<<"Secant Stiffness Method:  Equilibrium Iteration = "<<CurrentLoadStep<<"\n";
			cout<<" Error = "<<sqrt((L-LNext).SumSquare())/sqrt(L.SumSquare())<<"\n";
			flush(cout);
			fstream loads;
			loads.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
			printvector(loads,LNext,CurrentLoadStep);
			flush(loads);
			loads.close();
			fstream displ;
			displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
			printvector(displ,U,CurrentLoadStep);
			flush(displ);
			displ.close();
			this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
			this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
			this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
			this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
		}//Μετα την πρώτη βόλτα που τα στοιχεία εχουν αποκτήσει έντάσεις
		//Γράφω στο αρχείο τα στοιχεία κάθε στοιχείου
		switch (SolverType)
		{
		case 1 : //Gauss Elimination
			U=GaussSolver(K,L);
			break;
		case 2:
			//κωδικας για αλλους solvers
			break;
		}//end switch
		this->UpdateNodes(U);//Ενημερώνω τους κόμβους
		K= this->AssembleStiffnessMatrix();//Υπολογίζω τον καινούργιο πίνακα δυσκαμψίας
		LNext=K*U;//Υπόλογίζω το νέο φορτίο
		// Τωρα αφαιρώ απο το ανω διανυσμα τις αντιδράσεις των στηρίξεων
		//Ετσι ώστε να μπορώ να το συγκρινω με το διανυσμα των φορτίων
		int tempindex=0;//βοηθητικός μετρητης
		for(int k=1;k<=this->numberofNodes;k++)
		{
			for(int l=1;l<=GetVariables();l++)
			{
				if(nodes[k-1]->GetExist(l)==1)
				{
					tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
					if(nodes[k-1]->GetFixed(l)==1) //Ο βαθμός ελευθερίας ειναι δεσμευμένος
					{
						LNext(tempindex)=nodes[k-1]->GetLoad(l);
					}//endif
				}//endif
			}//endfor l
		} //endfor k
		Penalty(K);
	}//end while


}

void FEM::StressAndStainToFile(char *filename,bool principal)
{
	fstream elem;
	elem.open(filename,ios::app);
	for(int i=1;i<=this->numberofelements;i++)
	{
		elements[i-1]->StressAndStrainToFile(filename,principal,CurrentLoadStep);
	}
	elem.close();
}

void FEM::OrthogonalResidualMethod()
{
	//Solve Method = 4 ΟRTHOGONAL RESIDUAL METHOD 
	//Solver Type = 1 ----> GaussElimination
	long iteration=0; //for testing
	CurrentLoadStep=0;
	int iterationsdesired=this->IterationsDesired;
	float C=1.0;double Labs,Lmax,ro,ksi,psi=1;
	StiffnessMatrix K;
	Matrix Kinv; // Stiffness Matrix Vector
	ColumnVector L; //Load Vector
	L=this->AssembleLoads();
	int VectorLength = L.Nrows();//μέγεθος του διανυσματος L
	ColumnVector LNext(VectorLength);//Φορτίο που προσπαθώ να φτάσω στο επόμενο βήμα ισορροπίας
	ColumnVector LPrevious(VectorLength);//Φορτίο στο προηγούμενο βήμα ισορροπίας ή στο τρέχον φορτίο
	ColumnVector U(VectorLength); //Displacement Vector
	ColumnVector DU(VectorLength);
	ColumnVector DUOld(VectorLength);
	ColumnVector DUStart(VectorLength);
	ColumnVector LoadImbalance(VectorLength);//Η διαφορα εσωτερικών και εξωτερικών φορτίων = residual
	ColumnVector DQ(VectorLength);//Εσωτερικό φορτίο - Φορτίο Ισορροπίας
	ColumnVector DF(VectorLength);//Διαφορά βήματος φόρτισης
	ColumnVector Delta(VectorLength);//Διαφορά βήματος φόρτισης
	ColumnVector Temp(VectorLength);//Διαφορά βήματος φόρτισης
	bool NewStep;
	LNext=0;LPrevious=0;U=0;DU=0;
	double icrpercent = this->IncrLoadFactor;//Pososto auksisis fortioy
	double totalLoadpercent=0;
	int SolverType = this->SolverType;
	double ErrorTolerance = this->ErrorTolerance;
	//Θέτω τις τεταγμένες που αντιστοιχούν σε fixαρισμένους βαθμούς ελευθερίας
	//ίσες με το φορτίο που ασκείται σε αυτές γιατί δεν παίζει ρόλο η τμηματική αυξηση του φορτίου
	// σε αυτές
	int tempindex=0;//βοηθητικός μετρητης
	for(int k=1;k<=this->numberofNodes;k++)
	{
		for(int l=1;l<=GetVariables();l++)
		{
			if(nodes[k-1]->GetExist(l)==1)
			{
				tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
				if(nodes[k-1]->GetFixed(l)==1) //Ο βαθμός ελευθερίας ειναι δεσμευμένος
				{
					LPrevious(tempindex)=nodes[k-1]->GetLoad(l);
					LNext(tempindex)=nodes[k-1]->GetLoad(l);
				}//endif
			}//endif
		}//endfor l
	} //endfor k
	while ( (sqrt( (L-LNext).SumSquare() ) >=this->ErrorTolerance)  && CurrentLoadStep<= this->IterationsDesired  && totalLoadpercent<=1) //Αν δεν εχω φτασει ακόμη το φορτίο που επιβάλεται
	{
		CurrentLoadStep++;
		totalLoadpercent+=icrpercent;
		for (int i=1;i<=VectorLength;i++)//Ευρεση του φορτίου ισορροπίας
		{
			LNext(i)+=icrpercent*L(i);
	
		}//endfor i  Βρήκα το νέο φορτίο ισορροπίας
		DF=LNext-LPrevious;
		if (DF.Norm1()==0)
		{
			break;
		}
		NewStep=true;
		do
		{
			if (NewStep==true)
			{
				K= this->AssembleStiffnessMatrix();
				Penalty(K);
				if (SolverType==1) { Kinv=K.i();}
				DUOld=DU;
				switch (SolverType)
				{
				case 1 : //Gauss Elimination
					DU=Kinv*DF;
					break;
				case 2:
					//κωδικας για αλλους solvers
					break;
				}//end switch
				if (dot(DUOld,DU)<0)
				{
					DF=-DF;
					DU=-DU;
				}
				double l=sqrt(DU.SumSquare());//νόρμα DU
				if (CurrentLoadStep==1)
				{
					Labs=C*l;
					Lmax=Labs;
				}
				if(iteration<iterationsdesired)
				{
					Lmax=Min(Labs,2*Lmax);
				}
				ro=Min(1.0,Lmax/l);
				DUStart=DU;
				DU=ro*DU;
				Temp=U;
			}//END of New Step
			iteration=0;
			do ///Equilibrium iterations
			{
				iteration++;
				U=Temp+DU;
				UpdateNodes(U);
				DQ=AssembleInternalLoads(true,Method)-LPrevious;
//					fstream tt;
//					tt.open("tt.txt",ios::app);//Εξωτερικά φορτία
//   				printvector(tt,DQ,CurrentLoadStep);
//					printvector(tt,DF,CurrentLoadStep);

				ksi=dot(DQ,DU)/dot(DF,DU);
				LoadImbalance=-DQ+ksi*DF;
				switch (SolverType)
				{
				case 1 : //Gauss Elimination
					Delta=Kinv*LoadImbalance;
					break;
				case 2:
					//κωδικας για αλλους solvers
					break;
				}//end switch
//				double eta = -dot(DQ,Delta)/dot(DQ,DU);
//				Delta = (1/(1+eta))*Delta;
				DU=DU+Delta;
				cout<<iteration<<"\t"<<sqrt(LoadImbalance.SumSquare())/sqrt(DF.SumSquare())<<"\n";
				flush(cout);
			}while ( sqrt(LoadImbalance.SumSquare())>ErrorTolerance*sqrt(DF.SumSquare()) && iteration<=this->MaxIterations);
			//End of equilibrium iterations//
			if(iteration>this->IterationsDesired) //Δεν εγινε σωστή σύγκλιση
			{
				NewStep=false;
				psi=0.5*psi;
				DU=psi*ro*DUStart;
				U=Temp+DU;
				UpdateNodes(U);
			}
			else
			{
				NewStep=true;
				Lmax=Lmax*psi;
				psi=1;
				U=Temp+DU;
				UpdateNodes(U);
				LPrevious+=ksi*DF;//Νέο φορτίο ισορροπίας
				totalLoadpercent+=(ksi-1)*icrpercent;
				//Γράφω τα στοιχεία στο αρχείο//
				fstream loads;
				loads.open(this->LoadVectorFileOut,ios::app);//Εξωτερικά φορτία
				printvector(loads,LPrevious,CurrentLoadStep);
				flush(loads);
				loads.close();
				fstream displ;
				displ.open(this->DisplacementsFileOut,ios::app);//Μετακινήσεις
				printvector(displ,U,CurrentLoadStep);
				displ.close();
				this->ReactionToFile(this->ReactionFileOut);//Αντιδράσεις
				this->InternalForceToFile(this->ForceFileOut);//Εσωτερικές Δυνάμεις
				this->StressAndStainToFile(this->ElementFileOut,true);//Κύριες Τάσεις και παραμορφώσεις
				this->StressAndStainToFile("Stresses.txt",false);//Τάσεις και παραμορφώσεις
			}//endif
		}while (NewStep==false); //end do while
		LNext=LPrevious;
	}//end while
}

void FEM::GeometryToFile(char *filename)
{
	for(int i=1;i<=this->numberofNodes;i++)
	{
		nodes[i-1]->GeometryToFile(filename);
	}//Output Nodes
	for(int iot=1;iot<=this->numberofelements;iot++)
	{
		elements[iot-1]->GeometryToFile(filename);
	}//output Elements
}

ColumnVector FEM::ReactionVector(bool update,int method)
{
	//Την Ρουτίνα αυτη την χρησιμοποιώ για να φτίαξω το διάνυσμα των αντιδράσεων 
	//που αναπτύσονται στην κατασκευή 
	int Vectorlength;
	Vectorlength=(GetVariables()*this->numberofNodes)-(this->LostDofsVector(this->numberofNodes+1));
	ColumnVector loads(Vectorlength);
	loads=0;
	for(int i=1;i<=this->numberofelements;i++)
	{
		int numberofnodes=elements[i-1]->dofs()*elements[i-1]->GetNumberOfNodes();
		ColumnVector elementforce=elements[i-1]->InternalForce(update,method);
		for(int j=1;j<=numberofnodes;j++)
		{
			loads( elements[i-1]->ConnectivityVector()(j) )+=elementforce(j);
		}//end for j
	}//end for i
	// Τωρα αφαιρώ απο το ανω διανυσμα τις εσωτερικές δυνάμεις
	int tempindex=0;//βοηθητικός μετρητης
	for(int k=1;k<=this->numberofNodes;k++)
	{
		for(int l=1;l<=GetVariables();l++)
		{
			if(nodes[k-1]->GetExist(l)==1)
			{
				tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
				if(nodes[k-1]->GetFixed(l)==0) //Ο βαθμός ελευθερίας δεν ειναι δεσμευμένος
				{
					loads(tempindex)=0;
				}//endif
			}//endif
		}//endfor l
	} //endfor k
	return loads;

}


void FEM::ReactionToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	printvector(kos,ReactionVector(false,Method),CurrentLoadStep);
	kos.close();
}


void FEM::InternalForceToFile(char *filename)
{
	fstream elem;
	elem.open(filename,ios::app);
	for(int i=1;i<=this->numberofelements;i++)
	{
		printvector(elem,i,elements[i-1]->InternalForce(false,Method),CurrentLoadStep);
	}
	elem.close();
}

void FEM::InitialiseGaussLibrary()
{

	this->GausspointsLibrary.ReSize(15);
	this->GausspointsLibrary(1) = 0.;				// GaussPoint Order One
	this->GausspointsLibrary(2) = 2.;				// GaussPoint Order One Weight
	this->GausspointsLibrary(3) = -1/sqrt(3);		// GaussPoint Order Two
	this->GausspointsLibrary(4) = 1/sqrt(3);		// GaussPoint Order Two
	this->GausspointsLibrary(5) = 1.;				// GaussPoint Order Two Weight
	this->GausspointsLibrary(6) = 1.;				// GaussPoint Order Two Weight
	this->GausspointsLibrary(7) = -sqrt(0.6);   	// GaussPoint Order Three
	this->GausspointsLibrary(8) = 0.0;				// GaussPoint Order Three
	this->GausspointsLibrary(9) = sqrt(0.6);	    // GaussPoint Order Three
	this->GausspointsLibrary(10) = 5/9.;	        // GaussPoint Order Three Weight
	this->GausspointsLibrary(11) = 8/9.;	            // GaussPoint Order Three Weight
	this->GausspointsLibrary(12) = 5/9.;			// GaussPoint Order Three Weight

}

ColumnVector FEM::GaussPoint(int OrderOfIntegration, int ElementDimension, int GaussPointNumber)
{

  ColumnVector GP(ElementDimension + 1); //Επιστρέφει τις συντεταγμένες του σημείου + το βάρος
  int startposition=0;
  startposition = (OrderOfIntegration-1)*OrderOfIntegration; // (Sum i)*2 for i< Orderofintegration = n*(n+1)/2 *2
	// where n = OderOfintegration -1
  if (ElementDimension == 1) 
  {
	  GP(1) = GausspointsLibrary(startposition + GaussPointNumber);
  }
  else if(ElementDimension == 2)
  {
	  GP(1) = GausspointsLibrary(startposition + ( (GaussPointNumber-1)%OrderOfIntegration) +1 );
	  GP(2) = GausspointsLibrary(startposition + (floor( (GaussPointNumber-1)/OrderOfIntegration ) +1));
	  GP(3) = GausspointsLibrary(startposition + OrderOfIntegration + ( (GaussPointNumber-1)%OrderOfIntegration) +1 ) *
		  GausspointsLibrary(startposition + OrderOfIntegration + (floor((GaussPointNumber-1)/OrderOfIntegration) +1 ));
  }
  return GP;

}

int FEM::countfixes()
{
	// this routine just counts how many fixes are in the structure
	int tempindex=0;
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			if( (nodes[i-1]->GetExist(j)==1) && (nodes[i-1]->GetFixed(j)==1) )
			{
				tempindex++;
			}//endif
		}//end for(j)
	}//end for(i)
	return tempindex;
}

int FEM::countkinimatic(int stage)
{
	int tempsum=0;// total number of kinimatic conditions
	for(int i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli kinimatikou periorismou....
	{
		FixCommand* TempCommand= this->Fixes[i-1];
		if(stage>=TempCommand->startingstage && stage<=TempCommand->endingstage)//Process Fix Command
		{
			tempsum+=TempCommand->shiftingtimes ;
		}
	}//Next Command
	return tempsum;
}

int FEM::GetGlobalDof(int node,int dof)
{
		//οπου i-1 οι κόμβοι που προηγήθηκαν του εξεταζόμενου κόμβου
		//(DIM-1)*3 οι βαθμοί ελευθερίας καθε κόμβου ασχετα αν μερικοί είναι ανενεργοί
		//δηλαδή 3 για DIM=2 και 6 για DIM=3
		// και fem->LostDosfs(i) επιστρέφει τους βαθμούς ελευθερίας που ΔΕΝ
		// υπάρχουν απο τον 1 μεχρι και τον i-1 κόμβο
		int initialdofnumber=(node-1)*this->GetVariables()-this->LostDofsVector(node);
		if(dof>this->GetVariables()){return -1;}
		for(int k=1;k<=dof;k++)
		{
			if( nodes[node-1]->GetExist(k)==1 )
			{
				//Αν ο βαθμός ελευθερίας είναι ενεργός στο σύστημα
				//τότε απλά αυξάνω τον initialdofnumber
				initialdofnumber+=1;
			}
		}//end for(k)
		return initialdofnumber;
}


Matrix FEM::CreateLagrangeMatrix(int stage)
{
	int m,n,l;
	m=countfixes();n=countkinimatic(stage);
	l=GetNumberOfActiveSystemDofs();
	Matrix B(m+n,l); //this matrix has rows as many as kinematic conditions
	B=0;
	//First i process fixes.
	int tempindex=0;
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			if( (nodes[i-1]->GetExist(j)==1) && (nodes[i-1]->GetFixed(j)==1) )
			{
				tempindex++;
				B(tempindex,GetGlobalDof(i,j))=1;
			}//endif
		}//end for(j)
	}//end for(i)

	//Now proccessing the kinimatic conditions
	for(i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= this->Fixes[i-1];
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
					int dofj = GetGlobalDof(TempCommand->dofs[2*(j-1)]+svalue*(l-1),TempCommand->dofs[2*(j-1)+1]);
					B(tempindex,dofj)= TempCommand->coeff[j-1];
				}//next coeff
			}//next  shift
		}// end if
	}//Next Command

	return B;
}


ColumnVector FEM::CreateLagrangeVector(int stage)
{
	int m,n;
	m=countfixes();n=countkinimatic(stage);
	ColumnVector B(m+n); //this vector rows as many as kinematic conditions
	B=0;
	//First i process fixes.
	int tempindex=0;
	for(int i=1;i<=this->numberofNodes;i++)
	{
		for(int j=1;j<=GetVariables();j++)
		{
			if( (nodes[i-1]->GetExist(j)==1) && (nodes[i-1]->GetFixed(j)==1) )
			{
				tempindex++;
				B(tempindex)=0;// dof is fixed
			}//endif
		}//end for(j)
	}//end for(i)

	//Now proccessing the kinimatic conditions
	for(i=1;i<=GetFixCommandsNumber();i++)//gia kathe entoli....
	{
		FixCommand* TempCommand= this->Fixes[i-1];
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


FEM::LagrangeStiffnessMatrix(StiffnessMatrix &pinakas)
{

	double penalty=GetPenalty();
	Matrix B=CreateLagrangeMatrix(CurrentStage);
	int m=B.Ncols(); int n=B.Nrows();
	Matrix Zero(n,n); Zero=0;
	pinakas = pinakas + penalty*B.t()*B;
	pinakas = (pinakas | B.t()) & (B | Zero) ; //& concatenates vertically and | horizontally
}


FEM::LagrangeLoadVector(ColumnVector &Load)
{
	double penalty=GetPenalty();
	ColumnVector V=CreateLagrangeVector(CurrentStage);
	Matrix B=CreateLagrangeMatrix(CurrentStage);
	Load = Load + penalty*B.t()*V;
	Load = Load & V;
}


void FEM::CalibrateVector(ColumnVector& vector,bool fixed,bool lagrange)
{
	// if fixed = true then it returns the initial vector with zeros in all  fixed dofs
	// if fixed = false then it returns the initial vector with zeros in all non fixed dofs
//*	
	int tempindex=0;
	for(int k=1;k<=this->numberofNodes;k++)
	{
		for(int l=1;l<=GetVariables();l++)
		{
			if(nodes[k-1]->GetExist(l)==1)
			{
				tempindex++;//ο βαθμός ελευθερίας υπάρχει αρα πιάνει μια θέση στο διάνυσμα
				if(nodes[k-1]->GetFixed(l)==fixed) //Ο βαθμός ελευθερίας  ειναι του ζητούμενου τύπου
				{
					vector(tempindex)=0;
				}//endif
			}//endif
		}//endfor l
	} //endfor k 
//*/
/*
	Matrix B;B=ConstrainCenter->CreateConstrainsMatrix(CurrentStage);
	int m=B.Ncols(); int n=B.Nrows();
	for(int i=1;i<=m;i++)
	{
		for(int j=1;j<=n;j++)
		{
			if(B(j,i)!=0)
			{
				vector(i)=0;
			}
		}
	}
	
	if(lagrange)// zero lagrange multipliers
	{
		vector.Rows(m+1,vector.Nrows())=0;
	}
	*/
}








/*
bool FEM::CreateMessForGid(char* filename)
{
	int NofNodes = numberofNodes;
	int NofElem = numberofelements ;
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"MESH dimension "<<DIM<<" ElemType "<<"Linear"<<"  Nnode 8";
	kos<<"\nCoordinates";
	kos<<"\n# node number   coordinate_x  coordinate_y";
	for(int i=1;i<=NofNodes;i++)
	{
		kos<<"\n";
		kos.width(5);kos<<i;
		for(int j=1;j<=DIM;j++)
		{
			kos.width(15);
			kos<<nodes[i-1]->GetCoord(j);
		}
	}//end i
	kos.width(0);
	kos<<"\nend coordinates";
	kos<<"\nElements";
	kos<<"\n# element  node_1   node_2  material_number";
	for(i=1;i<=NofElem;i++)
	{
//		int NumberOfElementNodes=this->GetElement(i)->GetNumberOfNodes();
		kos<<"\n";
		kos.width(5);kos<<i;
		elements[i-1]->MeshToGid(kos);
//		kos.width(5);
//		kos<<1;
	}
	kos.width(0);
	kos<<"\nend elements\n";
	kos.close();
	return true;
}*/


bool FEM::CreateMessForGid(char* filename)
{
	int NofNodes = this->GetNumberofNodes();
	int NofElem = numberofelements;
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"MESH \"Quads8\" dimension "<<this->GetGeometry()<<" ElemType "<<"Quadrilateral"<<"  Nnode 8";
	kos<<"\nCoordinates";
	kos<<"\n# node number   coordinate_x  coordinate_y";
	for(int i=1;i<=NofNodes;i++)
	{
		kos<<"\n";
		kos.width(5);kos<<i;
		for(int j=1;j<=this->GetGeometry();j++)
		{
			kos.width(15);
			kos<<this->GetNode(i)->GetCoord(j);
		}
	}//end i
	kos.width(0);
	kos<<"\nend coordinates";
	kos<<"\nElements";
	kos<<"\n# element  node_1   node_2  material_number";
	for(i=1;i<=NofElem;i++)
	{
		int NumberOfElementNodes=this->GetElement(i)->GetNumberOfNodes();
		if(NumberOfElementNodes==8)
		{
			kos<<"\n";
			kos.width(5);kos<<i;
			this->GetElement(i)->MeshToGid(kos);
//			kos.width(5);
//			kos<<1;
		}
	}
	kos.width(0);
	kos<<"\nend elements\n";
	kos.close();
	return true;
}
