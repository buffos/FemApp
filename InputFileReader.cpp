#include "include.h"
#include "includeS.h"
#include "FEM.h"
#define _MATRICES_
#define _IOHEADERS_
#define _STRINGS_
#define _DEFINITIONS_

//#include "time.h"
//#include "Newmatio.h"
//#include <fstream.h>
//#include <conio.h>
//#include <ctype.h>
//#include "Str.h"

#include "FemIncludes.h"


#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

//void waitt();
void ProcessNodeFile(const char* filename,FEM& fem);
void ProcessControlFile(const char* filename,FEM& fem);
void ProcessElementFile(const char* filename,FEM& fem);
void ProcessPropertyFile(const char* filename,FEM& fem);
void ProcessMaterialFile(const char* filename,FEM& fem);
void copystring(char* s1,char* s2);
int stringsize(char* s1);

void main()
{
	FEM appl;
	char ControlFile[30];
	cout<<"\n Give the ControlFileName  ";
	cin>>ControlFile;
	String temp;
	cout<<"\n Processing Control File  "<<ControlFile;
	ProcessControlFile(ControlFile,appl);
	temp="Node";temp+=ControlFile;
	cout<<"\n Processing Node File  "<<temp;
	ProcessNodeFile(temp.data(),appl);
	temp="Elem";temp+=ControlFile;
	cout<<"\n Processing Element File  "<<temp;
	ProcessElementFile(temp.data(),appl);
	temp="Prop";temp+=ControlFile;
	cout<<"\n Processing Properies File  "<<temp;
	ProcessPropertyFile(temp.data(),appl);
	temp="Mat";temp+=ControlFile;
	cout<<"\n Processing Material File  "<<temp<<"\n";
	ProcessMaterialFile(temp.data(),appl);
	appl.Initialize();
	appl.GeometryToFile(appl.GeometryFileOut);
	appl.CreateMessForGid("problem.msh");
	if(appl.SolvingMethod==0)
	{
		appl.ElasticSolve(appl.SolverType);
	}
	else
	{
		appl.NonLinearSolve();
	}	
}

void ProcessControlFile(const char* filename,FEM& fem)
{
	String token;
	String EchoFilename="Echo";//Θα κανει echo τις επιλογές που διαβάζουμε για να μπορώ να κάνω έλεγχο
	double numdouble;
	EchoFilename+=filename;
	EchoFilename+='\0';
	unsigned long num;//Μεταβλητή που θα κρατάει τις τιμές που θα διαβάζω απο το αρχείο
	char temptoken[30];//Ομοια με πάνω
	ofstream echo(EchoFilename.data());//To echo αρχειο
	ifstream kos;
	kos.open(filename);
	while (!kos.eof() )
	{
		kos>>token;
		if (token=="NELEM")
		{
			kos>>num;
			fem.numberofelements=num;
//			fem.elements=new Element*[fem.numberofelements];
			echo<<"\nΟ αριθμός των στοιχείων τέθηκε ίσος με "<<num;
			flush(echo);
		}
		else if (token=="NNODE")
		{
			kos>>num;
			fem.SetNumberofNodes(num);
//			fem.nodes=new Node*[fem.GetNumberofNodes()];
			echo<<"\nΟ αριθμός των κόμβων τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="NMATL")
		{
			kos>>num;
			fem.numberofmaterials=num;
//			fem.Ylika=new Material*[fem.numberofmaterials];
			echo<<"\nΟ αριθμός των υλικών τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="NPROP")
		{
			kos>>num;
			fem.numberofproperties=num;
//			fem.Diatomes=new Property*[fem.numberofproperties];
			echo<<"\nΟ αριθμός των Ιδιοτήτων τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="NSTAG")
		{
			kos>>num;
			fem.NumberOfStages=num;
			echo<<"\nΟ αριθμός των Σταδίων Φόρτισης τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="NFIXS")
		{
			kos>>num;
			fem.SetFixCommandsNumber(num);
//			fem.Fixes=new FixCommand*[num];
			echo<<"\nΟ αριθμός των Περιορισμων μετακινήσεων τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="NSTCO") // Number of Stage Commands
		{
			kos>>num;
			fem.NumberOfStageCommands=num;
//			fem.Stages=new StageCommand*[fem.NumberOfStageCommands];
			echo<<"\nΟ αριθμός των εντολών φορτίου τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="PENSC") // PenaltyScalefactor
		{
			kos>>num;
			fem.SetPenaltyScaleFactor(num);
			echo<<"\nΟ αριθμός PenaltyScalefactor τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="DMVAR")//Αριθμός βαθμών ελευθερίας
		{
			kos>>num;
			fem.SetVariables(num);
			echo<<"\nΟ αριθμός των βαθμών ελευθερίας του προβλήματος τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="GEOME") // Διάσταση χώρου γεωμετρίας
		{
			kos>>num;
			fem.SetGeometry(num);
			echo<<"\nΗ διάσταση του χώρου τέθηκε ίση με  "<<num;
			flush(echo);
		}
		else if (token=="NITEG")
		{
			kos>>num;
			fem.NumberofItegrationLoops=num;
			echo<<"\nΟ αριθμός των βημάτων ολοκλήρωσης τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="MAXIT")
		{
			kos>>num;
			fem.MaxIterations=num;
			echo<<"\nΟ αριθμός των βημάτων φόρτισης σε κάθε βήμα ισσοροπίας τέθηκε ίσος με  "<<num;
			flush(echo);
		}
		else if (token=="ITERD")
		{
			kos>>num;
			fem.IterationsDesired=num;
			echo<<"\nΗ μεταβλητή IterationsDesired τέθηκε ίση με  "<<num;
			flush(echo);
		}
		else if (token=="METHD")
		{
			kos>>num;
			fem.Method=num;
			echo<<"\nΕπιλέχθηκε η μέθοδος "<<num;
			flush(echo);
		}
		else if (token=="FELEM")
		{
			kos>>temptoken;
			fem.ElementFileOut=new char[stringsize(temptoken)];
			copystring(fem.ElementFileOut,temptoken);
			echo<<"\nΤο αρχείο εξόδου αποτελεσμάτων είναι το "<<fem.ElementFileOut;
			flush(echo);
		}
		else if (token=="FGEOM")
		{
			kos>>temptoken;
			fem.GeometryFileOut=new char[stringsize(temptoken)];
			echo<<"\nΤο αρχείο εξόδου Γεωμετρίας είναι το "<<fem.GeometryFileOut;
			flush(echo);
		}
		else if (token=="FDISP")
		{
			kos>>temptoken;
			fem.DisplacementsFileOut=new char[stringsize(temptoken)];
			echo<<"\nΤο αρχείο εξόδου Μετακινήσεων είναι το "<<fem.DisplacementsFileOut;
			flush(echo);
		}
		else if (token=="FREAC")
		{
			kos>>temptoken;
			fem.ReactionFileOut=new char[stringsize(temptoken)];
			echo<<"\nΤο αρχείο εξόδου Αντιδράσεων Στηρίξεων είναι το "<<fem.ReactionFileOut;
			flush(echo);
		}		
		else if (token=="FFORC")
		{
			kos>>temptoken;
			fem.ForceFileOut=new char[stringsize(temptoken)];
			echo<<"\nΤο αρχείο εξόδου Εσωτερικών Δυνάμεων είναι το "<<fem.ForceFileOut;
			flush(echo);
		}		
		else if (token=="SOLVT")
		{
			kos>>num;
			fem.SolverType=num;
			echo<<"\nΗ μέθοδος επίλυσης του συστήματος εξισώσεων είναι η  "<<fem.SolverType;
			flush(echo);
		}		
		else if (token=="SMETH")
		{
			kos>>num;
			fem.SolvingMethod=num;
			echo<<"\nΗ μέθοδος επίλυσης του προβλήματος είναι η  "<<fem.SolvingMethod;
			flush(echo);
		}		
		else if (token=="LINCR")
		{
			kos>>numdouble;
			fem.LoadIncreament=numdouble;
			echo<<"\nΗ αύξηση του φορτίου σε κάθε βήμα επίλυσης είναι  "<<fem.LoadIncreament;
			flush(echo);
		}		
		else if (token=="ERRTL")
		{
			kos>>numdouble;
			fem.ErrorTolerance=numdouble;
			echo<<"\nΤο επιτρεπόμενο σφάλμα είναι  "<<fem.ErrorTolerance;
			flush(echo);
		}		
		else if (token=="ERRDTL")
		{
			kos>>numdouble;
			fem.ErrorDispTolerance=numdouble;
			echo<<"\nΤο επιτρεπόμενο σφάλμα για τις μετακινήσεις είναι  "<<fem.ErrorDispTolerance;
			flush(echo);
		}		
		else if (token=="WINFO")
		{
			kos>>num;
			fem.WriteExtraInfo=num;
			echo<<"\nΗ τιμή της Write Info To file είναι  "<<fem.WriteExtraInfo ;
			flush(echo);
		}		
		else if (token=="AUTOI")
		{
			kos>>num;
			fem.AutoIncrement=num;
			echo<<"\nΗ τιμή της AutoIncrement είναι  "<<fem.AutoIncrement;
			flush(echo);
		}		
		else if (token=="ENARC")
		{
			kos>>num;
			fem.EnableArcLengthMethod=num;
			echo<<"\nΗ τιμή της EnableArcLengthMethod είναι  "<<fem.EnableArcLengthMethod;
			flush(echo);
		}		
		else if (token=="SWARC")
		{
			kos>>num;
			fem.SwitchToArcLength=num;
			echo<<"\nΗ τιμή της SwitchToArcLength είναι  "<<fem.SwitchToArcLength;
			flush(echo);
		}		
		else if (token=="ACCEL")
		{
			kos>>num;
			fem.Accelerate=num;
			echo<<"\nΗ τιμή της Accelerate είναι  "<<fem.Accelerate;
			flush(echo);
		}		
		else if (token=="RSTRT")
		{
			kos>>num;
			fem.Restart=num;
			echo<<"\nΗ τιμή της Restart είναι  "<<fem.Restart;
			flush(echo);
		}		
		else if (token=="MXRSD")
		{
			kos>>num;
			fem.MaxResidualCalculations =num;
			echo<<"\nΗ τιμή της MaxResidualCalculations είναι  "<<fem.MaxResidualCalculations;
			flush(echo);
		}		
		else if (token=="MXLNS")
		{
			kos>>num;
			fem.MaxNoLineSearches =num;
			echo<<"\nΗ τιμή της MaxNoLineSearches είναι  "<<fem.MaxNoLineSearches;
			flush(echo);
		}		
		else if (token=="TOLRT")
		{
			kos>>numdouble;
			fem.ToleranceOnRatio =numdouble;
			echo<<"\nΗ τιμή της ToleranceOnRatio είναι  "<<fem.ToleranceOnRatio;
			flush(echo);
		}		
		else if (token=="MXAMP")
		{
			kos>>numdouble;
			fem.MaxAmpAtAnyStemp =numdouble;
			echo<<"\nΗ τιμή της MaxAmpAtAnyStemp είναι  "<<fem.MaxAmpAtAnyStemp;
			flush(echo);
		}		
		else if (token=="MXTSL")
		{
			kos>>numdouble;
			fem.MaxTotalStepLength =numdouble;
			echo<<"\nΗ τιμή της MaxTotalStepLength είναι  "<<fem.MaxTotalStepLength;
			flush(echo);
		}		
		else if (token=="MNTSL")
		{
			kos>>numdouble;
			fem.MinTotalStepLength =numdouble;
			echo<<"\nΗ τιμή της MinTotalStepLength είναι  "<<fem.MinTotalStepLength;
			flush(echo);
		}		
		else if (token=="DNITR")
		{
			kos>>num;
			fem.DesiredNumberOfIterations =num;
			echo<<"\nΗ τιμή της DesiredNumberOfIterations είναι  "<<fem.DesiredNumberOfIterations;
			flush(echo);
		}		
		else if (token=="MXLDI")
		{
			kos>>numdouble;
			fem.MaxLoadIncrement =numdouble;
			echo<<"\nΗ τιμή της MaxLoadIncrement είναι  "<<fem.MaxLoadIncrement;
			flush(echo);
		}		
		else if (token=="MNLDI")
		{
			kos>>numdouble;
			fem.MinLoadIncrement =numdouble;
			echo<<"\nΗ τιμή της MinLoadIncrement είναι  "<<fem.MinLoadIncrement;
			flush(echo);
		}		
		else if (token=="LSTIF")
		{
			kos>>numdouble;
			fem.LowestStiffnessSwitch =numdouble;
			echo<<"\nΗ τιμή της LowestStiffnessSwitch είναι  "<<fem.LowestStiffnessSwitch;
			flush(echo);
		}		
		else if (token=="DLINC")
		{
			kos>>numdouble;
			fem.DesiredLengthIncrement =numdouble;
			echo<<"\nΗ τιμή της DesiredLengthIncrement είναι  "<<fem.DesiredLengthIncrement;
			flush(echo);
		}		
		else if (token=="MXLIN")
		{
			kos>>numdouble;
			fem.MaxLengthIncrement =numdouble;
			echo<<"\nΗ τιμή της MaxLengthIncrement είναι  "<<fem.MaxLengthIncrement;
			flush(echo);
		}		
		else if (token=="MNLIN")
		{
			kos>>numdouble;
			fem.MinLengthIncrement =numdouble;
			echo<<"\nΗ τιμή της MinLengthIncrement είναι  "<<fem.MinLengthIncrement;
			flush(echo);
		}		
		else if (token=="ICRLF")
		{
			kos>>numdouble;
			fem.IncrLoadFactor =numdouble;
			echo<<"\nΗ τιμή της IncrLoadFactor είναι  "<<fem.IncrLoadFactor;
			flush(echo);
		}		
		else if (token=="ARMTH")
		{
			kos>>num;
			fem.ArcLengthAlgorithm =num;
			echo<<"\nΕπιλέχθηκε η "<<fem.ArcLengthAlgorithm<<"η  Μέθοδος ";
			flush(echo);
		}		
		else if (token=="FULNR")
		{
			kos>>num;
			fem.FullNR =num;
			echo<<"\nΗ τιμή της IncrLoadFactor είναι  "<<fem.FullNR;
			flush(echo);
		}		
		else if (token=="NORTY")
		{
			kos>>num;
			fem.NormType =num;
			echo<<"\nΗ τιμή της NormType είναι  "<<fem.NormType;
			flush(echo);
		}		
	}//end while
	kos.close();
	echo.close();
}//end of ProcessControlFile




void ProcessNodeFile(const char* filename,FEM& fem)
{
	int geometry = fem.GetGeometry();
	int variables = fem.GetVariables();

	String token;
	unsigned long index;//Μεταβλητή που θα κρατάει τον index του κομβου
	unsigned long StartingNode,EndingNode,FromStage,ToStage,keep;
	double coord[3],Dcoord[3];//Κρατά τις συντεταγμένς του κόμβου
	int exists,step,dof,i;
	double fix,Load,coeff,dofs;
	int xnodes,ynodes;//Number of grid nodes in x and y direction
	double *xgrid,*ygrid;
	int indexoffixcommand=0;
	Node* TempNode;
	StageCommand* TempStage;
	FixCommand* TempFix;
	ifstream kos;
	kos.open(filename);
	while (!kos.eof() )
	{
		kos>>token;
		if (token=="C")
		{
			kos>>index;
			for(i=0;i<=geometry-1;i++)
			{
				kos>>coord[i];//Διαβάζω τις συντεταγμένες
			}
			if (geometry==2)
			{
				TempNode=new Node(&fem,index,coord[0],coord[1]);
			}
			else
			{
				TempNode=new Node(&fem,index,coord[0],coord[1],coord[2]);
			}
			fem.SetNode(TempNode,index);//fem.nodes[index-1]=TempNode;//Θέτω τον κόμβο που έφτιαξα στον πίνακα
		}
		else if (token=="S")
		{

			kos>>StartingNode;
			kos>>EndingNode;
			kos>>step;
			for(i=0;i<=geometry-1;i++)
			{
				kos>>coord[i];//Διαβάζω τις συντεταγμένες
			}
			for(i=0;i<=geometry-1;i++)
			{
				kos>>Dcoord[i];//Διαβάζω τις διαφορές τωνσυντεταγμένων μεταξύ των κόμβων
			}
			int counter=0;
			for(i=StartingNode;i<=EndingNode;i+=step)
			{
				if (geometry==2)
				{
					TempNode=new Node(&fem,i,coord[0]+counter*Dcoord[0],coord[1]+counter*Dcoord[1]);
				}
				else
				{
					TempNode=new Node(&fem,i,coord[0]+counter*Dcoord[0],coord[1]+counter*Dcoord[1],coord[2]+counter*Dcoord[2]);
				}
				fem.SetNode(TempNode,i);
				counter++;
			}//end of Node creation
		}//end of Serial Node creation
		else if (token=="W")
		{

			kos>>StartingNode;
			kos>>EndingNode;
			kos>>step;
			for(i=0;i<=geometry-1;i++)
			{
				kos>>coord[i];//Διαβάζω τις συντεταγμένες
			}
			for(i=0;i<=geometry-1;i++)
			{
				kos>>Dcoord[i];//Διαβάζω τις διαφορές τωνσυντεταγμένων μεταξύ των κόμβων
			}
			int counter=0;
			for(i=StartingNode+step;i<EndingNode;i+=step)
			{
				if (geometry==2)
				{
					TempNode=new Node(&fem,i,coord[0]+counter*Dcoord[0],coord[1]+counter*Dcoord[1]);
				}
				else
				{
					TempNode=new Node(&fem,i,coord[0]+counter*Dcoord[0],coord[1]+counter*Dcoord[1],coord[2]+counter*Dcoord[2]);
				}
				fem.SetNode(TempNode,i);
				counter++;
			}//end of Node creation
		}//end of Serial Node creation without ends
		else if (token=="G")
		{

			kos>>StartingNode;
			kos>>EndingNode;
			kos>>step;
			kos>>xnodes>>ynodes;
			xgrid = new double[xnodes];
			ygrid = new double[ynodes];
			for(i=0;i<=xnodes-1;i++)
			{
				kos>>xgrid[i];//Διαβάζω τις συντεταγμένες
			}
			for(i=0;i<=ynodes-1;i++)
			{
				kos>>ygrid[i];//Διαβάζω τις συντεταγμένες
			}
			int counter=StartingNode;
			if (xnodes<=ynodes)
			{
				for(i=0;i<=ynodes-1;i++)
				{
					for(int k=0;k<=xnodes-1;k++)
					{
						if (geometry==2)
						{
							TempNode=new Node(&fem,counter,xgrid[k],ygrid[i]);
							fem.SetNode(TempNode,counter);
							counter=counter+step;
						}
						else
						{
							//To write code
						}
					}//end for k
				}//end for i
			}
			else
			{
				for(i=0;i<=xnodes-1;i++)
				{
					for(int k=0;k<=ynodes-1;k++)
					{
						if (geometry==2)
						{
							TempNode=new Node(&fem,counter,xgrid[i],ygrid[k]);
							fem.SetNode(TempNode,counter);
							counter=counter+step;
						}
						else
						{
							//To write code
						}
					}//end for k
				}//end for i
			}//end of Node creation
		}//end of Grid Node creation
		else if (token=="L")//Define Loads at Nodes
		{
			kos>>dof;
			kos>>Load;
			kos>>StartingNode;
			kos>>EndingNode;
			kos>>step;
			for(int i=StartingNode;i<=EndingNode;i+=step)
			{
				fem.GetNode(i)->SetLoad(dof,Load);
			}
		}//end if
		else if (token=="E")//Define Exist at Nodes
		{
			kos>>dof;
			kos>>exists;
			kos>>StartingNode;
			kos>>EndingNode;
			kos>>step;
			for(int i=StartingNode;i<=EndingNode;i+=step)
			{
				fem.GetNode(i)->SetExist(dof,exists);
			}
		}//end if
		else if (token=="FE")//Define Cinimatic Conditions at Nodes
		{
			TempFix=new FixCommand;
			kos>>TempFix->startingstage>>TempFix->endingstage;
			kos>>dof;// arithmos twn bathmwn eleutherias poy tha einai stin eksisosi...
			for(int i= 1;i<=dof;i++)
			{
				kos>>coeff;
				TempFix->coeff.push_back(coeff);//reading coefficients
			}
			for(i= 1;i<=2*dof;i++)
			{
				kos>>dofs;
				TempFix->dofs.push_back(dofs);//reading dofindices
			}
			kos>>fix;TempFix->result=fix;
			kos>>TempFix->shiftingtimes;
			kos>>TempFix->shiftvalue;
			indexoffixcommand++;
			fem.SetFixCommand(TempFix,indexoffixcommand);
		}
		else if (token=="FS")//Define Fix Conditions at Nodes
		{
			kos>>dof;
			kos>>fix;
			kos>>StartingNode;
			kos>>EndingNode;
			kos>>step;
			for(int i=StartingNode;i<=EndingNode;i+=step)
			{
				fem.GetNode(i)->SetFixed(dof,fix);
			}
		}//end if
		else if (token=="SC")
		{
			kos>>index;
			kos>>StartingNode>>EndingNode;
			kos>>step>>dof>>Load>>FromStage>>ToStage>>keep;
			TempStage=new StageCommand;
			TempStage->dof=dof;TempStage->fromnode=StartingNode;
			TempStage->fromstage=FromStage;TempStage->keep=keep;
			TempStage->load=Load;TempStage->step=step;
			TempStage->tonode=EndingNode;TempStage->tostage=ToStage;
			fem.SetStageCommand(TempStage,index);
		}//end if

	}//end while
}//end of Reading Node File



void ProcessElementFile(const char* filename,FEM& fem)
{
	int index,node1,node2,node3,node4,matnumber,propnumber,i,counter,dof;
	int node5,node6,node7,node8;
	int startingindex,endingindex,stepindex,dnode1,dnode2,dnode3,dnode4,dnode5,dnode6,dnode7,dnode8;
	int dynode1,dynode2,dynode3,dynode4,dynode5,dynode6,dynode7,dynode8;
	int xgridelements,ygridelements; //Number of elements in the x and y direction
	int integrationorder;
	Element* elem;
	String token;
	ifstream kos;
	kos.open(filename);

//	int geometry = fem.GetGeometry();
//	int variables = fem.GetVariables();


	while (!kos.eof() )
	{
		kos>>token;
		if (token=="CBEAM2D")
		{
			kos>>index;
			kos>>node1>>node2>>matnumber>>propnumber;
			elem=new Beam2D(&fem,index,node1,node2);
			elem->SetMaterial(matnumber);
			elem->SetProperty(propnumber);
			fem.SetElement(elem,index);
		}
		if (token=="CSPRG2D")
		{
			kos>>index;
			kos>>node1>>node2>>matnumber>>propnumber;
			elem=new Spring2D(&fem,index,node1,node2);
			elem->SetMaterial(matnumber);
			elem->SetProperty(propnumber);
			fem.SetElement(elem,index);
		}
		if (token=="CSSPR2D")
		{
			kos>>index>>dof;
			kos>>node1>>node2>>matnumber>>propnumber;
			elem=new SpringXYZ2D(&fem,index,dof,node1,node2);
			elem->SetMaterial(matnumber);
			elem->SetProperty(propnumber);
			fem.SetElement(elem,index);
		}
		if (token=="CTRSS2D")
		{
			kos>>index;
			kos>>node1>>node2>>matnumber>>propnumber;
			elem=new Truss2D(&fem,index,node1,node2);
			elem->SetMaterial(matnumber);
			elem->SetProperty(propnumber);
			fem.SetElement(elem,index);
		}
		if (token=="CQU4N2D")
		{
			kos>>index;
			kos>>node1>>node2>>node3>>node4>>matnumber>>propnumber>>integrationorder;
			elem=new FourNode2DElasticityElement(&fem,index,node1,node2,node3,node4,integrationorder);
			elem->SetMaterial(matnumber);
			elem->SetProperty(propnumber);
			fem.SetElement(elem,index);
		}
		if (token=="CQU8N2D")
		{
			kos>>index;
			kos>>node1>>node2>>node3>>node4>>node5>>node6>>node7>>node8>>matnumber>>propnumber>>integrationorder;
			elem=new EightNode2DElasticityElement(&fem,index,node1,node2,node3,node4,node5,node6,node7,node8,integrationorder);
			elem->SetMaterial(matnumber);
			elem->SetProperty(propnumber);
			fem.SetElement(elem,index);
		}
		if (token=="SBEAM2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>node1>>node2>>dnode1>>dnode2>>matnumber>>propnumber;
			counter=0;
			for(i=startingindex;i<=endingindex;i+=stepindex)
			{
				elem=new Beam2D(&fem,i,node1+counter*dnode1,node2+counter*dnode2);
				elem->SetMaterial(matnumber);
				elem->SetProperty(propnumber);
				fem.SetElement(elem,i);
				counter++;
			}//end for
		}//endif
		if (token=="SSPRG2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>node1>>node2>>dnode1>>dnode2>>matnumber>>propnumber;
			counter=0;
			for(i=startingindex;i<=endingindex;i+=stepindex)
			{
				elem=new Spring2D(&fem,i,node1+counter*dnode1,node2+counter*dnode2);
				elem->SetMaterial(matnumber);
				elem->SetProperty(propnumber);
				fem.SetElement(elem,i);
				counter++;
			}//end for
		}//endif
		if (token=="SSSPR2D")
		{
			kos>>startingindex>>endingindex>>stepindex>>dof;
			kos>>node1>>node2>>dnode1>>dnode2>>matnumber>>propnumber;
			counter=0;
			for(i=startingindex;i<=endingindex;i+=stepindex)
			{
				elem=new SpringXYZ2D(&fem,i,dof,node1+counter*dnode1,node2+counter*dnode2);
				elem->SetMaterial(matnumber);
				elem->SetProperty(propnumber);
				fem.SetElement(elem,i);
				counter++;
			}//end for
		}//endif
		if (token=="STRSS2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>node1>>node2>>dnode1>>dnode2>>matnumber>>propnumber;
			counter=0;
			for(i=startingindex;i<=endingindex;i+=stepindex)
			{
				elem=new Truss2D(&fem,i,node1+counter*dnode1,node2+counter*dnode2);
				elem->SetMaterial(matnumber);
				elem->SetProperty(propnumber);
				fem.SetElement(elem,i);
				counter++;
			}//end for
		}//endif
		if (token=="SQU4N2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>node1>>node2>>node3>>node4>>dnode1>>dnode2>>dnode3>>dnode4>>matnumber>>propnumber>>integrationorder;
			counter=0;
			for(i=startingindex;i<=endingindex;i+=stepindex)
			{
				elem=new FourNode2DElasticityElement(&fem,i,node1+counter*dnode1,node2+counter*dnode2,node3+counter*dnode3,node4+counter*dnode4,integrationorder);
				elem->SetMaterial(matnumber);
				elem->SetProperty(propnumber);
				fem.SetElement(elem,i);
				counter++;
			}//end for
		}//endif
		if (token=="SQU8N2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>node1>>node2>>node3>>node4>>node5>>node6>>node7>>node8;
			kos>>dnode1>>dnode2>>dnode3>>dnode4>>dnode5>>dnode6>>dnode7>>dnode8>>matnumber>>propnumber>>integrationorder;
			counter=0;
			for(i=startingindex;i<=endingindex;i+=stepindex)
			{
				elem=new EightNode2DElasticityElement(&fem,i,node1+counter*dnode1,node2+counter*dnode2,node3+counter*dnode3,node4+counter*dnode4,node5+counter*dnode5,node6+counter*dnode6,node7+counter*dnode7,node8+counter*dnode8,integrationorder);
				elem->SetMaterial(matnumber);
				elem->SetProperty(propnumber);
				fem.SetElement(elem,i);
				counter++;
			}//end for
		}//endif
		if (token=="GQU4N2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>xgridelements>>ygridelements>>node1>>node2>>node3>>node4;
			kos>>dnode1>>dnode2>>dnode3>>dnode4;
			kos>>dynode1>>dynode2>>dynode3>>dynode4>>matnumber>>propnumber>>integrationorder;
			counter=startingindex;
			for(i=0;i<=xgridelements-1;i++)
			{
				for(int j=0;j<=ygridelements-1;j++)
				{
					elem=new FourNode2DElasticityElement(&fem,counter,node1+i*dnode1+j*dynode1,node2+i*dnode2+j*dynode2,node3+i*dnode3+j*dynode3,node4+i*dnode4+j*dynode4,integrationorder);
					elem->SetMaterial(matnumber);
					elem->SetProperty(propnumber);
					fem.SetElement(elem,counter);
					counter=counter+stepindex;
				}
			}//end for
		}//endif
		if (token=="GQU8N2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>xgridelements>>ygridelements>>node1>>node2>>node3>>node4>>node5>>node6>>node7>>node8;
			kos>>dnode1>>dnode2>>dnode3>>dnode4>>dnode5>>dnode6>>dnode7>>dnode8;
			kos>>dynode1>>dynode2>>dynode3>>dynode4>>dynode5>>dynode6>>dynode7>>dynode8>>matnumber>>propnumber>>integrationorder;
			counter=startingindex;
			for(i=0;i<=xgridelements-1;i++)
			{
				for(int j=0;j<=ygridelements-1;j++)
				{
					elem=new EightNode2DElasticityElement(&fem,counter,node1+i*dnode1+j*dynode1,node2+i*dnode2+j*dynode2,node3+i*dnode3+j*dynode3,node4+i*dnode4+j*dynode4,node5+i*dnode5+j*dynode5,node6+i*dnode6+j*dynode6,node7+i*dnode7+j*dynode7,node8+i*dnode8+j*dynode8,integrationorder);
					elem->SetMaterial(matnumber);
					elem->SetProperty(propnumber);
					fem.SetElement(elem,counter);
					counter=counter+stepindex;
				}
			}//end for
		}//endif
		if (token=="GTRSS2D")
		{
			kos>>startingindex>>endingindex>>stepindex;
			kos>>xgridelements>>ygridelements>>node1>>node2;
			kos>>dnode1>>dnode2;
			kos>>dynode1>>dynode2>>matnumber>>propnumber;
			counter=startingindex;
			for(i=0;i<=xgridelements-1;i++)
			{
				for(int j=0;j<=ygridelements-1;j++)
				{
					elem=new Truss2D(&fem,counter,node1+i*dnode1+j*dynode1,node2+i*dnode2+j*dynode2);
					elem->SetMaterial(matnumber);
					elem->SetProperty(propnumber);
					fem.SetElement(elem,counter);
					counter=counter+stepindex;
				}
			}//end for
		}//endif
	}//end while
	kos.close();
}


void copystring(char* s1,char* s2)
{
	while(*s2)
	{
		*s1++=*s2++;
	}
	*s1='\0';
	
}

int stringsize(char* s1)
{
	int index=0;
	while(*s1)
	{
		index++;
		s1++;
	}
	return index;
}

void ProcessPropertyFile(const char* filename,FEM& fem)
{
	double CrossSection,Inertia,Thickness;
	int index=0;
	Property* prop;
	String token;
	ifstream kos;
	kos.open(filename);
	while (!kos.eof() )
	{
		kos>>token;
		if (token=="BEAM2D")
		{
			index++;
			kos>>CrossSection>>Inertia;
			prop=new BeamProperty(CrossSection,Inertia);
			fem.SetProperty(prop,index);
		}
		if (token=="PLAN2D")
		{
			index++;
			kos>>Thickness;
			prop=new PlaneElementProperty(Thickness);
			fem.SetProperty(prop,index);
		}
		if (token=="TRSS2D")
		{
			index++;
			kos>>CrossSection;
			prop=new TrussProperty(CrossSection);
			fem.SetProperty(prop,index);
		}
	}//end while
	kos.close();
}

void ProcessMaterialFile(const char* filename,FEM& fem)
{
	double num1,num2,num3,num4,num5,num6,num7;
	double fc,ec,ksi,ita;
	int confined,lockingcondition;
	int index=0;
	Material* mat;
	String token;
	ifstream kos;
	kos.open(filename);
	while (!kos.eof() )
	{
		kos>>token;
		if (token=="ELTRUS2D")
		{
			index++;
			kos>>num1;
			mat=new TrussMaterialS(num1);
			fem.SetMaterial(mat,index);
		}
		if (token=="NLTRUS2D")
		{
			index++;
			kos>>num1>>num2>>num3>>num4>>num5>>num6>>num7;
			mat=new TrussMaterialS(num1,num2,num3,num4,num5,num6,num7);
			fem.SetMaterial(mat,index);
		}
		if (token=="ELTRUT2D")
		{
			index++;
			kos>>num1;
			mat=new TrussMaterialS(num1);
			fem.SetMaterial(mat,index);
		}
		if (token=="NLTRUT2D")
		{
			index++;
			kos>>num1>>num2>>num3>>num4>>num5>>num6>>num7;
			mat=new TrussMaterial(num1,num2,num3,num4,num5,num6,num7);
			fem.SetMaterial(mat,index);
		}
		if (token=="PLANES2D")
		{
			index++;
			kos>>num1>>num2;
			mat=new PlaneStress(num1,num2);
			fem.SetMaterial(mat,index);
		}
		if (token=="VECIO2DT")
		{
			index++;
			kos>>fc>>ec>>confined;
			mat=new Veccio2DT(fc,ec,confined);
			fem.SetMaterial(mat,index);
		}
		if (token=="VECIO2DS")
		{
			index++;
			kos>>fc>>ec>>confined;
			mat=new Veccio2DS(fc,ec,confined);
			fem.SetMaterial(mat,index);
		}
		if (token=="VECI2DMS")
		{
			index++;
			kos>>fc>>ec>>confined>>lockingcondition>>ksi>>ita;
			mat=new Veccio2DModS(fc,ec,confined,lockingcondition,ksi,ita);
			fem.SetMaterial(mat,index);
		}
		if (token=="TRUSPL2D")
		{
			index++;
 			kos>>num1>>num2>>num3>>num4>>num5>>num6;
			mat=new TrussPolyonimial(num1,num2,num3,num4,num5,num6);
			fem.SetMaterial(mat,index);
		}
		if (token=="INHSSPR2D")
		{
			index++;
 			kos>>num1>>num2>>num3>>num4;
			mat=new InhomSerialSprings (num1,num2,num3,num4);
			fem.SetMaterial(mat,index);
		}
	}//end while
	kos.close();
}

