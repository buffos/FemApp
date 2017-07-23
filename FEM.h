// FEM.h: interface for the FEM class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEM_H__F4361C16_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_FEM_H__F4361C16_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define _CONSTRAINS_
#define _CONTOLPANEL_
#define _PARAMETERSPANEL_
#define _MATRICES_
#define _FEMIO_
#define _STAGES_
#define _EXTRAFUNCTIONS_
#define _DEFINITIONS_


///Elements
#include "Element.h"
#include "IsoparametricElement.h"
#include "FourNode2DElasticityElement.h"
#include "EightNode2DElasticityElement.h"
#include "Truss2D.h"
#include "Beam2D.h"
#include "Spring2D.h"
#include "SpringXYZ2D.h"
//Materials
#include "Material.h"
#include "TrussMaterial.h"
#include "TrussMaterialS.h"
#include "TrussPolyonimial.h"
#include "Veccio2DT.h"
#include "Veccio2DS.h"
#include "Veccio2DModS.h"
#include "PlaneStress.h"
#include "SpringMaterial.h"
#include "SpringMaterialS.h"
#include "InhomSerialSprings.h"
//Properties
#include "Property.h"
#include "TrussProperty.h"
#include "BeamProperty.h"
#include "PlaneElementProperty.h"
//Node
#include "Node.h"

#include "Element.h"

// Analysis
#include "ElasticSolver.h"
#include "ArcLength.h"


class FEM  
{
public:
	ControlCenter*  CCenter;
	ParametersCenter* PCenter;
	Constrains*     ConstrainCenter;
	//Gerenal Stuff//
	int NumberofItegrationLoops;
	int Method;//Method 1 Dε//Dσ Method 2 απο την σχεση Stress-Strain
	int SolverType;
	int SolvingMethod;//0=ElasticSolve 1=TangentStiffnessN_RComplete ....
	int NumberOfStageCommands; //Αριθμός εντολών σταδίου φόρτισης
	int CurrentLoadStep;
	bool ValidateStress;//Αν είναι true τοτε δεν επιτρέπω ενα σημείο να αλλάξει πρόσημο στην τάση
	double IncrLoadFactor;
	int CurrentStage;
	int	NumberOfStages;
	int IterationsDesired;//Ο επιθυμητος Ελάχιστος αριθμός επαναλήψεων στην Orthogonal Residual Method και στην Arc-Length
	//Πολύ χρήσιμο κατα την μη γραμμική αναλυση.
	//Αν εχω δυναμική φόρτιση τοτε πρέπει στη φάση του unloading να θέσω την τιμή του ίση με false
	char* ElementFileOut;
	char* GeometryFileOut;
	char* DisplacementsFileOut;
	char* ReactionFileOut;//Αρχείο με τις αντιδράσεις στηρίξεων
	char* ForceFileOut;// Με τις εσωτερικές δυνάεμις κάθε στοιχείου
	char* LoadVectorFileOut;
	char* CurrentStiffnessFile;
	char* InfoFile;
	char* GeneralDataFileName;
	double ErrorTolerance;
	double ErrorDispTolerance;
	double LoadIncreament;
	double InitialDeltaL;
	bool WriteExtraInfo;
	bool AutoIncrement;
	bool EnableArcLengthMethod;//If true then arc length is in use//
	bool Accelerate;
	int Restart;//restarting loadstep.. default =0;
	bool FullNR;//Use Full or modified Newton Raphson method
	int MaxIterations;//Ο μεγαλύτερος αριθμος επαναλήψεων που μπορει να γίνουν ωστε να μην θεωρηθεί μη συγκίνουσα η διαδικασία
	int MaxResidualCalculations; //Max number of residual calculations during line search
	int MaxNoLineSearches;
	bool DisplacementControl; // Use Load Control if false (Default) or Displ Control if true
	int NormType; // 1 if normal Euclidean Norm , 2 if Scaling is Used : Default = 1
	//End of General Parameters//
	//Line search parameters//
	double ToleranceOnRatio;
	double MaxAmpAtAnyStemp;
	double MaxTotalStepLength;
	double MinTotalStepLength;
	//End of Line Search Stuff//
	//AutoIncrementation Parameters//
	int DesiredNumberOfIterations;
	double MaxLoadIncrement;
	double MinLoadIncrement;
	bool SwitchToArcLength;
	double LowestStiffnessSwitch;
	//End of AutoIncrementation Parameters// 
	//Arc Length Method Parameters//
	double DesiredLengthIncrement;
	double ArcHint;
	double MaxLengthIncrement;
	double MinLengthIncrement;
	int ArcLengthAlgorithm;
	//End of Arc Length Parameters//
	//Gausspoint Library
	ColumnVector GausspointsLibrary;
	ColumnVector DisplacementAtPreviousStage;


	//Functions//
	ColumnVector GaussPoint(int OrderOfIntegration ,int ElementDimension,int GaussPointNumber);
	void WriteFileheaders();
	void InitialiseGaussLibrary();
	void CalculatePenalty();
	void ArcLengthMethod();
	void InternalForceToFile(char* filename);
	void ReactionToFile(char* filename);
	ColumnVector ReactionVector(bool update,int method);
	void GeometryToFile(char* filename);
	void OrthogonalResidualMethod();
	void StressAndStainToFile(char* filename,bool principal);
	void SecantStiffnessMethod();
	void SecantStiffnessStepMethod();
	void TangentStiffnessN_RIncomplete();
	void TangentStiffnessN_RComplete();
	ColumnVector AssembleInternalLoads(bool update,int method,bool withpenaltyloads=true);
	void NonLinearSolve();
	void CalibrateVector(ColumnVector& vector,bool fixed,bool lagrange);
	ColumnVector ReactionAtNode(int nodenumber);
	//Ενημερώνειτα Nodes με τα νέα displacements
	UpdateNodes(const ColumnVector& disp);
	//Συνάρτηση που εφαρμόζει την penalty method
	Penalty(StiffnessMatrix&  pinakas);
	//Use Lagrange multipliers to impose boundary conditions
	LagrangeStiffnessMatrix(StiffnessMatrix&  pinakas);
	LagrangeLoadVector(ColumnVector&  Load);
	//Επιλύει το  σύστημα και ενημερώνει τα Nodes με τα νέα displacements
	ElasticSolve(int solvertype);
	ColumnVector AssembleLoads();
	ColumnVector AssembleStageLoads(int fromstage,int tostage);
	ColumnVector AssembleStageLoads(int stage);
	ColumnVector AssemblePenaltyStageLoads(int stage);
	ColumnVector LostDofsVector;
	ColumnVector LostDofs();
	Matrix CreateLagrangeMatrix(int stage);
	ColumnVector CreateLagrangeVector(int stage);
	Initialize();
	int band;
	int systemband();
	int countfixes();
	int countkinimatic(int stage);
	int GetGlobalDof(int node,int dof);
	int GetNumberOfActiveSystemDofs(){return (GetVariables()*this->numberofNodes)-LostDofsVector(this->numberofNodes+1);}
	bool CreateMessForGid(char* filename);
	StiffnessMatrix AssembleStiffnessMatrix();

	FEM();
	virtual ~FEM();
	unsigned long numberofelements;
	unsigned long numberofmaterials;
	unsigned long numberofproperties;

	//Get Set Functions
	void SetGeometry(int num){geometry=num;}
	int  GetGeometry(){return geometry;}
	void SetVariables(int num){variables=num;}
	int  GetVariables(){return variables;}
	void SetPenalty(double num){PenaltyNumber=num;}
	double  GetPenalty(){return PenaltyNumber;}
	void SetFixCommandsNumber(int num){FixCommandsNumber=num;}
	int  GetFixCommandsNumber(){return FixCommandsNumber;}
	void IncrFixCommandsNumber(){FixCommandsNumber++;}
	void SetPenaltyScaleFactor(int num){PenaltyScaleFactor=num;}
	int  GetPenaltyScaleFactor(){return PenaltyScaleFactor;}

	unsigned long GetNumberofNodes(){return numberofNodes;}
	void SetNumberofNodes(unsigned long num){numberofNodes=num;}

	void SetNode(Node* n,int position){(nodes.size()<position)? nodes.resize(position,n) : nodes[position-1]=n;}
	void SetMaterial(Material* n,int position){(Ylika.size()<position)? Ylika.resize(position,n) : Ylika[position-1]=n;}
	void SetElement(Element* n,int position){(elements.size()<position)? elements.resize(position,n) : elements[position-1]=n;}
	void SetProperty(Property* n,int position){(Diatomes.size()<position)? Diatomes.resize(position,n) : Diatomes[position-1]=n;}
	void SetStageCommand(StageCommand* n,int position){(Stages.size()<position)? Stages.resize(position,n) : Stages[position-1]=n;}
	void SetFixCommand(FixCommand* n,int position){(Fixes.size()<position)? Fixes.resize(position,n) : Fixes[position-1]=n;}

	void SetNodeSize(int s){nodes.resize(s,0);}
	void SetMaterialSize(int s){Ylika.resize(s,0);}
	void SetElementSize(int s){elements.resize(s,0);}
	void SetPropertySize(int s){Diatomes.resize(s,0);}
	void SetStageCommandSize(int s){Stages.resize(s,0);}
	void SetFixCommandSize(int s){Fixes.resize(s,0);}

	Node* GetNode(int position){return nodes[position-1];}
	Material* GetMaterial(int position){return Ylika[position-1];}
	Element* GetElement(int position){return elements[position-1];}
	Property* GetProperty(int position){return Diatomes[position-1];}
	StageCommand* GetStageCommand(int position){return Stages[position-1];}
	FixCommand* GetFixCommand(int position){return Fixes[position-1];}

private:
	int geometry;
	int variables;
	int PenaltyScaleFactor;
	int FixCommandsNumber;
	unsigned long numberofNodes;
	double PenaltyNumber;


	std::vector<Node*> nodes;
	std::vector<Material*> Ylika;
	std::vector<Element*> elements;
	std::vector<Property*> Diatomes;
	std::vector<StageCommand*> Stages;
	std::vector<FixCommand*> Fixes;

};

#endif // !defined(AFX_FEM_H__F4361C16_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
