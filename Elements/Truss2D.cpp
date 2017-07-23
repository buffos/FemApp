// Truss2D.cpp: implementation of the Truss2D class.
//
//////////////////////////////////////////////////////////////////////

#include "Truss2D.h"
#include <math.h>
#include <iostream.h>
#include <fstream.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Truss2D::~Truss2D()
{

}

Truss2D::Truss2D(FEM* fem,int index, int firstnodenumber, int secondnodenumber)
{
	this->SetNumberOfGaussPoints(1);
	//Δυο κόμβους εχει το στοιχείο
	this->SetNumberofNodes(2);
	//Εξασφαλιζω θέση για αυτούς τους κόμβους
	this->SetNodesSize(2);
	this->fem=fem;
	//Τοποθετώ τον αριθμό του στοιχείου στο πρόβλημα
	this->SetIndex(index);
	this->SetNode(1,firstnodenumber);
	this->SetNode(2,secondnodenumber);
	this->Setdofcode(1,1);
	this->Setdofcode(2,1);
//	this->Setdofcode(3,0);
	this->Strain.ReSize(1);
	this->Stress.ReSize(1);
	this->Strain=0;
	this->Stress=0;
	fem->GetNode(firstnodenumber)->ElementAdd(this->GetIndex() );
	fem->GetNode(secondnodenumber)->ElementAdd(this->GetIndex() );

}

Matrix Truss2D::CreateStiffnessMatrix()
{
	double c;//cosine of truss angle
	double s;//sine of truss angle
	double l;//length of truss
	double A;//cross section of truss
	double E;//Μετρο Ελαστικότητας
	E=this->fem->GetMaterial(this->GetMaterial())->C(this,0)(1,1);
	A=this->fem->GetProperty(this->GetProperty())->GetParameter(1);
	double x1,y1,x2,y2;
	x1=fem->GetNode(this->GetNode(1))->GetCoord(1);
	y1=fem->GetNode(this->GetNode(1))->GetCoord(2);
	x2=fem->GetNode(this->GetNode(2))->GetCoord(1);
	y2=fem->GetNode(this->GetNode(2))->GetCoord(2);
	l=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) );
	c=(x2-x1)/l;
	s=(y2-y1)/l;
	SymmetricMatrix k(4);
	k(1,1)=c*c;
	k(1,2)=c*s;
	k(1,3)=-c*c;
	k(1,4)=-c*s;
	k(2,2)=s*s;
	k(2,3)=-c*s;
	k(2,4)=-s*s;
	k(3,3)=c*c;
	k(3,4)=c*s;
	k(4,4)=s*s;
	k=(E*A/l)*k;
	return k;

}

Truss2D::Truss2D()
{
	this->SetNodesSize(2);
	this->SetNumberofNodes(2);

}

ColumnVector Truss2D::GetStrain(int gausspointnumber)
{
	double c;//cosine of truss angle
	double s;//sine of truss angle
	double l;//length of truss
	double x1,y1,x2,y2;
	x1=fem->GetNode(this->GetNode(1))->GetCoord(1);
	y1=fem->GetNode(this->GetNode(1))->GetCoord(2);
	x2=fem->GetNode(this->GetNode(2))->GetCoord(1);
	y2=fem->GetNode(this->GetNode(2))->GetCoord(2);
	l=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) );
	c=(x2-x1)/l;
	s=(y2-y1)/l;
	ColumnVector displacements=GetDisplacementOfNodes();
	ColumnVector L(4);
	L(1)=-c;L(2)=-s;L(3)=c;L(4)=s;
	return (1./l)*L.t()*displacements;
}


ColumnVector Truss2D::GetStress(int gausspointnumber,bool update,int method)
{
	//Ενημερώνει την τάση στο Gausspoint κατα δεδομένο Fraction (χρησιμεύει στην ολοκήρωση)
	ColumnVector CurrentStress;
	ColumnVector PreviousStrain;
	ColumnVector PreviousStress;
	PreviousStrain=this->Strain;
	PreviousStress=this->Stress;
	ColumnVector CurrentStrain;
	CurrentStrain=GetStrain(gausspointnumber);
	if (method==1)
	{
		ColumnVector StrainDifference=(1.0/fem->NumberofItegrationLoops)*(CurrentStrain-PreviousStrain);
		Matrix C;
		for(int i=1;i<=fem->NumberofItegrationLoops;i++)
		{
			C=this->fem->GetMaterial(this->GetMaterial())->C(this,gausspointnumber);
			this->Strain=(PreviousStrain+i*StrainDifference);
			//Ανανεώνω τα strains sta gausspoints
			CurrentStress=C*StrainDifference + this->Stress;
//			ValidateStress(PreviousStress,CurrentStress);
			this->Stress=CurrentStress;
		}
		//Ειναι Δσ=C*Δε οπου Δε η μεταβολή των παραμορφώσεων
	}
	if (method==2)
	{
		CurrentStress=fem->GetMaterial(this->GetMaterial())->Stress(this,gausspointnumber);
		this->Stress=CurrentStress;
		this->Strain=CurrentStrain;
	}
	if (update==false)
	{
		this->Strain=PreviousStrain;
		this->Stress=PreviousStress;
	}
	return CurrentStress;

}



ColumnVector Truss2D::InternalForce(bool update,int method)
{
	double c;//cosine of truss angle
	double s;//sine of truss angle
	double l;//length of truss
	double x1,y1,x2,y2;
	x1=fem->GetNode(this->GetNode(1))->GetCoord(1);
	y1=fem->GetNode(this->GetNode(1))->GetCoord(2);
	x2=fem->GetNode(this->GetNode(2))->GetCoord(1);
	y2=fem->GetNode(this->GetNode(2))->GetCoord(2);
	l=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) );
	c=(x2-x1)/l;
	s=(y2-y1)/l;
	int VectorLength=dofs()*GetNumberOfNodes();
	ColumnVector IntForce(VectorLength);
	IntForce=0;//Το διανυσμα των φορτίων στο γενικό σύστημα συντεταγμένων
	double Force=0;//Η αξονική δύναμη της ράβδου
	Force=Force + (fem->GetProperty(this->GetProperty())->GetParameter(1) )* GetStress(1,update,method)(1);
	IntForce(1)=-Force*c;IntForce(2)=-Force*s;IntForce(3)=Force*c;IntForce(4)=Force*s;
	return IntForce;

}

void Truss2D::GeometryToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"\nTRSS2D   ";
	kos<<"\t"<<GetIndex();
	for(int l=1;l<=GetNumberOfNodes();l++)
	{
		kos<<"\t"<<GetNode(l);
	}
	kos<<"\t"<<GetMaterial();
	kos<<"\t"<<GetProperty();
	flush(kos);
	kos.close();

}
