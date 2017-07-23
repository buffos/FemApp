// SpringXYZ2D.cpp: implementation of the SpringXYZ2D class.
//
//////////////////////////////////////////////////////////////////////

#include "SpringXYZ2D.h"
#include <math.h>
#include <iostream.h>
#include <fstream.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SpringXYZ2D::~SpringXYZ2D()
{

}

SpringXYZ2D::SpringXYZ2D(FEM* fem,int index,int dof, int firstnodenumber, int secondnodenumber)
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
//	this->Setdofcode(1,0);
//	this->Setdofcode(2,0);
//	this->Setdofcode(3,0);
	this->Setdofcode(dof,1);// το dof είναι ο βαθμός ελευθερίας που ενεργοποιεί το ελατήριο
	this->Strain.ReSize(1);
	this->Stress.ReSize(1);
	this->Strain=0;
	this->Stress=0;
	fem->GetNode(firstnodenumber)->ElementAdd(this->GetIndex() );
	fem->GetNode(secondnodenumber)->ElementAdd(this->GetIndex() );

}

Matrix SpringXYZ2D::CreateStiffnessMatrix()
{
	double ke;//Δυσκαμψία ελατηρίου
	ke=this->fem->GetMaterial(this->GetMaterial())->C(this,0)(1,1);
	SymmetricMatrix k(2);
	k(1,1)=1;
	k(1,2)=-1;
	k(2,2)=1;
	k=ke*k;
	return k;

}

SpringXYZ2D::SpringXYZ2D()
{
	this->SetNodesSize(2);
	this->SetNumberofNodes(2);

}

ColumnVector SpringXYZ2D::GetStrain(int gausspointnumber)
{
	// Η συνάρτηση αυτή επιστρέφει την τιμή u2-u1=u
	//οπου u1=μετακίνηση του κόμβου ένα κατα μήκος του άξονα του ελατηρίου και
	//u2=μετακίνηση του κόμβου δύο κατα μήκος του άξονα του ελατηρίου 
	//Ετσι γνωρίζοντας το u μπορώ να υπολογίσω τόσο το k όσο και την δύναμη F
	//που αντιστοιχεί στην μετακίνηση u
	ColumnVector displacements=GetDisplacementOfNodes();
	ColumnVector L(2);
	L(1)=-1;L(2)=1;
	return L.t()*displacements;
}


ColumnVector SpringXYZ2D::GetStress(int gausspointnumber,bool update,int method)
{
	//Ενημερώνει την τάση στο Gausspoint κατα δεδομένο Fraction (χρησιμεύει στην ολοκήρωση)
	ColumnVector CurrentStress(1);
	CurrentStress=0;
	return CurrentStress;

}



ColumnVector SpringXYZ2D::InternalForce(bool update,int method)
{
	int VectorLength=dofs()*GetNumberOfNodes();
	ColumnVector IntForce(VectorLength);
	IntForce=0;//Το διανυσμα των φορτίων στο γενικό σύστημα συντεταγμένων
	double Force=0;//Η αξονική δύναμη της ράβδου
	Force=Force + fem->GetMaterial(this->GetMaterial())->Stress(this,1 )(1);
	//Στη περίπτωση του ελατηρίου η συνάρτηση επιστρέφει την δύναμη που προκαλεί η παραμόρφωση  strain
	IntForce(1)=-Force;IntForce(2)=Force;
	return IntForce;

}

void SpringXYZ2D::GeometryToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"\nSSPR2D   ";
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
