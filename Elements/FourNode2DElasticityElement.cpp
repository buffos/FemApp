// FourNode2DElasticityElement.cpp: implementation of the FourNode2DElasticityElement class.
//
//////////////////////////////////////////////////////////////////////

#include "FourNode2DElasticityElement.h"
#include "Node.h"
#include "Newmat.h"
#include "Newmatap.h"
#include "Statheres.h"
#include "FEM.h"
#include <math.h>
#include "fstream.h"  //for testing
#include "Newmatio.h" // for testing



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FourNode2DElasticityElement::FourNode2DElasticityElement(FEM* fem,int index,int node1,int node2,int node3,int node4,int IntegrationOrder)
{
	this->fem=fem;
	this->SetNodeCoordSize(4);
	this->SetNumberofNodes(4);
	this->SetNumberOfGaussPoints(pow(IntegrationOrder,2)); // for Order 1 -> 1 GP -> UnderIntegration
	this->SetIndex(index);
	this->SetOrderOfIntegration(IntegrationOrder);
	this->gausspoints=new GaussPoint[pow(IntegrationOrder,2)]; //Allocate Memory for GaussPoints
	for(int i=1;i<=pow(IntegrationOrder,2);i++){gausspoints[i-1].Initialise(3);}
	//Η x συντεταγμένη υπάρχει
	this->Setdofcode(1,1);
	//Η y συντεταγμένη υπάρχει
	this->Setdofcode(2,1);
	//Η φ συντεταγμένη δεν υπάρχει
//	this->Setdofcode(3,0);
	//Θέτω τις φυσικές συντεταγμένες των κόμβων
	this->SetNodeCoord(1,1,-1);
	this->SetNodeCoord(1,2,-1);
	this->SetNodeCoord(2,1,1);
	this->SetNodeCoord(2,2,-1);
	this->SetNodeCoord(3,1,1);
	this->SetNodeCoord(3,2,1);
	this->SetNodeCoord(4,1,-1);
	this->SetNodeCoord(4,2,1);
	//Θέτω τους κόμβους του στοιχείου
	this->SetNodesSize(4);
	this->SetNode(1,node1);
	this->SetNode(2,node2);
	this->SetNode(3,node3);
	this->SetNode(4,node4);
	fem->GetNode(node1)->ElementAdd(this->GetIndex() );
	fem->GetNode(node2)->ElementAdd(this->GetIndex() );
	fem->GetNode(node3)->ElementAdd(this->GetIndex() );
	fem->GetNode(node4)->ElementAdd(this->GetIndex() );

}

FourNode2DElasticityElement::~FourNode2DElasticityElement()
{

}


double FourNode2DElasticityElement::N(int nodenumber, int gausspointnumber)
{
	//Τιμή της Συνάρτησης σχηματος του κόμβου ι στο gausspoint j
	int i=nodenumber;
	int j=gausspointnumber;
	int k = GetOrderOfIntegration();
	ColumnVector GP(3); //gausspoint
	GP = this->fem->GaussPoint(k,2,j);
	return (1/4.)*(1+this->GetNodeCoord(i,1)*GP(1))*(1+this->GetNodeCoord(i,2)*GP(2));
}


double FourNode2DElasticityElement::N(int nodenumber, double ksi , double ita)
{
	//Τιμή της Συνάρτησης σχηματος του κόμβου ι στο gausspoint j
	int i=nodenumber;
	return (1/4.)*(1+this->GetNodeCoord(i,1)* ksi)*(1+this->GetNodeCoord(i,2)* ita);
}


double FourNode2DElasticityElement::dN(int nodenumber, int gausspointnumber, int direction)
{
	int i=nodenumber;
	int j=gausspointnumber;
	int k=direction;
	double temp;
	int l = GetOrderOfIntegration();
	ColumnVector GP(3); //gausspoint
	GP = this->fem->GaussPoint(l,2,j);

	if (k==1)
	{
		temp=(1/4.)*(1+this->GetNodeCoord(i,2)*GP(2))*this->GetNodeCoord(i,1);
	}
	else
	{
		temp=(1/4.)*(1+this->GetNodeCoord(i,1)*GP(1))*this->GetNodeCoord(i,2);
	}
	return temp;
}

double FourNode2DElasticityElement::dN(int nodenumber, double ksi , double ita, int direction)
{
	int i=nodenumber;
	int k=direction;
	double temp;
	if (k==1)
	{
		temp=(1/4.)*(1+this->GetNodeCoord(i,2)* ita)*this->GetNodeCoord(i,1);
	}
	else
	{
		temp=(1/4.)*(1+this->GetNodeCoord(i,1)* ksi)*this->GetNodeCoord(i,2);
	}
	return temp;
}



Matrix FourNode2DElasticityElement::Jacombian(int gausspointnumber)
{
	Matrix J(2,2);
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,4);
	//πίνακας που περιέχει τις συντεταγμένες των κορυφών των nodes
	//κάθε γραμμή και αλλος κόμβος
	Matrix X(4,2);
	//βοηθητική μεταβλητή
	Node temp(fem);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=4;i++)
		{
			DN(j,i)=this->dN(i,gausspointnumber,j);
			X(i,j)=this->fem->GetNode(this->GetNode(i))->GetCoord(j);
		}
	}
	J=DN*X;
	return J;

}


Matrix FourNode2DElasticityElement::Jacombian(double ksi , double ita)
{
	Matrix J(2,2);
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,4);
	//πίνακας που περιέχει τις συντεταγμένες των κορυφών των nodes
	//κάθε γραμμή και αλλος κόμβος
	Matrix X(4,2);
	//βοηθητική μεταβλητή
	Node temp(fem);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=4;i++)
		{
			DN(j,i)=this->dN(i,ksi,ita,j);
			X(i,j)=this->fem->GetNode(this->GetNode(i))->GetCoord(j);
		}
	}
	J=DN*X;
	return J;

}


Matrix FourNode2DElasticityElement::B(int gausspointnumber)
{
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,4);
	Matrix J(2,2);
	J=Jacombian(gausspointnumber);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=4;i++)
		{
			DN(j,i)=this->dN(i,gausspointnumber,j);
		}
	}
	return J.i()*DN;

}

Matrix FourNode2DElasticityElement::B(double ksi , double ita)
{
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,4);
	Matrix J(2,2);
	J=Jacombian(ksi,ita);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=4;i++)
		{
			DN(j,i)=this->dN(i,ksi,ita,j);
		}
	}
	return J.i()*DN;

}

Matrix FourNode2DElasticityElement::CreateStiffnessMatrix()
{
	Matrix Stiffness(8,8);
	Stiffness=0;
	Matrix C(3,3);
	double weights;
	int k = GetOrderOfIntegration();
	ColumnVector GP(3); //gausspoint

	for(int l=1;l<=this->GetNumberOfGaussPoints();l++)
	{
		//weight ειναι το βαρος του gausspoint
		GP = this->fem->GaussPoint(k,2,l);
		weights=GP(3);
		//H παρακάτω σειράεπιστρεφει τον πίνακα C που συνδέει τάσεις με παραμορφώσεις
		C=this->fem->GetMaterial(this->GetMaterial())->C(this,l);
		Stiffness=Stiffness+weights*Bita(l).t()*C*Bita(l)*dV(l);
	}
	return Stiffness;
}


Matrix FourNode2DElasticityElement::Bita(int gausspointnumber)
{
	Matrix Bita(3,8);
	Matrix Bita1(4,8);
	Matrix Temp(2,4);
	Temp=B(gausspointnumber);
	Bita=0;
	Bita1=0;
	Matrix A(3,4);
	A=0;
	A(1,1)=1;
	A(2,4)=1;
	A(3,2)=1;
	A(3,3)=1;
	for(int k=1;k<=this->GetNumberOfNodes();k++)
	{
		Bita1(1,2*k-1)=Temp(1,k);
		Bita1(2,2*k-1)=Temp(2,k);
		Bita1(3,2*k)=Temp(1,k);
		Bita1(4,2*k)=Temp(2,k);
	}
	Bita=A*Bita1;
	return Bita;
}


Matrix FourNode2DElasticityElement::Bita(double ksi , double ita)
{
	Matrix Bita(3,8);
	Matrix Bita1(4,8);
	Matrix Temp(2,4);
	Temp=B(ksi,ita);
	Bita=0;
	Bita1=0;
	Matrix A(3,4);
	A=0;
	A(1,1)=1;
	A(2,4)=1;
	A(3,2)=1;
	A(3,3)=1;
	for(int k=1;k<=this->GetNumberOfNodes();k++)
	{
		Bita1(1,2*k-1)=Temp(1,k);
		Bita1(2,2*k-1)=Temp(2,k);
		Bita1(3,2*k)=Temp(1,k);
		Bita1(4,2*k)=Temp(2,k);
	}
	Bita=A*Bita1;
	return Bita;
}

double FourNode2DElasticityElement::dV(int gausspointnumber)
{
	//te ειναι το πάχος του στοιχείου και βρισκεται στα properties του υλικού
	double te=fem->GetProperty(this->GetProperty())->GetParameter(1);
	return Jacombian(gausspointnumber).LogDeterminant().Value()*te;
}


double FourNode2DElasticityElement::dV(double ksi , double ita)
{
	//te ειναι το πάχος του στοιχείου και βρισκεται στα properties του υλικού
	double te=fem->GetProperty(this->GetProperty())->GetParameter(1);
	return Jacombian(ksi,ita).LogDeterminant().Value()*te;
}


void FourNode2DElasticityElement::GeometryToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"\nQU4N2D   ";
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
