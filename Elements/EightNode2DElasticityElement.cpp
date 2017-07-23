// EightNode2DElasticityElement.cpp: implementation of the EightNode2DElasticityElement class.
//
//////////////////////////////////////////////////////////////////////

#include "EightNode2DElasticityElement.h"
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

EightNode2DElasticityElement::EightNode2DElasticityElement(FEM* fem,int index,int node1,int node2,int node3,int node4,int node5,int node6,int node7,int node8,int IntegrationOrder)
{
	this->fem=fem;
	this->SetNodeCoordSize(8);
	this->SetNumberofNodes(8);
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
	this->SetNodeCoord(1,1,-1.);
	this->SetNodeCoord(1,2,-1.);
	this->SetNodeCoord(2,1,1.);
	this->SetNodeCoord(2,2,-1.);
	this->SetNodeCoord(3,1,1.);
	this->SetNodeCoord(3,2,1.);
	this->SetNodeCoord(4,1,-1.);
	this->SetNodeCoord(4,2,1.);
	// Ενδιάμεση Κόμβοι
	this->SetNodeCoord(5,1,0.);
	this->SetNodeCoord(5,2,-1.);
	this->SetNodeCoord(6,1,1.);
	this->SetNodeCoord(6,2,0.);
	this->SetNodeCoord(7,1,0.);
	this->SetNodeCoord(7,2,1.);
	this->SetNodeCoord(8,1,-1.);
	this->SetNodeCoord(8,2,0.);

	//Θέτω τους κόμβους του στοιχείου
	this->SetNodesSize(8);
	this->SetNode(1,node1);
	this->SetNode(2,node2);
	this->SetNode(3,node3);
	this->SetNode(4,node4);
	this->SetNode(5,node5);
	this->SetNode(6,node6);
	this->SetNode(7,node7);
	this->SetNode(8,node8);

	fem->GetNode(node1)->ElementAdd(this->GetIndex() );
	fem->GetNode(node2)->ElementAdd(this->GetIndex() );
	fem->GetNode(node3)->ElementAdd(this->GetIndex() );
	fem->GetNode(node4)->ElementAdd(this->GetIndex() );
	fem->GetNode(node5)->ElementAdd(this->GetIndex() );
	fem->GetNode(node6)->ElementAdd(this->GetIndex() );
	fem->GetNode(node7)->ElementAdd(this->GetIndex() );
	fem->GetNode(node8)->ElementAdd(this->GetIndex() );


}

EightNode2DElasticityElement::~EightNode2DElasticityElement()
{

}


double EightNode2DElasticityElement::N(int nodenumber, int gausspointnumber)
{
	//Τιμή της Συνάρτησης σχηματος του κόμβου ι στο gausspoint j
	int i=nodenumber;
	int j=gausspointnumber;
	int k = GetOrderOfIntegration();
	double result;
	ColumnVector GP(3); //gausspoint
	GP = this->fem->GaussPoint(k,2,j);
	double Nksi,Nita;// Node coordinates
	Nksi = this->GetNodeCoord(i,1);
	Nita = this->GetNodeCoord(i,2);
	if (i>=1 && i<=4)
	{
		result = (1/4.)*(1+Nksi*GP(1))*(1+Nita*GP(2))*(Nksi*GP(1)+Nita*GP(2) - 1);
	}
	else if (i==5 || i==7)
	{
		result = (1/2.)*(1-GP(1)*GP(1))*(1+Nita*GP(2));
	}
	else if (i==6 || i==8)
	{
		result = (1/2.)*(1+Nksi*GP(1))*(1-GP(2)*GP(2));
	}
	return result;
}


double EightNode2DElasticityElement::N(int nodenumber, double ksi , double ita)
{
	//Τιμή της Συνάρτησης σχηματος του κόμβου ι στο gausspoint j
	int i=nodenumber;
	double result;
	double Nksi,Nita;// Node coordinates
	Nksi = this->GetNodeCoord(i,1);
	Nita = this->GetNodeCoord(i,2);
	if (i>=1 && i<=4)
	{
		result = (1/4.)*(1+Nksi*ksi)*(1+Nita*ita)*(Nksi*ksi+Nita*ita - 1);
	}
	else if (i==5 || i==7)
	{
		result = (1/2.)*(1-ksi*ksi)*(1+Nita*ita);
	}
	else if (i==6 || i==8)
	{
		result = (1/2.)*(1+Nksi*ksi)*(1-ita*ita);
	}
	return result;
}


double EightNode2DElasticityElement::dN(int nodenumber, int gausspointnumber, int direction)
{
	int i=nodenumber;
	int j=gausspointnumber;
	int k=direction;
	int l = GetOrderOfIntegration();
	double result ;
	double Nksi,Nita;// Node coordinates
	ColumnVector GP(3); //gausspoint
	GP = this->fem->GaussPoint(l,2,j);
	Nksi = this->GetNodeCoord(i,1);
	Nita = this->GetNodeCoord(i,2);

	if (i>=1 && i<=4)
	{
		if (k==1)
		{
			result=(1/4.)*(1+Nita*GP(2))*(Nita*GP(2)*Nksi +2*Nksi*Nksi*GP(1));
		}
		else
		{
			result=(1/4.)*(1+Nksi*GP(1))*(Nksi*GP(1)*Nita+2*Nita*Nita*GP(2));
		}
	}
	else if (i==5 || i==7)
	{
		if (k==1)
		{
			result = -(1+Nita*GP(2))*GP(1);
		}
		else
		{
			result = (1/2.)*(1-GP(1)*GP(1))*Nita;
		}
	}
	else if (i==6 || i==8)
	{
		if (k==1)
		{
			result = (1/2.)*(1-GP(2)*GP(2))*Nksi;
		}
		else
		{
			result = -(1+Nksi*GP(1))*GP(2);
		}
	}

	return result;
}

double EightNode2DElasticityElement::dN(int nodenumber, double ksi , double ita, int direction)
{
	int i=nodenumber;
	int k=direction;
	double result ;
	double Nksi,Nita;// Node coordinates
	Nksi = this->GetNodeCoord(i,1);
	Nita = this->GetNodeCoord(i,2);

	if (i>=1 && i<=4)
	{
		if (k==1)
		{
			result=(1/4.)*(1+Nita*ita)*(Nita*ita*Nksi +2*Nksi*Nksi*ksi);
		}
		else
		{
			result=(1/4.)*(1+Nksi*ksi)*(Nksi*ksi*Nita+2*Nita*Nita*ita);
		}
	}
	else if (i==5 || i==7)
	{
		if (k==1)
		{
			result = -(1+Nita*ita)*ksi;
		}
		else
		{
			result = (1/2.)*(1-ksi*ksi)*Nita;
		}
	}
	else if (i==6 || i==8)
	{
		if (k==1)
		{
			result = (1/2.)*(1-ita*ita)*Nksi;
		}
		else
		{
			result = -(1+Nksi*ksi)*ita;
		}
	}

	return result;
}



Matrix EightNode2DElasticityElement::Jacombian(int gausspointnumber)
{
	Matrix J(2,2);
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,8);
	//πίνακας που περιέχει τις συντεταγμένες των κορυφών των nodes
	//κάθε γραμμή και αλλος κόμβος
	Matrix X(8,2);
	//βοηθητική μεταβλητή
	Node temp(fem);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=8;i++)
		{
			DN(j,i)=this->dN(i,gausspointnumber,j);
			X(i,j)=this->fem->GetNode(this->GetNode(i))->GetCoord(j);
		}
	}
	J=DN*X;
	return J;

}


Matrix EightNode2DElasticityElement::Jacombian(double ksi , double ita)
{
	Matrix J(2,2);
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,8);
	//πίνακας που περιέχει τις συντεταγμένες των κορυφών των nodes
	//κάθε γραμμή και αλλος κόμβος
	Matrix X(8,2);
	//βοηθητική μεταβλητή
	Node temp(fem);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=8;i++)
		{
			DN(j,i)=this->dN(i,ksi,ita,j);
			X(i,j)=this->fem->GetNode(this->GetNode(i))->GetCoord(j);
		}
	}
	J=DN*X;
	return J;

}


Matrix EightNode2DElasticityElement::B(int gausspointnumber)
{
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,8);
	Matrix J(2,2);
	Matrix Jinv(2,2);
	J=Jacombian(gausspointnumber);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=8;i++)
		{
			DN(j,i)=this->dN(i,gausspointnumber,j);
		}
	}
	Jinv = J.i();
	return Jinv*DN;

}

Matrix EightNode2DElasticityElement::B(double ksi , double ita)
{
	//πινακαςμε τις μερικές παραγώγους των συναρτήσεων σχήματος
	//κάθε γραμμή περιέχει τις παραγώγους κατα μία κατεύθυνση
	Matrix DN(2,8);
	Matrix J(2,2);
	J=Jacombian(ksi,ita);
	Matrix Jinv(2,2);
	for(int j=1;j<=2;j++)
	{
		for(int i=1;i<=8;i++)
		{
			DN(j,i)=this->dN(i,ksi,ita,j);
		}
	}
	Jinv = J.i();
	return Jinv*DN;
}

Matrix EightNode2DElasticityElement::CreateStiffnessMatrix()
{
	Matrix Stiffness(16,16);
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


Matrix EightNode2DElasticityElement::Bita(int gausspointnumber)
{
	Matrix Bita(3,16);
	Matrix Bita1(4,16);
	Matrix Temp(2,8);
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


Matrix EightNode2DElasticityElement::Bita(double ksi , double ita)
{
	Matrix Bita(3,16);
	Matrix Bita1(4,16);
	Matrix Temp(2,8);
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

double EightNode2DElasticityElement::dV(int gausspointnumber)
{
	//te ειναι το πάχος του στοιχείου και βρισκεται στα properties του υλικού
	double te=fem->GetProperty(this->GetProperty())->GetParameter(1);
	return Jacombian(gausspointnumber).LogDeterminant().Value()*te;
}


double EightNode2DElasticityElement::dV(double ksi , double ita)
{
	//te ειναι το πάχος του στοιχείου και βρισκεται στα properties του υλικού
	double te=fem->GetProperty(this->GetProperty())->GetParameter(1);
	return Jacombian(ksi,ita).LogDeterminant().Value()*te;
}


void EightNode2DElasticityElement::GeometryToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"\nQU8N2D   ";
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
