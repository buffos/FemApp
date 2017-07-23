// Element.cpp: implementation of the Element class.
//
//////////////////////////////////////////////////////////////////////

#include "Element.h"
#include <math.h>
#include "FEM.h"
#include <fstream.h>
#include "Newmatio.h" // for testing



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Element::Element()
{

}

Element::~Element()
{

}

void Element::ActivateExistingDofs()
{
	for(int i = 1 ; i<=this->fem->GetVariables();i++)
	{
		if(Getdofcode(i)==1)//activate
		{
			for(int j = 1;j<=this->NumberOfNodes;j++)
			{
				this->fem->GetNode(this->GetNode(j))->SetExist(i,1);
			}
		}

	}
}

int Element::GetMaterial()
{
	return this->Material;
}

void Element::SetMaterial(int material)
{
	this->Material=material;
}

int Element::GetProperty()
{
	return this->Property;
}

void Element::SetProperty(int property)
{
	this->Property=property;
}

int Element::GetNode(int nodenumber)
{
	return this->Nodes[nodenumber-1];
}

void Element::SetNode(int internalnodenumber,int externalnodenumber)
{
	this->Nodes[internalnodenumber-1]=externalnodenumber;
}

void Element::SetNumberOfGaussPoints(int number)
{
	this->NumberOfGaussPoints=number;
}

int Element::GetNumberOfGaussPoints()
{
	return this->NumberOfGaussPoints;
}

Element::Element(FEM *fem)
{
	this->fem=fem;
}

void Element::Setdofcode(int position, int dof)
{
	if (dofcode.size()==0)
	{ 
//		for(int i=1;i<=fem->GetVariables();i++)
//		{
//			dofcode.push_back(0);
//		}
		dofcode.resize(fem->GetVariables(),0);
	}

	this->dofcode[position-1]=dof;
}

int Element::Getdofcode(int position)
{
	return this->dofcode[position-1];
}

int Element::GetNumberOfNodes()
{
	return this->NumberOfNodes;
}

Element::SetNumberofNodes(int numberofnodes)
{
	this->NumberOfNodes=numberofnodes;
}

int Element::band()
{
	int bandmax=0;
	int nodei;
	int nodej;
	int temp;
	for(int i=1;i<=this->GetNumberOfNodes();i++)
	{
		for(int j=i;j<=this->GetNumberOfNodes();j++)
		{
			nodei=this->GetNode(i);
			nodej=this->GetNode(j);
			temp=fem->GetVariables()*( abs(nodei-nodej)+1 );
			if(bandmax<temp)
			{
				bandmax=temp;
			}
		}
	}
	return bandmax;
}

ColumnVector Element::ConnectivityVector()
{
	int tempindex=0;
	int i=0;
	int initialdofnumber=0;
	//Βρίσκω πόσους βαθμούς ελευθερίας εχει ο κάθε κόμβος μου
	int dof=this->dofs();
	//
	ColumnVector conn(dof*this->GetNumberOfNodes());
	for(int j=1;j<=this->GetNumberOfNodes();j++)
	{
		//επιστρέφει τον αριθμό του κόμβου στο γενικό σύστημα συντεταγμένων
		i=this->GetNode(j);
		//οπου i-1 οι κόμβοι που προηγήθηκαν του εξεταζόμενου κόμβου
		//(DIM-1)*3 οι βαθμοί ελευθερίας καθε κόμβου ασχετα αν μερικοί είναι ανενεργοί
		//δηλαδή 3 για DIM=2 και 6 για DIM=3
		// και fem->LostDosfs(i) επιστρέφει τους βαθμούς ελευθερίας που ΔΕΝ
		// υπάρχουν απο τον 1 μεχρι και τον i-1 κόμβο
		initialdofnumber=(i-1)*fem->GetVariables()-this->fem->LostDofsVector(i);
		for(int k=1;k<=fem->GetVariables();k++)
		{
			if(this->Getdofcode(k)==1)
			{
				//Αν ο βαθμός ελευθερίας αξιοποιείται απο το στοιχείο
				//τοτε θα υπάρχει και στον πινακα ακαμψίας
				//αρα αυξάνω τον tempindex κατα 1
				//αυξάνω τον initialdofnumber αφου ο βαθμός ελευθερίας
				//ειναι ενεργός και για το σύστημα
				tempindex+=1;
				initialdofnumber+=1;
				conn(tempindex)=initialdofnumber;
			}
			if( (this->Getdofcode(k)==0)  && (fem->GetNode(i)->GetExist(k)==1) )
			{
				//Αν ο βαθμός ελευθερίας δεν αξιοποιείται απο το στοιχέιο
				//αλλα είναι ενεργός στο σύστημα
				//τότε απλά αυξάνω τον initialdofnumber
				initialdofnumber+=1;
			}
		}//end for(k)
	}//enf for(j)

	return conn;
}

int Element::dofs()
{
	int sum=0;
	for(int i=1;i<=fem->GetVariables();i++)
	{
		if(this->Getdofcode(i)==1)
		{
			sum+=1;
		}
	}
	return sum;
}

void Element::SetIndex(int index)
{
	this->index=index;
}

ColumnVector Element::GetStrain(int gausspointnumber)
{
	ColumnVector Temp(1);
	Temp=0;
	return Temp;
}

ColumnVector Element::GetStrain(double ksi, double ita)
{
	ColumnVector Temp(1);
	Temp=0;
	return Temp;
}


ColumnVector Element::GetStress(int gausspointnumber,bool update,int method)
{
	ColumnVector Stress(fem->GetVariables());
	Stress= ( fem->GetMaterial(this->GetMaterial())->Stress(this,gausspointnumber));
	return Stress;
}



Matrix Element::CreateStiffnessMatrix()
{
	Matrix M;
	return M;
}


void Element::SetNodesSize(int size)
{
	this->Nodes=new int[size];
}

ColumnVector Element::GetDisplacementOfNodes()
{
	int dof=dofs()*(this->GetNumberOfNodes());
	int node;
	int tempindex=0;
	ColumnVector U(dof);//Εδω θα φυλάξω τις μετατοπίσεις
	for(int i=1;i<=this->GetNumberOfNodes();i++)
	{
		node=this->GetNode(i);//ο γενικός αριθμός του Node
		for(int j=1;j<=fem->GetVariables();j++)
		{
			if (this->Getdofcode(j)==1)
			{
				tempindex++;
				U(tempindex)=this->fem->GetNode(node)->GetDisplacement(j);
			}//endif
		}//endfor j
	}//enfor i
	return U;

}

ColumnVector Element::InternalForce(bool update,int method)
{
	return ( (this->CreateStiffnessMatrix()) * (this->GetDisplacementOfNodes()) );
}

void Element::StressAndStrainToFile(char* filename,bool principal,int LoadStep)
{
	fstream elem;
	elem.open(filename,ios::app);
	for(int i=1;i<=this->GetNumberOfGaussPoints();i++)
	{
		printstrain(elem,this->index,i,this->GetStrain(i),LoadStep);
		printstress(elem,this->index,i,this->GetStress(i,true,fem->Method),LoadStep);
	}
	elem.close();

}

ColumnVector Element::ReadStress(int gausspointnumber)
{
	ColumnVector Temp;
	return Temp;
}

void Element::GeometryToFile(char *filename)
{

}

int Element::GetIndex()
{
	return this->index;
}

void Element::MeshToGid(fstream& kos)
{
	for(int m=1;m<=GetNumberOfNodes();m++)
	{
		kos.width(11);
		kos<<this->GetNode(m);
	}
}

double Element::GetGlobalCoords(int gausspointnumber,int coord)
{
	double M;
	M=0.0;
	return M;
}