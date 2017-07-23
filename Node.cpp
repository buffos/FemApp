// Node.cpp: implementation of the Node class.
//
//////////////////////////////////////////////////////////////////////
#include "Node.h"
#include "FEM.h"
#include <stdio.h>
#include <stdarg.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Node::Node(FEM* fem)
{
	this->fem = fem;
	int variables = fem->GetVariables();
	for(int i = 1;i<=variables;i++)
	{
		exist.push_back(0);
		fixed.push_back(0);
	}
	this->LoadAtNodes.ReSize(variables);
	this->displacements.ReSize(variables);
	this->NumberofElements=0;
	for(int j=1;j<=variables;j++)
	{
		LoadAtNodes(j)=0;
		displacements(j)=0;
	}


}

Node::Node(FEM* fem,int index,double coord,...)
{
	this->fem = fem;
	int variables = fem->GetVariables();
	for(int k = 1;k<=variables;k++)
	{
		exist.push_back(0);
		fixed.push_back(0);
	}
	this->LoadAtNodes.ReSize(variables);
	this->displacements.ReSize(variables);
	this->NumberofElements=0;
	this->index=index;
	//initialisation
	for(int j=1;j<=variables;j++)
	{
		LoadAtNodes(j)=0;
		displacements(j)=0;
	}
	coordinates.push_back(coord);
	va_list sintetagmenes;	
    va_start( sintetagmenes, coord );
	for(int i=1;i<fem->GetGeometry();i++)
	{
		coordinates.push_back( va_arg( sintetagmenes, double) );
	}
	va_end( sintetagmenes );
}

Node::~Node()
{

}

void Node::ElementAdd(int elementnumber)
{	
	(this->NumberofElements)++;
	ColumnVector Temp(NumberofElements-1);
	if (NumberofElements>1)
	{
		Temp=ElementList;
	}
	ElementList.ReSize(NumberofElements);
	if (NumberofElements>1)
	{
		ElementList.Rows(1,NumberofElements-1)=Temp;
	}
	ElementList(NumberofElements)=(double) elementnumber;
}


double Node::GetCoord(int index)
{
	return this->coordinates[index-1];
}

int Node::GetIndex()
{
	return this->index;
}

double Node::GetLoad(int dof)
{
	return this->LoadAtNodes(dof);
}


void Node::SetCoord(double coord,int index)
{
	this->coordinates[index-1]=coord;

}

void Node::SetCoord(double* coord)
{
	for(int i=0;i<fem->GetGeometry();i++)
	{
		this->coordinates[i]=coord[i];
	}
}

void Node::SetLoad(int index, double load)
{
	this->LoadAtNodes(index)=load;
}

void Node::SetDisplacement(int index, double value)
{
	this->displacements(index)=value;
}

double Node::GetDisplacement(int index)
{
	return this->displacements(index);
}

int Node::GetExist(int index)
{
	return exist[index-1];
}

void Node::SetExist(int index, int exists)
{
	exist[index-1]=exists;
}

int Node::GetFixed(int index)
{
	return fixed[index-1];
}

void Node::SetFixed(int index, int fix)
{
	fixed[index-1]=fix;
}



bool Node::ReadNodeFromFile(char* filename,int index)
{
	return true;
}

bool Node::WriteNodeToFile(char* filename,int index)
{
	return true;
}

void Node::info()
{
	char* temp;
	cout<<"\t\t\t"<<"KOMBOS  "<<GetIndex()<<"\n\n";
	cout<<"Syntetagmenes  ";
	for(int j=1;j<=fem->GetGeometry();j++)
	{
		cout<<"\t\t"<<this->GetCoord(j);
	}
	cout<<"\n";
	for(int i=1;i<=3;i++)
	{
		cout<<"\t\t"<<"BATHMOS ELEUTHERIAS "<<"\t"<<i<<"\n";
		temp=(GetExist(i)==0)?("MH YPARKTOS"):("YPARKTOS");
		cout<<"Bathmos Eleutherias  "<<"\t"<<temp<<"\n";
//		temp=(GetFixed(i)==0)?("Eleutheros"):("Desmeumenos");
		cout<<"Bathmos Eleutherias "<<"\t"<<temp<<"\n";
		cout<<"Fortio  "<<"\t\t"<<GetLoad(i)<<"\n";
		cout<<"Metakinisi  "<<"\t\t"<<GetDisplacement(i)<<"\n\n";
	}


}

Node& Node::operator =(Node& n)
{
	this->fem = n.fem;
	int variables = fem->GetVariables();

	this->LoadAtNodes.ReSize(variables);
	this->displacements.ReSize(variables);

	exist.clear();coordinates.clear();
	for(int j=1;j<=variables;j++)
	{
		exist.push_back(n.exist[j-1]);
		fixed.push_back(n.fixed[j-1]);

	}

	this->LoadAtNodes=n.LoadAtNodes;
	this->displacements=n.displacements;

	for(int i=0;i<fem->GetGeometry();i++)
	{
		coordinates.push_back(n.coordinates[i]);
	}
	index=n.index;
	NumberofElements=n.NumberofElements;
	return *this;
}

int Node::ElementGet(int elementnumber)
{
	return this->ElementList(elementnumber);
}


int Node::GetNumberOfElements()
{
	return (int)this->NumberofElements;
}

void Node::GeometryToFile(char *filename)
{	
	fstream kos;
	kos.open(filename,ios::app);

	int variables = fem->GetVariables();
	int geometry = fem->GetGeometry();

	if (geometry==2)
	{
		kos<<"\nNODE2D  ";
	}
	else
	{
		kos<<"NODE3D  ";
	}
	kos<<"\t"<<GetIndex();
	for(int j=1;j<=geometry;j++)
	{
		kos<<"\t"<<GetCoord(j);
	}
	for(int jj=1;jj<=variables;jj++)
	{
		kos<<"\t"<<GetLoad(jj);
	}
	for(int kk=1;kk<=variables;kk++)
	{
		kos<<"\t"<<GetExist(kk);
	}
	flush(kos);
	kos.close();
}
