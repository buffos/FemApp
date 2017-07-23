// Spring2D.cpp: implementation of the Spring2D class.
//
//////////////////////////////////////////////////////////////////////

#include "Spring2D.h"
#include <math.h>
#include <fstream.h>
#include <iostream.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Spring2D::~Spring2D()
{

}

Spring2D::Spring2D(FEM* fem,int index, int firstnodenumber, int secondnodenumber)
{
	this->SetNumberOfGaussPoints(1);
	//��� ������� ���� �� ��������
	this->SetNumberofNodes(2);
	//���������� ���� ��� ������ ���� �������
	this->SetNodesSize(2);
	this->fem=fem;
	//�������� ��� ������ ��� ��������� ��� ��������
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

Matrix Spring2D::CreateStiffnessMatrix()
{
	double c;//cosine of truss angle
	double s;//sine of truss angle
	double l;//length of truss
	double ke;//��������� ���������
	ke=this->fem->GetMaterial(this->GetMaterial())->C(this,0)(1,1);
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
	k=ke*k;
	return k;

}

Spring2D::Spring2D()
{
	this->SetNodesSize(2);
	this->SetNumberofNodes(2);

}

ColumnVector Spring2D::GetStrain(int gausspointnumber)
{
	// � ��������� ���� ���������� ��� ���� u2-u1=u
	//���� u1=���������� ��� ������ ��� ���� ����� ��� ����� ��� ��������� ���
	//u2=���������� ��� ������ ��� ���� ����� ��� ����� ��� ��������� 
	//���� ����������� �� u ����� �� ��������� ���� �� k ��� ��� ��� ������ F
	//��� ����������� ���� ���������� u
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
	return L.t()*displacements;
}


ColumnVector Spring2D::GetStress(int gausspointnumber,bool update,int method)
{
	//���������� ��� ���� ��� Gausspoint ���� �������� Fraction (���������� ���� ���������)
	ColumnVector CurrentStress(1);
	CurrentStress=0;
	return CurrentStress;

}



ColumnVector Spring2D::InternalForce(bool update,int method)
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
	IntForce=0;//�� �������� ��� ������� ��� ������ ������� �������������
	double Force=0;//� ������� ������ ��� ������
	Force=Force + fem->GetMaterial(this->GetMaterial())->Stress(this,1 )(1);
	//��� ��������� ��� ��������� � ��������� ���������� ��� ������ ��� �������� � �����������  strain
	IntForce(1)=-Force*c;IntForce(2)=-Force*s;IntForce(3)=Force*c;IntForce(4)=Force*s;
	return IntForce;

}

void Spring2D::GeometryToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"\nSPRG2D   ";
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
