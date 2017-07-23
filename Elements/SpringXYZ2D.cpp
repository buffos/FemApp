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
	//��� ������� ���� �� ��������
	this->SetNumberofNodes(2);
	//���������� ���� ��� ������ ���� �������
	this->SetNodesSize(2);
	this->fem=fem;
	//�������� ��� ������ ��� ��������� ��� ��������
	this->SetIndex(index);
	this->SetNode(1,firstnodenumber);
	this->SetNode(2,secondnodenumber);
//	this->Setdofcode(1,0);
//	this->Setdofcode(2,0);
//	this->Setdofcode(3,0);
	this->Setdofcode(dof,1);// �� dof ����� � ������ ���������� ��� ����������� �� ��������
	this->Strain.ReSize(1);
	this->Stress.ReSize(1);
	this->Strain=0;
	this->Stress=0;
	fem->GetNode(firstnodenumber)->ElementAdd(this->GetIndex() );
	fem->GetNode(secondnodenumber)->ElementAdd(this->GetIndex() );

}

Matrix SpringXYZ2D::CreateStiffnessMatrix()
{
	double ke;//��������� ���������
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
	// � ��������� ���� ���������� ��� ���� u2-u1=u
	//���� u1=���������� ��� ������ ��� ���� ����� ��� ����� ��� ��������� ���
	//u2=���������� ��� ������ ��� ���� ����� ��� ����� ��� ��������� 
	//���� ����������� �� u ����� �� ��������� ���� �� k ��� ��� ��� ������ F
	//��� ����������� ���� ���������� u
	ColumnVector displacements=GetDisplacementOfNodes();
	ColumnVector L(2);
	L(1)=-1;L(2)=1;
	return L.t()*displacements;
}


ColumnVector SpringXYZ2D::GetStress(int gausspointnumber,bool update,int method)
{
	//���������� ��� ���� ��� Gausspoint ���� �������� Fraction (���������� ���� ���������)
	ColumnVector CurrentStress(1);
	CurrentStress=0;
	return CurrentStress;

}



ColumnVector SpringXYZ2D::InternalForce(bool update,int method)
{
	int VectorLength=dofs()*GetNumberOfNodes();
	ColumnVector IntForce(VectorLength);
	IntForce=0;//�� �������� ��� ������� ��� ������ ������� �������������
	double Force=0;//� ������� ������ ��� ������
	Force=Force + fem->GetMaterial(this->GetMaterial())->Stress(this,1 )(1);
	//��� ��������� ��� ��������� � ��������� ���������� ��� ������ ��� �������� � �����������  strain
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
