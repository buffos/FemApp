// Beam2D.cpp: implementation of the Beam2D class.
//
//////////////////////////////////////////////////////////////////////

#include "Beam2D.h"
#include <fstream.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Beam2D::Beam2D()
{

}

Beam2D::~Beam2D()
{

}

Beam2D::Beam2D(FEM *fem, int index, int firstnodenumber, int secondnodenumber)
{
	//��� ������� ���� �� ��������
	this->SetNumberofNodes(2);
	//���������� ���� ��� ������ ���� �������
	this->SetNodesSize(2);
	this->fem=fem;
	//�������� ��� ������ ��� ��������� ��� ��������
	this->SetNumberOfGaussPoints(1);
	this->SetIndex(index);
	this->SetNode(1,firstnodenumber);
	this->SetNode(2,secondnodenumber);
	this->Setdofcode(1,1);
	this->Setdofcode(2,1);
	this->Setdofcode(3,1);
	fem->GetNode(firstnodenumber)->ElementAdd(this->GetIndex() );
	fem->GetNode(secondnodenumber)->ElementAdd(this->GetIndex() );

}

Matrix Beam2D::CreateStiffnessMatrix()
{
	double c;//cosine of beam angle
	double s;//sine of beam angle
	double l;//length of beam
	double A;//cross section of beam
	double E;//����� �������������
	double I;//���� ��������� �����
	E=this->fem->GetMaterial(this->GetMaterial())->C(this,1)(1,1);
	A=this->fem->GetProperty(this->GetProperty())->GetParameter(1);
	I=this->fem->GetProperty(this->GetProperty())->GetParameter(2);
	double x1,y1,x2,y2;
	x1=fem->GetNode(this->GetNode(1))->GetCoord(1);
	y1=fem->GetNode(this->GetNode(1))->GetCoord(2);
	x2=fem->GetNode(this->GetNode(2))->GetCoord(1);
	y2=fem->GetNode(this->GetNode(2))->GetCoord(2);
	l=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) );
	c=(x2-x1)/l;
	s=(y2-y1)/l;
	Matrix L(6,6);//Transformation matrix
	L=0;//initialisation
	L(1,1)=c;L(2,2)=c;L(3,3)=1;L(4,4)=c;L(5,5)=c;L(6,6)=1;
	L(1,2)=s;L(2,1)=-s;L(4,5)=s;L(5,4)=-s;
	SymmetricMatrix k(6);
	k=0;
	k(1,1)=E*A/l;
	k(1,4)=-E*A/l;
	k(2,2)=12*E*I/(l*l*l);
	k(2,3)=6*E*I/(l*l);
	k(2,5)=-12*E*I/(l*l*l);
	k(2,6)=6*E*I/(l*l);
	k(3,3)=4*E*I/l;
	k(3,6)=2*E*I/l;
	k(3,5)=-6*E*I/(l*l);
	k(4,4)=E*A/l;
	k(5,5)=12*E*I/(l*l*l);
	k(5,6)=-6*E*I/(l*l);
	k(6,6)=4*E*I/l;
	L=L.t()*k*L;//��� ����� k=L.t()*k*L ����� � k ����� ����������� ��� ���� ��� ���������� ��� �������
	//��� ����� �� ��������� �.� k*L � ������� ����� ����� �� ����������� ����� ������ ������ ������
	return L;

}

void Beam2D::GeometryToFile(char *filename)
{
	fstream kos;
	kos.open(filename,ios::app);
	kos<<"\nBEAM2D   ";
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
