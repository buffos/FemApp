// Element.h: interface for the Element class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ELEMENT_H__F4361C0B_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_ELEMENT_H__F4361C0B_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
//#include "newmat.h"
//#include "newmatap.h"
//#include "statheres.h"
#define _MATRICES_
#define _EXTRAFUNCTIONS_
#include "FemIncludes.h"
#include <fstream.h>




class FEM ;
class Element  
{
public:
	int GetIndex();
	virtual void GeometryToFile(char* filename);
	ColumnVector ReadStress(int gausspointnumber);
	virtual void StressAndStrainToFile(char* filename,bool principal,int Loadstep);
	virtual ColumnVector InternalForce(bool update,int method);
	ColumnVector GetDisplacementOfNodes();
	//������ �� ������� ��� ������ ��� ������� �� Nodes
	void SetNodesSize(int size);
	void SetIndex(int index);
	//� ��������� ���� ���������� ���� ������� ���������� ��� ���������
	//�� �������� �� ���� �����.���������� �� ���� ������� �dofcode ��� ���������
	int dofs();
	//���������� ��� �������� ��� ������������� ���� �������� ������� ����������
	//����� �������� ������� ��� ����������
	void ActivateExistingDofs();//energopoiei tous katallilous bathmous eleutherias analoga me ton dofcode tou programmatos
	ColumnVector ConnectivityVector();
	int band();
	SetNumberofNodes(int numberofnodes);
	int GetNumberOfNodes();
	int Getdofcode(int position);
	void Setdofcode(int position,int dof);
	void MeshToGid(fstream& kos);
	Element(FEM* fem);
	//���������� ��� ������ ��� ��� ������������� �� ������������ GaussPoint
	int GetNumberOfGaussPoints();
	void SetNumberOfGaussPoints(int number);
	//��������� ��� �������������� ������ ��� ������������� �� ������������ ������
	ColumnVector GetStress(int gausspointnumber,bool update,int method);
	virtual ColumnVector GetStrain(int gausspointnumber);
	virtual ColumnVector GetStrain(double ksi, double ita);
	virtual Matrix CreateStiffnessMatrix();
	virtual double GetGlobalCoords(int gausspointnumber,int coord);//Gausspointglobal coords
	void SetNode(int internalnodenumber,int externalnodenumber);
	int GetNode(int nodenumber);
	void SetProperty(int property);
	int GetProperty();
	void SetMaterial(int material);
	int GetMaterial();
	Element();
	virtual ~Element();
	FEM* fem;


private:
	//� ����� ������� ��� ���������
	int index;
	//� ������� ��� ������ ��� ���������
	int NumberOfNodes;
	//O ������� ��� ������ ��� ���������������
	int Material;
	//O ������� ��� ��������� ��� ���������������
	int Property;
	//������ ��� ���� ������� ���������� ��� ������
	// �������� , ���� � ���� ������ ��� ���������
	intvector dofcode;
	//� ������� ��� Gauss Points ��� ���������
	int NumberOfGaussPoints;
	//� ������������ �������� ��� �� Nodes ���� ��� index ��� node ���� ���� �� ��� ������� ������
	int* Nodes;
	//� �������� ������� ���� ������� ��� ��������� ��� master ���������, ���� ������
	// ���� ���������� ���� ����� �� ������� ��� Nodes,��� Elements,��� Materials
	// ������ �� ������� ����� ������ ���� ���� constructor �lement(FEM* fem))

};

#endif // !defined(AFX_ELEMENT_H__F4361C0B_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
