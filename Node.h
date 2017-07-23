// Node.h: interface for the Node class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NODE_H__F4361C0A_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_NODE_H__F4361C0A_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define _MATRICES_
#define _DEFINITIONS_
#define _EXTRAFUNCTIONS_

#include "FemIncludes.h"


class FEM;
class Node  
{
public:
	void GeometryToFile(char* filename);
	int GetNumberOfElements();
	int ElementGet(int elementnumber);
	Node& operator = (Node& n);
	//���������� ����������� ������� �� ��� �����
	void info();
	void SetExist(int index,int exists);
	int GetExist(int index);
	double GetDisplacement(int index);
	void SetDisplacement(int index,double value);
	//�������� ���� ���� index ��� ���� load
	void SetLoad(int index,double load);
	void ElementAdd(int elementnumber);
	//���� � ��������� �������� ��� �� ������ ��� ����� index
	//���������� true �� �� �������� ����� ��������
	virtual bool ReadNodeFromFile(char* filename,int index);
	virtual bool WriteNodeToFile(char* filename,int index);
	//������ ��� ������ �� ��� �������������
	void SetCoord(double* coord);
	//������ ��� index-������������ ��� ������ coordinates
	void SetCoord(double coord,int index);
	//���������� ��� ������ �� ��� ������������� ��� ������
//	double* GetCoord();
	//���������� �� index ������������ ��� ������
	double GetCoord(int index);
	//���������� ��� ����� ������ ��� ������
	int GetIndex();
	//���������� �� ������ ��� ������������ ����� ����������
	double GetLoad(int dof);
	int  GetFixed(int index);
	void SetFixed(int index, int fix);
	Node(FEM* fem);
	Node(FEM* fem,int index,double coord,...);
	FEM* fem;
	virtual ~Node();

private:
	int NumberofElements;
	//� ����� ������� ��� ������
	int index;
	//� ������� ����� �� ���� ��� ������������� ��� ������
	doublevector coordinates;
	//� ������� ����� ��������� �������� ��� ���������� ���� �����
	//int ����� � �������� ��� ��� ������� ��� ���� ��� ��������� ���� ����� ��� �������� �� ������� FEM
	ColumnVector ElementList;
	//������� �� ������� ������.���� �������� ��� �� dofnode.
	//��� DIM=2 ��� ��� ���� ����������� ��� ������ ����������
	// ���������� �,���������� y ,������ �
	//��� DIM=3 ��� ��� ���� ����������� ��� ������ ����������
	//���������� � , ���������� � , ���������� z , ������ ����� � =������ �������� yz
	//������ ����� y = ������ �������� �z , ������ z= ������ �������� �y
	ColumnVector LoadAtNodes;
	//� ������� ����� ���������� ������ ��� ���� ������� ���������� ��������
	intvector exist;
	//� ������� ����� ���������� ������ ��� ���� ������� ���������� ����� ������������.
	intvector fixed;
	//��� ������������� �� ������������ ��� ������
	ColumnVector  displacements;
};

#endif // !defined(AFX_NODE_H__F4361C0A_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)