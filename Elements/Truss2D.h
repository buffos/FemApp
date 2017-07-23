// Truss2D.h: interface for the Truss2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRUSS2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
#define AFX_TRUSS2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Element.h"
#include "FEM.h"


class Truss2D : public Element  
{
public:
	void GeometryToFile(char* filename);
	ColumnVector InternalForce(bool update,int method);//���������� ��� ���������� �������� ���� ��������� 
	//��� ��� ��������� ������������ ��� �������� ����� �������
	//�� update=true ���� �������� ��� ������ ��� ��� ������������� ��� ��������
	//��� �� update=false ���� �������� ��� ���������� �������� ��� ��� �������� �������� ������������ ��� ������
	ColumnVector GetStress(int gausspointnumber,bool update,int method);
	//���������� ��� ����������� ��� ���������
	ColumnVector GetStrain(int gausspointnumber);
	Truss2D();
	Matrix CreateStiffnessMatrix();
	Truss2D(FEM* fem,int index,int firstnodenumber,int secondnodenumber);
	virtual ~Truss2D();
private:
	ColumnVector Strain;
	ColumnVector Stress;

};

#endif // !defined(AFX_TRUSS2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
