// SpringXYZ2D.h: interface for the SpringXYZ2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SpringXYZ2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
#define AFX_SpringXYZ2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Element.h"
#include "FEM.h"


class SpringXYZ2D : public Element  
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
	SpringXYZ2D();
	Matrix CreateStiffnessMatrix();
	SpringXYZ2D(FEM* fem,int index,int dof,int firstnodenumber,int secondnodenumber);
	virtual ~SpringXYZ2D();
private:
	ColumnVector Strain;
	ColumnVector Stress;

};

#endif // !defined(AFX_SpringXYZ2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
