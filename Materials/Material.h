// Material.h: interface for the Material class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATERIAL_H__F4361C14_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_MATERIAL_H__F4361C14_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "Newmat.h"
#include "Newmatap.h"
#include "Statheres.h"
#include "Element.h"


class Material  
{
public:
	//���������� �� �������� ��� ����� ��������� ��� ����������� �������
	virtual ColumnVector Stress(Element *elem,int gausspointnumber); 
	void SetParameterSize(int size);
	void info();
	//�� ���������� ��� ������ ��� ������� ��� ������ � �� ��� ������������� �
	//�� ����� � ������� ��� �� ���������� ��� �� ������������ ��� ������
	//���������� �� ���� ��� ��������� ��� ���������,������ ��� �� �� �������� �����,����
	//���� �� ������ �� �������� ��� �������� ��������� ��� ���������
	virtual Matrix C(Element* elem,int gausspointnumber);
	//���������� ��� index ��������� (� �������� ������� ��� �� 1)
	double GetParameter(int index);
	void SetParameter(int index,double value);
	int GetNumberOfParameters();
	void SetNumberOfParameters(int numb);
	char* GetDescription();
	void SetDescription(char* label);
	Material();
	virtual ~Material();

private:
	//�� �������� ��� ����������� ��� ������
	//���� ����� ������������� �,����.Poisson � ���
	double* ParamList;
	//�� �������� ��� �������� ��� ������ �� "���������-�� �������� �������"
	char* Description;
	int NumberOfParameters;
};

#endif // !defined(AFX_MATERIAL_H__F4361C14_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
