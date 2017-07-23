// IsoparametricElement.h: interface for the IsoparametricElement class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ISOPARAMETRICELEMENT_H__06628341_7E18_11D2_B4AC_444553540000__INCLUDED_)
#define AFX_ISOPARAMETRICELEMENT_H__06628341_7E18_11D2_B4AC_444553540000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Newmat.h"
#include "Newmatap.h"
#include "Statheres.h"
#include "Element.h"
#include "GaussPoint.h"	// Added by ClassView

class IsoparametricElement : public Element  
{
public:
	void SetOrderOfIntegration(int intgorder);
	int GetOrderOfIntegration();
	void GeometryToFile(char* filename);
	ColumnVector ReadStress(int gausspointnumber);
	virtual void StressAndStrainToFile(char* filename,bool principal,int LoadStep);
	virtual double dV(int gausspointnumber);
	virtual double dV(double ksi, double ita);
	ColumnVector InternalForce(bool update,int method);
	ColumnVector GetStress(int gausspointnumber,bool update,int method);
	ColumnVector GetStrain(int gausspointnumber);
	ColumnVector GetStrain(double ksi, double ita);
	virtual Matrix Bita(int gausspointnumber);
	virtual Matrix Bita(double ksi, double ita);
	void SetNodeCoordSize(int nodes);
	virtual Matrix B(int gausspointnumber);
	virtual Matrix B(double ksi, double ita);
	double GetNodeCoord(int nodenumber,int coordnumber);
	//����� ���� coordnumber ������������ ��� nodenumber ������ ��� ���� value
	void SetNodeCoord(int nodenumber,int coordnumber,double value);
	virtual Matrix Jacombian(int gausspointnumber);
	virtual Matrix Jacombian(double ksi, double ita);
	//��������� �������� ��� ������ � ��� gausspoint j
	virtual double N(int nodenumber, int gausspointnumber);
	virtual double N(int nodenumber, double ksi, double ita);
	//��������� ��� ���������� �������� ��� ������ � ��� gausspoint j
	//���� ��������� k
	virtual double dN(int nodenumber, int gausspointnumber, int direction);
	virtual double dN(int nodenumber, double ksi, double ita, int direction);
	//���������� ��� ������ �������� ��� ���������
	virtual Matrix CreateStiffnessMatrix();
	//���������� �� ���� ��� gausspoint �� index gausspointnumber
	//virtual double weight(int gausspointnumber); Not Needed Anymore Implemented As A LIbrary
	//���������� ��� ������������� ��� gausspoint �� index gausspointnumber
	//virtual double* xi(int gausspointnumber);Not Needed Anymore Implemented As A LIbrary
	//� ��������� ���������� ��� ������������� ��� �� ����
	//��� Gausspoints
	//virtual void InitGauss();
	double GetGlobalCoords(int gausspointnumber,int coord);
	IsoparametricElement();
	IsoparametricElement(FEM* fem);
	virtual ~IsoparametricElement();
	 GaussPoint* gausspoints; 

private:
	//� ������� ����� �������� ��� ������� ������������� ��� ������ ��� ���������
	//�� ������� ��� ������ ����� �� ������ ��� �� ������ �� �������������
	Matrix nodecoord;
	int OrderOfIntegration;

};

#endif // !defined(AFX_ISOPARAMETRICELEMENT_H__06628341_7E18_11D2_B4AC_444553540000__INCLUDED_)
