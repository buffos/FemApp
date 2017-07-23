//� �������� ��� ����������� 2 ��� 2-D ����������

#ifndef _STATHERES_
#define _STATHERES_


#define _MATRICES_
#define _DEFINITIONS_
#include "FemIncludes.h"

double Min(double x,double y); //��������� ��� �������
double Max(double x,double y); //�������� ��� �������
ColumnVector GaussSolver(StiffnessMatrix& K,ColumnVector& b);
//���������� ��� ������ ������ ��� ��� ����� �������
ColumnVector PrincipalVectorStrain(ColumnVector& vector);
//���� � ��������� ���������� ��� ������ ��������������� ���
//������������� ��� ��� ������� �� ��� ���� �� ������ ���� angle
//To �������� ������������� ����� �={��,��,���/2}
ColumnVector PrincipalVectorStress(ColumnVector& vector);
Matrix TransposeStrain(double angle);
double sign(double x);
double dot(ColumnVector& X,ColumnVector& Y);
double dot(ColumnVector& X,ColumnVector& Y,int firstn);
int Irreg(const Matrix& X);
void ValidateStress(const ColumnVector& previousstress ,ColumnVector& stress);
double norm(ColumnVector &V,int type,StiffnessMatrix K,int firstn);


#endif