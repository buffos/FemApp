// Material.cpp: implementation of the Material class.
//
//////////////////////////////////////////////////////////////////////

#include "Material.h"
#include <stdio.h>
#include <stdarg.h>
#include <iostream.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Material::Material()
{

}

Material::~Material()
{

}

void Material::SetDescription(char* label)
{
	this->Description=label;
}

char* Material::GetDescription()
{
	return this->Description;

}

void Material::SetNumberOfParameters(int numb)
{
	this->NumberOfParameters=numb;
}

int Material::GetNumberOfParameters()
{
	return this->NumberOfParameters;
}

void Material::SetParameter(int index, double value)
{
	this->ParamList[index-1]=value;
}

double Material::GetParameter(int index)
{
	return this->ParamList[index-1];
}

Matrix Material::C(Element* elem,int gausspointnumber)
{
	Matrix C;
	return C;
}

void Material::info()
{
	cout<<"Yliko "<<"\t"<<this->GetDescription()<<"\t";
}

void Material::SetParameterSize(int size)
{
	this->ParamList=new double[size];
}

ColumnVector Material::Stress(Element *elem, int gausspointnumber)
{
	ColumnVector stress;
	return stress;
}
