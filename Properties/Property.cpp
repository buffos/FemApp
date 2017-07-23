// Property.cpp: implementation of the Property class.
//
//////////////////////////////////////////////////////////////////////

#include "Property.h"
#include "include.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Property::Property()
{

}

Property::~Property()
{

}

double Property::GetParameter(int index)
{
	return this->ParamList[index-1];
}

void Property::SetParameter(int index, double value)
{
	this->ParamList[index-1]=value;
}

int Property::GetNumberOfParameters()
{
	return this->NumberOfParameters;
}

char* Property::GetDescription()
{
	return this->Description;
}

void Property::SetDescription(char *label)
{
	this->Description=label;
}

void Property::SetNumberOfParameters(int number)
{
	this->NumberOfParameters=number;
}

void Property::SetParameterSize(int size)
{
	this->ParamList=new double[size];
}
