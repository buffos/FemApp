// PlaneElementProperty.cpp: implementation of the PlaneElementProperty class.
//
//////////////////////////////////////////////////////////////////////

#include "PlaneElementProperty.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlaneElementProperty::PlaneElementProperty(double thickness)
{
	this->SetParameterSize(1);
	this->SetNumberOfParameters(1);
	this->SetDescription("Idiotites Diatomis gia Epipeda stoixeia");
	this->SetParameter(1,thickness);

}

PlaneElementProperty::~PlaneElementProperty()
{

}
