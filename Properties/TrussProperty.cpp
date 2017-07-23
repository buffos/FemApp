// TrussProperty.cpp: implementation of the TrussProperty class.
//
//////////////////////////////////////////////////////////////////////

#include "TrussProperty.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TrussProperty::TrussProperty(double CrossSection)
{
	this->SetNumberOfParameters(1);
	this->SetParameterSize(1);
	this->SetDescription("Idiotites Diatomis gia Truss");
	this->SetParameter(1,CrossSection);
}

TrussProperty::~TrussProperty()
{

}
