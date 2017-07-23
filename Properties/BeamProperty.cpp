// BeamProperty.cpp: implementation of the BeamProperty class.
//
//////////////////////////////////////////////////////////////////////

#include "BeamProperty.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BeamProperty::BeamProperty(double CrossSection,double Inertia)
{
	this->SetNumberOfParameters(2);
	this->SetParameterSize(2);
	this->SetDescription("Idiotites Diatomis gia Beam");
	this->SetParameter(1,CrossSection);
	this->SetParameter(2,Inertia);
}

BeamProperty::~BeamProperty()
{

}
