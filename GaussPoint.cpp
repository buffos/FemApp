// GaussPoint.cpp: implementation of the GaussPoint class.
//
//////////////////////////////////////////////////////////////////////

#include "GaussPoint.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GaussPoint::GaussPoint()
{

}

GaussPoint::~GaussPoint()
{

}

void GaussPoint::Initialise(int size)
{
	stress.ReSize( size );
	strain.ReSize( size );
	strain=0;
	stress=0;

}

ColumnVector GaussPoint::GetStrain()
{
	return strain;
}

ColumnVector GaussPoint::GetStress()
{
	return stress;
}

void GaussPoint::SetStrain(ColumnVector Strain)
{
	this->strain=Strain;
}

void GaussPoint::SetStress(ColumnVector Stress)
{
	this->stress=Stress;
}
