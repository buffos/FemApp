// IsoparametricElement.cpp: implementation of the IsoparametricElement class.
//
//////////////////////////////////////////////////////////////////////

#include "IsoparametricElement.h"
#include "Element.h"
#include "FEM.h"
#include "iostream.h"  //for testing
#include "Newmatio.h"  //for testing
#include "fstream.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

IsoparametricElement::IsoparametricElement()
{

}

IsoparametricElement::~IsoparametricElement()
{

}

IsoparametricElement::IsoparametricElement(FEM* fem):Element(fem)
{

}


void IsoparametricElement::SetNodeCoord(int nodenumber,int coordnumber,double value)
{
	this->nodecoord(nodenumber,coordnumber)=value;
}

double IsoparametricElement::GetNodeCoord(int nodenumber, int coordnumber)
{
	return this->nodecoord(nodenumber,coordnumber);
}

//double IsoparametricElement::weight(int gausspointnumber)
//{
//	return this->gausspoints[gausspointnumber-1].GetWeight();
//}

double IsoparametricElement::dN(int nodenumber, int gausspointnumber, int direction)
{
	return 0;
}

double IsoparametricElement::dN(int nodenumber, double ksi , double ita, int direction)
{
	return 0;
}

double IsoparametricElement::N(int nodenumber, int gausspointnumber)
{
	return 0;
}

double IsoparametricElement::N(int nodenumber,double ksi , double ita)
{
	return 0;
}

Matrix IsoparametricElement::B(int gausspointnumber)
{
	Matrix M;
	return M;
}

Matrix IsoparametricElement::B(double ksi , double ita)
{
	Matrix M;
	return M;
}

Matrix IsoparametricElement::CreateStiffnessMatrix()
{
	Matrix M;
	return M;
}


Matrix IsoparametricElement::Jacombian(int gausspointnumber)
{
	Matrix M;
	return M;
}

Matrix IsoparametricElement::Jacombian(double ksi , double ita)
{
	Matrix M;
	return M;
}


//double* IsoparametricElement::xi(int gausspointnumber)
//{
//	return 0;
//}

//void IsoparametricElement::InitGauss()
//{
//	
//}

void IsoparametricElement::SetNodeCoordSize(int nodes)
{
	this->nodecoord.ReSize(nodes,fem->GetGeometry());
	this->nodecoord=0;
}

Matrix IsoparametricElement::Bita(int gausspointnumber)
{
	Matrix Temp;
	return Temp;
}

Matrix IsoparametricElement::Bita(double ksi , double ita)
{
	Matrix Temp;
	return Temp;
}

ColumnVector IsoparametricElement::GetStrain(int gausspointnumber)
{
	ColumnVector Displacement;
	Displacement = this->GetDisplacementOfNodes();
	return Bita(gausspointnumber)*Displacement;

}

ColumnVector IsoparametricElement::GetStrain(double ksi , double ita)
{
	ColumnVector Displacement;
	Displacement = this->GetDisplacementOfNodes();
	return Bita(ksi,ita)*Displacement;

}


ColumnVector IsoparametricElement::GetStress(int gausspointnumber,bool update,int method)
{
	//Ενημερώνει την τάση στο Gausspoint κατα δεδομένο Fraction (χρησιμεύει στην ολοκήρωση)
	ColumnVector CurrentStrain;
	ColumnVector CurrentStress;
	ColumnVector PreviousStrain;
	ColumnVector PreviousStress;
	Matrix C;
	CurrentStrain=GetStrain(gausspointnumber);
	PreviousStrain=gausspoints[gausspointnumber-1].GetStrain();
	PreviousStress=gausspoints[gausspointnumber-1].GetStress();
	if (method==1)
	{
		ColumnVector StrainDifference=(1.0/fem->NumberofItegrationLoops)*(CurrentStrain-PreviousStrain);
		for(int i=1;i<=fem->NumberofItegrationLoops;i++)
		{
			C=this->fem->GetMaterial(this->GetMaterial())->C(this,gausspointnumber);
			gausspoints[gausspointnumber-1].SetStrain(PreviousStrain+i*StrainDifference);
			//Ανανεώνω τα strains sta gausspoints
			CurrentStress=C*StrainDifference + gausspoints[gausspointnumber-1].GetStress();
			ValidateStress(PreviousStress,CurrentStress);
			gausspoints[gausspointnumber-1].SetStress(CurrentStress);
		}
		//Ειναι Δσ=C*Δε οπου Δε η μεταβολή των παραμορφώσεων
	}
	if (method==2)
	{
		CurrentStress=fem->GetMaterial(this->GetMaterial())->Stress(this,gausspointnumber);
		gausspoints[gausspointnumber-1].SetStress(CurrentStress);
		gausspoints[gausspointnumber-1].SetStrain(CurrentStrain);
	}
	if (update==false)
	{
		gausspoints[gausspointnumber-1].SetStrain(PreviousStrain);
		gausspoints[gausspointnumber-1].SetStress(PreviousStress);
	}
	return CurrentStress;
}

ColumnVector IsoparametricElement::InternalForce(bool update,int method)
{
	int VectorLength=dofs()*GetNumberOfNodes();
	ColumnVector IntForce(VectorLength);
	IntForce=0;
	int l = GetOrderOfIntegration();
	ColumnVector GP(3); //this need to change for 3D isoparametric Elements

	for(int i=1;i<=GetNumberOfGaussPoints();i++)
	{
		GP = this->fem->GaussPoint(l,2,i);
		double weights=GP(3);
		IntForce=IntForce + weights * Bita(i).t() * GetStress(i,update,method) * dV(i);
	}//endfor j
	return IntForce;
}

double IsoparametricElement::dV(int gausspointnumber)
{
	return 0.0;
}

double IsoparametricElement::dV(double ksi , double ita)
{
	return 0.0;
}

void IsoparametricElement::StressAndStrainToFile(char* filename,bool principal,int LoadStep)
{
	fstream elem;
	elem.open(filename,ios::app);
	for(int i=1;i<=this->GetNumberOfGaussPoints();i++)
	{
		if(principal == true )
		{
			printstrain(elem,GetIndex(),i,PrincipalVectorStrain(this->GetStrain(i)),LoadStep);
			printstress(elem,GetIndex(),i,PrincipalVectorStress(this->GetStress(i,true,fem->Method)),LoadStep);
		}
		else
		{
			printstrain(elem,GetIndex(),i,this->GetStrain(i),LoadStep);
			printstress(elem,GetIndex(),i,this->GetStress(i,true,fem->Method),LoadStep);
		}
	}
	elem.close();
}

ColumnVector IsoparametricElement::ReadStress(int gausspointnumber)
{
	return this->gausspoints[gausspointnumber-1].GetStress();
}

void IsoparametricElement::GeometryToFile(char *filename)
{

}

int IsoparametricElement::GetOrderOfIntegration()
{
	return this->OrderOfIntegration;
}

void IsoparametricElement::SetOrderOfIntegration(int intgorder)
{
	this->OrderOfIntegration = intgorder; 
}

double IsoparametricElement::GetGlobalCoords(int gausspointnumber,int coord)
{
	double thecoord=0;
	for(int i=1;i<=this->GetNumberOfNodes();i++)
	{
		thecoord += this->N(i,gausspointnumber)*this->fem->GetNode(this->GetNode(i))->GetCoord(coord);

	}
	return thecoord;

}