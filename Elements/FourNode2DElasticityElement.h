// FourNode2DElasticityElement.h: interface for the FourNode2DElasticityElement class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FOURNODE2DELASTICITYELEMENT_H__C921B52A_EC1A_11D2_B4AD_BA24EFEC680D__INCLUDED_)
#define AFX_FOURNODE2DELASTICITYELEMENT_H__C921B52A_EC1A_11D2_B4AD_BA24EFEC680D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "IsoparametricElement.h"
#include "Node.h"
#include "Newmat.h"
#include "Newmatap.h"
#include "Statheres.h"
#include "Element.h"
#include "GaussPoint.h"	// Added by ClassView
#include "FEM.h"
#include <math.h>


class FourNode2DElasticityElement : public IsoparametricElement  
{
public:
	void GeometryToFile(char* filename);
	double dV(int gausspointnumber);
	double dV(double ksi , double ita);
	Matrix Bita(int gausspointnumber);
	Matrix Bita(double ksi , double ita);
	Matrix CreateStiffnessMatrix();
	Matrix B(int gausspointnumber);
	Matrix B(double ksi , double ita);
	Matrix Jacombian(int gausspointnumber);
	Matrix Jacombian(double ksi , double ita);
	//Παράγωγος της συναρτησης σχήματος του κόμβου ι στο gausspoint j
	//στην διευθυνση k
	double dN(int nodenumber, int gausspointnumber,int direction);
	double dN(int nodenumber, double ksi , double ita, int direction);
	//Συναρτηση σχήματος του κόμβου ι στο gausspoint j
	double N(int nodenumber, int gausspointnumber);
	double N(int nodenumber, double ksi , double ita);
	//Θέτονται όλα τα απαραίτητα στοιχεία των gausspoints
//	void InitGauss();
	FourNode2DElasticityElement(FEM* fem,int index,int node1,int node2,int node3,int node4,int IntegrationOrder);
	virtual ~FourNode2DElasticityElement();

};

#endif // !defined(AFX_FOURNODE2DELASTICITYELEMENT_H__C921B52A_EC1A_11D2_B4AD_BA24EFEC680D__INCLUDED_)
