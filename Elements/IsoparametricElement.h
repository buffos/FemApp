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
	//θέτει στην coordnumber συντετεγμένη του nodenumber κομβου την τιμή value
	void SetNodeCoord(int nodenumber,int coordnumber,double value);
	virtual Matrix Jacombian(int gausspointnumber);
	virtual Matrix Jacombian(double ksi, double ita);
	//Συναρτηση σχήματος του κόμβου ι στο gausspoint j
	virtual double N(int nodenumber, int gausspointnumber);
	virtual double N(int nodenumber, double ksi, double ita);
	//Παράγωγος της συναρτησης σχήματος του κόμβου ι στο gausspoint j
	//στην διευθυνση k
	virtual double dN(int nodenumber, int gausspointnumber, int direction);
	virtual double dN(int nodenumber, double ksi, double ita, int direction);
	//υπολογίζει τον πίνακα ακαμψίας του στοιχείου
	virtual Matrix CreateStiffnessMatrix();
	//επιστρεφει τα βάρη του gausspoint με index gausspointnumber
	//virtual double weight(int gausspointnumber); Not Needed Anymore Implemented As A LIbrary
	//επιστρεφει τις συντεταγμενες του gausspoint με index gausspointnumber
	//virtual double* xi(int gausspointnumber);Not Needed Anymore Implemented As A LIbrary
	//Η συναρτηση υπολογίζει τις συντεταγμένες και τα βάρη
	//των Gausspoints
	//virtual void InitGauss();
	double GetGlobalCoords(int gausspointnumber,int coord);
	IsoparametricElement();
	IsoparametricElement(FEM* fem);
	virtual ~IsoparametricElement();
	 GaussPoint* gausspoints; 

private:
	//Ο πίνακας αυτός περιέχει τις φυσικές συντεταγμένες των κόμβων του στοιχείου
	//οι γραμμές του πίνακα ειναι οι κόμβοι και οι στήλες οι συντεταγμένες
	Matrix nodecoord;
	int OrderOfIntegration;

};

#endif // !defined(AFX_ISOPARAMETRICELEMENT_H__06628341_7E18_11D2_B4AC_444553540000__INCLUDED_)
