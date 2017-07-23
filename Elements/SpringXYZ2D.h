// SpringXYZ2D.h: interface for the SpringXYZ2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SpringXYZ2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
#define AFX_SpringXYZ2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Element.h"
#include "FEM.h"


class SpringXYZ2D : public Element  
{
public:
	void GeometryToFile(char* filename);
	ColumnVector InternalForce(bool update,int method);//Επιστρέφει τις εσωτερικές δυνάμεις ενός στοιχείου 
	//για τις δεδομένες μετακινήσεις που υπάρχουν στους κόμβους
	//Αν update=true τοτε ανανεώνω τις τάσεις και τις παραμορφώσεις στο στοιχείο
	//ενω αν update=false απλα προβλεπω τις εσωτερικές δυναμεις για την δεδομένη μεταβολή μετακινήσεων των κόμβων
	ColumnVector GetStress(int gausspointnumber,bool update,int method);
	//Υπολογίζει την παραμόρφωση του στοιχείου
	ColumnVector GetStrain(int gausspointnumber);
	SpringXYZ2D();
	Matrix CreateStiffnessMatrix();
	SpringXYZ2D(FEM* fem,int index,int dof,int firstnodenumber,int secondnodenumber);
	virtual ~SpringXYZ2D();
private:
	ColumnVector Strain;
	ColumnVector Stress;

};

#endif // !defined(AFX_SpringXYZ2D_H__4225D8E1_F8A3_11D2_B4AD_F2638DAEA37F__INCLUDED_)
