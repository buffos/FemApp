// Element.h: interface for the Element class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ELEMENT_H__F4361C0B_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_ELEMENT_H__F4361C0B_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
//#include "newmat.h"
//#include "newmatap.h"
//#include "statheres.h"
#define _MATRICES_
#define _EXTRAFUNCTIONS_
#include "FemIncludes.h"
#include <fstream.h>




class FEM ;
class Element  
{
public:
	int GetIndex();
	virtual void GeometryToFile(char* filename);
	ColumnVector ReadStress(int gausspointnumber);
	virtual void StressAndStrainToFile(char* filename,bool principal,int Loadstep);
	virtual ColumnVector InternalForce(bool update,int method);
	ColumnVector GetDisplacementOfNodes();
	//Ορίζει το μέγεθος του πίνακα που κρατάει τα Nodes
	void SetNodesSize(int size);
	void SetIndex(int index);
	//Η συνάρτηση αυτή επιστρέφει τους βαθμούς ελευθερίας που αξιοποιεί
	//το στοιχείο σε κάθε κόμβο.Προυπόθεση να εχει οριστεί οdofcode του στοιχείου
	int dofs();
	//Επιστρέφει ενα διάνυσμα που αντίστοιχίζει τους τοπικούς βαθμούς ελευθερίας
	//στούς γενικούς βαθμούς του συστήματος
	void ActivateExistingDofs();//energopoiei tous katallilous bathmous eleutherias analoga me ton dofcode tou programmatos
	ColumnVector ConnectivityVector();
	int band();
	SetNumberofNodes(int numberofnodes);
	int GetNumberOfNodes();
	int Getdofcode(int position);
	void Setdofcode(int position,int dof);
	void MeshToGid(fstream& kos);
	Element(FEM* fem);
	//υπολογίζει τις τάσεις και τις παραμορφώσεις σε συγκεκριμένο GaussPoint
	int GetNumberOfGaussPoints();
	void SetNumberOfGaussPoints(int number);
	//επιτρέφει τις αναπτυσσόμενες τάσεις και παραμορφώσεις σε συγκεκριμένο σημείο
	ColumnVector GetStress(int gausspointnumber,bool update,int method);
	virtual ColumnVector GetStrain(int gausspointnumber);
	virtual ColumnVector GetStrain(double ksi, double ita);
	virtual Matrix CreateStiffnessMatrix();
	virtual double GetGlobalCoords(int gausspointnumber,int coord);//Gausspointglobal coords
	void SetNode(int internalnodenumber,int externalnodenumber);
	int GetNode(int nodenumber);
	void SetProperty(int property);
	int GetProperty();
	void SetMaterial(int material);
	int GetMaterial();
	Element();
	virtual ~Element();
	FEM* fem;


private:
	//Ο αύξων αριθμός του στοιχείου
	int index;
	//Ο αριθμός των κόμβων του στοιχείου
	int NumberOfNodes;
	//O Αριθμός του Υλικού που χρησιμοποιείται
	int Material;
	//O Αριθμός της Ιδιότητας που χρησιμοποιείται
	int Property;
	//Ποιούς απο τους βαθμούς ελευθερίας που γενικά
	// υπάρχουν , εχει ο κάθε κόμβος του στοιχείου
	intvector dofcode;
	//Ο αριθμός των Gauss Points του στοιχείου
	int NumberOfGaussPoints;
	//Ο Πίνακαςαυτός περιέχει ΟΧΙ τα Nodes αλλα τον index του node έτσι ώστε να τον βρίσκει εύκολα
	int* Nodes;
	//ο παρακάτω δείκτης εχει πάντοτε την διεύθηνση της master εφαρμογής, εκει δηλαδή
	// οπου κρατούνται στην μνήμη οι πίνακες των Nodes,των Elements,των Materials
	// πρεπει να διδεται σαυτό αρχική τιμή στον constructor Εlement(FEM* fem))

};

#endif // !defined(AFX_ELEMENT_H__F4361C0B_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
