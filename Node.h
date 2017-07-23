// Node.h: interface for the Node class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NODE_H__F4361C0A_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_NODE_H__F4361C0A_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define _MATRICES_
#define _DEFINITIONS_
#define _EXTRAFUNCTIONS_

#include "FemIncludes.h"


class FEM;
class Node  
{
public:
	void GeometryToFile(char* filename);
	int GetNumberOfElements();
	int ElementGet(int elementnumber);
	Node& operator = (Node& n);
	//Επιστρεφει πληροφοριες σχετικα με τον Κομβο
	void info();
	void SetExist(int index,int exists);
	int GetExist(int index);
	double GetDisplacement(int index);
	void SetDisplacement(int index,double value);
	//Τοποθετω στην θεση index την τιμη load
	void SetLoad(int index,double load);
	void ElementAdd(int elementnumber);
	//Αυτή η συνάρτηση διαβάζει απο το αρχείο τον κόμβο index
	//Επιστρέφει true αν το διάβασμα έγινε επιτυχώς
	virtual bool ReadNodeFromFile(char* filename,int index);
	virtual bool WriteNodeToFile(char* filename,int index);
	//ορίζει τον πίνακα με τις συντεταγμένες
	void SetCoord(double* coord);
	//ορίζει την index-συντεταγμένη του πίνακα coordinates
	void SetCoord(double coord,int index);
	//επιστρέφει τον πίνακα με τις συντεταγμένες του κόμβου
//	double* GetCoord();
	//επιστρέφει τη index συντεταγμένη του κόμβου
	double GetCoord(int index);
	//Επιστρέφει τον αυξων αριθμό του κόμβου
	int GetIndex();
	//επιστρέφει το φορτίο για συγκεκριμένο βαθμό ελευθερίας
	double GetLoad(int dof);
	int  GetFixed(int index);
	void SetFixed(int index, int fix);
	Node(FEM* fem);
	Node(FEM* fem,int index,double coord,...);
	FEM* fem;
	virtual ~Node();

private:
	int NumberofElements;
	//Ο αύξων αριθμός του κόμβου
	int index;
	//Ο πίνακας αυτός θα έχει τις συντεταγμένες του κόμβου
	doublevector coordinates;
	//Ο πίνακας αυτός κρατάειτα στοιχεία που συντρέχουν στον κόμβο
	//int είναι ο ακέραιος που μας δείχνει την θέση του στοιχείου στην λίστα που περιέχει τα διάφορα FEM
	ColumnVector ElementList;
	//Κρατάει τα κομβικά Φορτία.Εχει διάσταση ίση με dofnode.
	//Για DIM=2 εχω την εξης αντιστοιχία των βαθμών ελευθερίας
	// μετακίνηση χ,μετακίνηση y ,στροφή φ
	//Για DIM=3 έχω την εξής αντιστοιχία των βαθμών ελευθερίας
	//μετακίνηση χ , μετακίνηση ψ , μετακίνηση z , στροφη αξονα χ =στροφή επιπέδου yz
	//στροφή αξονα y = στροφή επιπέδου χz , στροφή z= στροφή επιπέδου χy
	ColumnVector LoadAtNodes;
	//Ο πίνακας αυτός περιγράφει ποιούς απο τους βαθμούς ελευθερίας υπάρχουν
	intvector exist;
	//Ο πίνακας αυτός περιγράφει ποιούς απο τους βαθμούς ελευθερίας είναι φιξαρισμένοι.
	intvector fixed;
	//Εδω αποθηκεύονται οι μετακινήσεις των κόμβων
	ColumnVector  displacements;
};

#endif // !defined(AFX_NODE_H__F4361C0A_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)