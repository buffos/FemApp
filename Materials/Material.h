// Material.h: interface for the Material class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATERIAL_H__F4361C14_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_MATERIAL_H__F4361C14_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "Newmat.h"
#include "Newmatap.h"
#include "Statheres.h"
#include "Element.h"


class Material  
{
public:
	//Επιστρέφει το διάνυσμα της τάσης δεδομένου του διανύσματος έντασης
	virtual ColumnVector Stress(Element *elem,int gausspointnumber); 
	void SetParameterSize(int size);
	void info();
	//Θα επιστρεφει τον Πινακα που συνδεει τις τασεις σ με τισ παραμορφώσεις ε
	//Θα ειναι ο πυρήνας που θα χειρίζεται την μή γραμμικότητα του υλικού
	//Χρειάζεται να εχει την διεύθυνση του στοιχείου,κυριως για τα μή γραμμικά υλικα,έτσι
	//ώστε να μπορεί να διαβάσει την εντατική κατάσταση τοθ στοιχείου
	virtual Matrix C(Element* elem,int gausspointnumber);
	//Επιστρέφει την index παραμετρο (η αρίθμηση γινεται απο το 1)
	double GetParameter(int index);
	void SetParameter(int index,double value);
	int GetNumberOfParameters();
	void SetNumberOfParameters(int numb);
	char* GetDescription();
	void SetDescription(char* label);
	Material();
	virtual ~Material();

private:
	//Θα περιέχει τις παραμέτρους του υλικου
	//οπως μέτρο ελαστικότητας Ε,συντ.Poisson ν κλπ
	double* ParamList;
	//Θα περιεχει την ονομασία του υλικού πχ "Σκυρόδεμα-Μη γραμμικη αναλυση"
	char* Description;
	int NumberOfParameters;
};

#endif // !defined(AFX_MATERIAL_H__F4361C14_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
