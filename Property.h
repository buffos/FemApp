// Property.h: interface for the Property class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PROPERTY_H__F4361C15_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
#define AFX_PROPERTY_H__F4361C15_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
//#include <afx.h>
//#include "FemIncludes.h"


class Property  
{
public:
	void SetParameterSize(int size);
	void SetNumberOfParameters(int number);
	void SetDescription(char* label);
	char* GetDescription();
	int GetNumberOfParameters();
	void SetParameter(int index,double value);
	double GetParameter(int index);
	Property();
	virtual ~Property();

private:
	//Θα περιέχει τις ιδιοτητες του υλικου
	//οπως ροπή αδράνειας Ι,εμβαδον διατομής Α κλπ
	double* ParamList;
	//Θα περιεχει την ονομασία της διατομής"
	char* Description;
	int NumberOfParameters;
};

#endif // !defined(AFX_PROPERTY_H__F4361C15_5EA4_11D2_B4AC_9582EC2F2A4D__INCLUDED_)
