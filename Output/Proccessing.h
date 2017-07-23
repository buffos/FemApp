#ifndef _FEMPROCCESSING_
#define _FEMPROCCESSING_


#define _IOHEADERS_
#define _MATRICES_
#include "FemIncludes.h"

void printvector(fstream& f,const ColumnVector& V,int LoadStep);
void printvector(fstream& f,const ColumnVector& V,int LoadStep,int VecEntries);
void printvector(fstream& f,int elementnumber,const ColumnVector& V,int LoadStep);
void printvector(fstream& f,const DiagonalMatrix& V,int LoadStep);
void printstrain(fstream& f,int elementnumber,int gausspoint,const ColumnVector& V,int LoadStep);
void printstress(fstream& f,int elementnumber,int gausspoint,const ColumnVector& V,int LoadStep);
void printvector(char* name,const ColumnVector& V,int LoadStep,int VecEntries);
void printvector(char* name,const ColumnVector& V,int LoadStep);
void printmatrix(char* name,Matrix& V,int LoadStep);

#endif