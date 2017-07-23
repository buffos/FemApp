//Η διάσταση του προβλήματος 2 για 2-D προβλήματα

#ifndef _STATHERES_
#define _STATHERES_


#define _MATRICES_
#define _DEFINITIONS_
#include "FemIncludes.h"

double Min(double x,double y); //Ελαχιστος Δυο Αριθμών
double Max(double x,double y); //Μεγιστος Δυο Αριθμών
ColumnVector GaussSolver(StiffnessMatrix& K,ColumnVector& b);
//Επιστρέφει τις κύριες τάσεις και την γωνία στροφής
ColumnVector PrincipalVectorStrain(ColumnVector& vector);
//Αυτή η συνάρτηση επιστρέφει τον πίνακα μετασχηματισμού των
//παραμορφώσεων απο ενα συστημα σε ενα αλλο με στροφή κατα angle
//To διάνυσμα παραμορφώσεων ειναι ε={εχ,εψ,γχψ/2}
ColumnVector PrincipalVectorStress(ColumnVector& vector);
Matrix TransposeStrain(double angle);
double sign(double x);
double dot(ColumnVector& X,ColumnVector& Y);
double dot(ColumnVector& X,ColumnVector& Y,int firstn);
int Irreg(const Matrix& X);
void ValidateStress(const ColumnVector& previousstress ,ColumnVector& stress);
double norm(ColumnVector &V,int type,StiffnessMatrix K,int firstn);


#endif