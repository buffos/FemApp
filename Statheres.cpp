//////////////////////////////////////////////////////////////////////
///////////////////////Helpfull Functions/////////////////////////////
//////////////////////////////////////////////////////////////////////
//Εδω θα γραφονται οι χρησιμες συναρτήσεις

#include "Statheres.h"
#include "iostream.h"
//#include <fstream.h>
#include "Newmatio.h"

double Min(double x,double y)
{
	return ((x<y) ? x : y);
}

double Max(double x,double y)
{
	return ((x>y) ? x : y);
}

ColumnVector GaussSolver(StiffnessMatrix& K,ColumnVector& b)
{
	return K.i()*b;
}


Matrix TransposeStrain(double angle)
{
	Matrix T(3,3);
	double c,s;
	c=cos(angle);s=sin(angle);
	T(1,1)=c*c;T(1,2)=s*s;T(1,3)=c*s;
	T(2,1)=s*s;T(2,2)=c*c;T(2,3)=-c*s;
	T(3,1)=-2*c*s;T(3,2)=2*c*s;T(3,3)=c*c-s*s;
	return T;
}

double sign(double x)
{
	return x/fabs(x);
}

ColumnVector PrincipalVectorStrain(ColumnVector& vector)
{
	//Μονο για διάνυσμα strain
	//γιατι υποθέτοιυμε οτι
	//vector={εχ,εy,γχψ}
	//οποτε κατα τον υπολογισμό του τανυστη Α
	int VectorLength;
	VectorLength=vector.Nrows();
	ColumnVector principal(3);
	if (VectorLength==3)
	{
		SymmetricMatrix A(2);
		DiagonalMatrix D(2);
		Matrix V(2,2);
		principal=0;
		A(1,1)=vector(1);
		A(2,2)=vector(2);
		A(1,2)=vector(3)/2;
		Jacobi(A,D,V);
		principal(1)=D(1);
		principal(2)=D(2);
		if (vector(1)==0 && vector(2)==0) //cosφ=0 & sinφ!=0
		{
			principal(3) =0;
		}
		else if(V(1,1)==0)
		{
			principal(3) = (V(2,1)==1 )? PI/2 : -PI/2;
		}
		else
		{
			principal(3)=atan( V(2,1)/V(1,1) );
		}
	}
	else //vectorlength = 6
	{

	}
	return principal;
}

ColumnVector PrincipalVectorStress(ColumnVector& vector)
{
	//Μονο για διάνυσμα strain
	//γιατι υποθέτοιυμε οτι
	//vector={εχ,εy,γχψ}
	//οποτε κατα τον υπολογισμό του τανυστη Α
	int VectorLength;
	VectorLength=vector.Nrows();
	ColumnVector principal(3);
	if (VectorLength==3)
	{
		SymmetricMatrix A(2);
		DiagonalMatrix D(2);
		Matrix V(2,2);
		principal=0;
		A(1,1)=vector(1);
		A(2,2)=vector(2);
		A(1,2)=vector(3);
		Jacobi(A,D,V);
		principal(1)=D(1);
		principal(2)=D(2);
		if (vector(1)==0 && vector(2)==0) //cosφ=0 & sinφ!=0
		{
			principal(3) =0;
		}
		else if(V(1,1)==0)
		{
			principal(3) = (V(2,1)==1 )? PI/2 : -PI/2;
		}
		else
		{
			principal(3)=atan( V(2,1)/V(1,1) );
		}
	}
	else //vectorlength = 6
	{

	}
	return principal;
}

double dot(ColumnVector& X,ColumnVector& Y)
{
	ColumnVector R;
	R=X.t()*Y;
	return R(1);
}


double dot(ColumnVector& X,ColumnVector& Y,int firstn)
{
	ColumnVector R;
	R=(X.t()).Columns(1,firstn)*Y.Rows(1,firstn);
	return R(1);
}

int Irreg(const Matrix& X)
{
	int nrow,count=0;
	nrow=X.Nrows();
	for(int i=1;i<=nrow;i++)
	{
		if( X(i,i)<0)
		{
			count++;
		}
	}
	return count;
}

void ValidateStress(const ColumnVector& previousstress ,ColumnVector& stress)
{
	for(int i=1;i<=previousstress.Nrows();i++)
	{
		if(previousstress(i)*stress(i)<0)
		{
			stress(i)=0;
		}
	}
}



double norm(ColumnVector &V,int type,StiffnessMatrix K,int firstn)
{
	double result;
	result = 0.0;
	double averagestiffness;
	averagestiffness = 0.0;
	if (type == 1)
	{
		result = sqrt( (V.Rows(1,firstn)).SumSquare());
		//normal euclidean norm
	}
	else if (type == 2)
	{
		for(int i=1;i<=firstn;i++)
		{
			double c=1/fabs(K(i,i));
			result += c*V(i)*V(i);
//			averagestiffness += K(i,i)*K(i,i);
		}
		result = sqrt(result);
	}
	return result;
}
