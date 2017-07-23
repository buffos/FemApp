
#include "Proccessing.h"

void printvector(fstream& f,const ColumnVector& V,int LoadStep)
{
	int l=V.Nrows();
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
}


void printvector(fstream& f,const ColumnVector& V,int LoadStep,int VecEntries)
{
	int l=VecEntries;
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
}

void printvector(char* name,const ColumnVector& V,int LoadStep,int VecEntries)
{
	fstream f;
	f.open("test.txt",ios::app);
	f<<name<<"\n";
	int l=VecEntries;
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
	f.close();
}

void printvector(char* name,const ColumnVector& V,int LoadStep)
{
	fstream f;
	f.open("test.txt",ios::app);
	f<<name<<"\n";
	int l=V.Nrows();
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
	f.close();
}


void printmatrix(char* name,Matrix& V,int LoadStep)
{
	fstream f;
	f.open("test.txt",ios::app);
	f<<name<<"\n";
	f.width(5);f.precision(6);f<<LoadStep<<"\n" ;
	for(int i=1;i<=V.Nrows();i++)
	{
		for(int j=1;j<=V.Ncols();j++)
		{
			f.setf(ios::fixed);
			f.width(15);f.precision(6);f<<V(i,j)<<"\t";
		}
		f<<"\n";
	}
	f.close();
}


void printvector(fstream& f,int elementnumber,const ColumnVector& V,int LoadStep)
{
	int l=V.Nrows();
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<elementnumber<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
}

void printvector(fstream& f,const DiagonalMatrix& V,int LoadStep)
{
	int l=V.Nrows();
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i,i)<<"\n";
	}
}


void printstrain(fstream& f,int elementnumber,int gausspoint,const ColumnVector& V,int LoadStep)
{
	int l=V.Nrows();
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<elementnumber<<" ";
		f.width(10);f.precision(6);f<<"strain"<<" ";
		f.width(5);f.precision(6);f<<gausspoint<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
}


void printstress(fstream& f,int elementnumber,int gausspoint,const ColumnVector& V,int LoadStep)
{
	int l=V.Nrows();
	for(int i=1;i<=l;i++)
	{
		f.setf(ios::fixed);
		f.width(5);f.precision(6);f<<LoadStep<<" ";
		f.width(5);f.precision(6);f<<elementnumber<<" ";
		f.width(10);f.precision(6);f<<"stress"<<" ";
		f.width(5);f.precision(6);f<<gausspoint<<" ";
		f.width(5);f.precision(6);f<<i<<" ";
		f.width(15);f.precision(6);f<<V(i)<<"\n";
	}
}
