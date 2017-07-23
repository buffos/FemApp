#include "newmatap.h"                // need matrix applications
#include "newmat.h"
#include "iostream.h"
#include "time.h"
#include "Newmatio.h"
#include "FEM.h"
#include "include.h"
#include <fstream.h>
#include <conio.h>
#include <ctype.h>


#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

void waitt();
void CreateNodes(double lx,double ly,int Xseg,int Yseg,FEM& fem);
void CreateElements(double lx,double ly,int Xseg,int Yseg,FEM& fem);
void main()
{
	double lx,ly,thickness,load,XSupportLeft,XSupportRight,fc,ec,LoadStep,ErrorTol,SteelA,SteelE;
	int Xsegments,Ysegments,solvemethod,MaxIter;
	int confined;
	char FileOut[30];
	cout<<"\nOutput Filename (e.g Results.txt)  ";
	cin>>FileOut;
	cout<<"\nGive Length X in mm  ";
	cin>>lx;
	cout<<"\nGive Length Y in mm ";
	cin>>ly;
	cout<<"\nGive number of X segments  ";
	cin>>Xsegments;
	cout<<"\nGive number of Y segments  ";
	cin>>Ysegments;
	cout<<"\nGive Thickness of Beam in mm ";
	cin>>thickness;
	cout<<"\nGive Load of Beam in N/mm    ";
	cin>>load;
	cout<<"\nGive SupportLeft Length of Beam in mm ";
	cin>>XSupportLeft;
	cout<<"\nGive SupportRight Length of Beam in mm ";
	cin>>XSupportRight;
	cout<<"\nGive Compressive strength in MPa ";
	cin>>fc;
	cout<<"\nGive Compressive strain in mm/mm  ";
	cin>>ec;
	cout<<"\nConfined Concrete??  1=Yes 0=No  ";
	cin>>confined;
	cout<<"\nSolving Method ??  1=Full N-R 2=Semi N-R 3=Secant 4=O.R ";
	cin>>solvemethod;
	cout<<"\nLoad Step per Iteration  in Newton  ";
	cin>>LoadStep;
	cout<<"\nError Tolerance    ";
	cin>>ErrorTol;
	cout<<"\nMaximum Number of Iterations till Equillibrum per LoadStep ";
	cin>>MaxIter;
	cout<<"\nSteel Cross Section in m2 ";
	cin>>SteelA;
	cout<<"\nSteel Elasticity meter ";
	cin>>SteelE;

	//Ιnitialisation of FEM
	FEM appl;
	appl.numberofmaterials=3;
	appl.numberofproperties=2;
	appl.numberofNodes=(Xsegments+1)*(Ysegments+1);
	appl.numberofelements=Xsegments*Ysegments+Xsegments;
	appl.ElementFileOut=FileOut;
	appl.IterationsDesired=MaxIter;
	//Δέσμευση Απαραίτητης Μνήμης
	appl.Ylika=new Material*[appl.numberofmaterials];
	appl.Diatomes=new Property*[appl.numberofproperties];
	appl.nodes=new Node*[appl.numberofNodes];
	appl.elements=new Element*[appl.numberofelements];
	/////////////////////////////////////////////////
	//////////Ορισμός Κόμβων και στοιχείων///////////
	/////////////////////////////////////////////////
	CreateNodes(lx,ly,Xsegments,Ysegments,appl);
	CreateElements(lx,ly,Xsegments,Ysegments,appl);
	//////////////////////////////////////
	////////////Τέλος Ορισμού Κομβων//////
	//////////////////////////////////////
	////////////Δημιουργία Υλικού/////////
	//////////////////////////////////////
	Veccio2DT mat(fc,ec,confined);
	appl.Ylika[0]=&mat;
	Veccio2DS mat2(fc,ec,confined);
	appl.Ylika[1]=&mat2;
	TrussMaterial mat3(SteelE);
	appl.Ylika[2]=&mat3;
	PlaneElementProperty prop(thickness);
	appl.Diatomes[0]=&prop;
	TrussProperty prop1(SteelA);
	appl.Diatomes[1]=&prop1;
	int matnumber;
	matnumber=(solvemethod==3)? 2 : 1;
	for(int k=1;k<=appl.numberofelements-Xsegments;k++)
	{
		appl.elements[k-1]->SetMaterial(matnumber);
		appl.elements[k-1]->SetProperty(1);
	}
//////////Setting Loads And fixing points///////////////////////
	for(int l=Ysegments+1;l<=appl.numberofNodes;l+=Ysegments+1)
	{
		if (l==Ysegments+1 ||l==appl.numberofNodes)
		{
			appl.nodes[l-1]->SetLoad(2,-load*(lx/Xsegments)/2 );
		}
		else
		{
			appl.nodes[l-1]->SetLoad(2,-load*(lx/Xsegments) );
		}
	}
	int XfixLeft,XfixRight;//Ποσα σημεια καταλαμβάνει η στίρηξη
	XfixLeft=floor(XSupportLeft/(lx/Xsegments));
	XfixRight=floor(XSupportRight/(lx/Xsegments));
	for(int ll=0;ll<=XfixLeft;ll++) 
	{
		appl.nodes[ll*(Ysegments+1)]->SetFix(1,1);
		appl.nodes[ll*(Ysegments+1)]->SetFix(2,1);
	}
	for(int mm=0;mm<=XfixRight;mm++) 
	{
		appl.nodes[appl.numberofNodes-Ysegments-mm*(Ysegments+1)-1]->SetFix(1,1);
		appl.nodes[appl.numberofNodes-Ysegments-mm*(Ysegments+1)-1]->SetFix(2,1);
	}
///////////////////Initialise APP AND SOlve/////////////////////
	appl.Initialize();
	ofstream kos;
	kos.open("stoixeia.txt");
	for(int iot=1;iot<=appl.numberofelements;iot++)
	{
		kos<<"\n"<<iot;
		for(int kk=1;kk<=appl.elements[iot-1]->GetNumberOfNodes();kk++)
		{
			kos<<"\t"<<appl.elements[iot-1]->GetNode(kk);
		}
		kos<<"\t"<<appl.elements[iot-1]->GetMaterial();
		kos<<"\t"<<appl.elements[iot-1]->GetProperty();
		flush(kos);
	}
	for(int iott=1;iott<=appl.numberofNodes;iott++)
	{
		kos<<"\n"<<iott;
		kos<<"\t"<<appl.nodes[iott-1]->GetCoord(1)<<"\t"<<appl.nodes[iott-1]->GetCoord(2);
		kos<<"\t"<<appl.nodes[iott-1]->GetLoad(1)<<"\t"<<appl.nodes[iott-1]->GetLoad(2);
		kos<<"\t"<<appl.nodes[iott-1]->GetFixed(1)<<"\t"<<appl.nodes[iott-1]->GetFixed(2);
		flush(kos);
	}
	kos.close();
	cout<<appl.elements[1]->ReadStress(1);
	appl.NonLinearSolve(solvemethod,1,LoadStep,ErrorTol);
//	appl.ElasticSolve(1);
}


void waitt()
{
	int t;
	cout<<"Press any key to continue "<<"\n";
	flush(cout);
	t=_getch();
}

void CreateNodes(double lx,double ly,int Xseg,int Yseg,FEM& fem)
{
	//Xseg,Yseg σε ποσα κομμάτια κόβω την δοκό κατα Χ και Ψ
	double XStep,YStep;
	XStep=lx/Xseg;
	YStep=ly/Yseg;
	double xi,yi;
	Node* Temp;
	for(int i=0;i<=Xseg;i++)
	{
		for(int j=0;j<=Yseg;j++)
		{
			xi=i*XStep;
			yi=j*YStep;
			Temp=new Node(i*(Yseg+1)+j+1,xi,yi);
			Temp->SetExist(3,0);
			fem.nodes[i*(Yseg+1)+j]=Temp;
		}
	}
}

void CreateElements(double lx,double ly,int Xseg,int Yseg,FEM& fem)
{
	//Xseg,Yseg σε ποσα κομμάτια κόβω την δοκό κατα Χ και Ψ
	double XStep,YStep;
	int Elementnumber=0;
	XStep=lx/Xseg;
	YStep=ly/Yseg;
	double xi=0,yi=0;
	FourNode2DElasticityElement* Temp;
	for(int i=0;i<=Xseg-1;i++)
	{
		for(int j=0;j<=Yseg-1;j++)
		{
			Elementnumber++;
			Temp=new FourNode2DElasticityElement(&fem,Elementnumber,i*(Yseg+1)+j+1,(i+1)*(Yseg+1)+j+1,(i+1)*(Yseg+1)+j+2,i*(Yseg+1)+j+2);
			fem.elements[i*Yseg+j]=Temp;
//			fem.nodes[i*(Yseg+1)+j]->ElementAdd(i*Yseg+j+1);
//			fem.nodes[(i+1)*(Yseg+1)+j]->ElementAdd(i*Yseg+j+1);
//			fem.nodes[(i+1)*(Yseg+1)+j+1]->ElementAdd(i*Yseg+j+1);
//			fem.nodes[i*(Yseg+1)+j+1]->ElementAdd(i*Yseg+j+1);
		}
	}
	Truss2D *Temp2;
	for(int k=0;k<Xseg;k++)
	{
		Elementnumber++;
		Temp2=new Truss2D(&fem,Elementnumber,k*(Yseg+1)+2,(k+1)*(Yseg+1)+2);
		fem.elements[Elementnumber-1]=Temp2;
		fem.nodes[k*(Yseg+1)+1]->ElementAdd(Elementnumber-1);
		fem.nodes[(k+1)*(Yseg+1)+1]->ElementAdd(Elementnumber-1);
		fem.elements[Elementnumber-1]->SetMaterial(3);
		fem.elements[Elementnumber-1]->SetProperty(2);
	}


}
