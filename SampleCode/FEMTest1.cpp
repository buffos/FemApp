#include "newmatap.h"                // need matrix applications
#include "newmat.h"
#include "iostream.h"
//#include "time.h"
#include "Newmatio.h"
#include "FEM.h"
#include "include.h"
#include <conio.h>
#include <ctype.h>


//#ifdef use_namespace
///using namespace NEWMAT;              // access NEWMAT namespace
//#endif

void waitt();
void main()
{
	//Ιnitialisation of FEM
	FEM appl;
	appl.numberofmaterials=1;
	appl.numberofproperties=1;
	appl.numberofNodes=4;
	appl.numberofelements=1;
	//Δέσμευση Απαραίτητης Μνήμης
	appl.Ylika=new Material*[1];//Δεσμεύω 2 δείκτες για την κλάση Material
	appl.Diatomes=new Property*[1];
	appl.nodes=new Node*[4];
	appl.elements=new Element*[1];//Δεσμεύω 4 δείκτες για την κλάση Element.Πρεπει για κάθε ενα στοιχείο να δεσμεύσω χωριστά μνήμη
	///////////////////////////////////
	//////////Ορισμός Κόμβων///////////
	///////////////////////////////////
	Node ena(1,2.,3.);
	ena.SetFix(1,1);
	ena.SetFix(2,1);
	ena.SetExist(3,0);
	ena.ElementAdd(1);
	ena.info();
	waitt();
	Node dio(2,3.,3.);
	dio.SetFix(2,1);
	dio.SetExist(3,0);
	dio.SetLoad(1,-500.0);
	dio.ElementAdd(1);
//	dio.ElementAdd(2);
	dio.info();
	waitt();
	Node tria(3,3.,4.);
	tria.SetExist(3,0);
	tria.SetLoad(2,-5000.0);
	tria.ElementAdd(1);
	tria.info();
	waitt();
	Node tess(4,2.,4.);
	tess.SetFix(1,1);
//	tess.SetFix(2,1);
	tess.SetExist(3,0);
	tess.ElementAdd(1);
	tess.info();
	waitt();
	/// Τοποθέτηση των κόμβων στον πίνακα
	appl.nodes[0]=&ena;
	appl.nodes[1]=&dio;
	appl.nodes[2]=&tria;
	appl.nodes[3]=&tess;
	//////////////////////////////////////
	////////////Τέλος Ορισμού Κομβων//////
	//////////////////////////////////////
	////////////Δημιουργία Υλικού/////////
	//////////////////////////////////////
//	PlaneStress mat(30000,0.23);
	Veccio2DT mat(-30000,-0.002,1);
	appl.Ylika[0]=&mat;
//	cout<<"\n"<<mat.C(0,0);
	PlaneElementProperty prop(1);
	appl.Diatomes[0]=&prop;
	FourNode2DElasticityElement first(&appl,1,1,2,3,4);
	first.SetMaterial(1);
	first.SetProperty(1);
	appl.elements[0]=&first;
//	cout<<"\n"<<first.CreateStiffnessMatrix();
	appl.Initialize();
//	appl.ElasticSolve(1);
//	ColumnVector s(3);
//	cout<<first.Jacombian(1);
//	appl.OrthogonalResidualMethod(1,10,.001);
//	cout<<mat.C(&first,1);
//	cout<<"\n"<<first.GetStress(1,false).t();
//	cout<<"\n"<<first.GetStrain(1).t();
//	cout<<"\n"<<mat.C(&first,1)*first.GetStrain(1);
//	flush(cout);
//	cout<<appl.AssembleInternalLoads();
//	cout<<first.InternalForce();
//	cout<<"\n"<<appl.ReactionAtNode(1);
//	cout<<"\n"<<first.GetStressAndUpdate(2);
//	cout<<"\n"<<first.gausspoints[1].GetStress();
	ColumnVector spri(3);
	ColumnVector s(3);
	s(1)=12;s(2)=-3;s(3)=4;
	spri=PrincipalVectorStrain(s);
	ColumnVector str(3);
	str(1)=5;str(2)=-1;str(3)=0;
	str=TransposeStrain(spri(3)).t()*str;
	str=PrincipalVectorStress(str);
	cout<<spri.t()<<"Strain Vector\n";
	cout<<str.t()<<"Stress Vector\n";
	flush(cout);
	

}


void waitt()
{
	int t;
	cout<<"Press any key to continue "<<"\n";
	flush(cout);
	t=_getch();
}