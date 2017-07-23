#include "newmatap.h"                // need matrix applications
#include "newmat.h"
#include "iostream.h"
#include "time.h"
#include "Newmatio.h"
#include "FEM.h"
#include "include.h"
#include <conio.h>
#include <ctype.h>

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif


void waitt();
void main()
{
	//Ιnitialisation of FEM
	FEM appl;
	appl.numberofmaterials=2;
	appl.numberofproperties=2;
	appl.numberofNodes=4;
	appl.numberofelements=4;
	//Δέσμευση Απαραίτητης Μνήμης
	appl.Ylika=new Material*[2];//Δεσμεύω 2 δείκτες για την κλάση Material
	appl.Diatomes=new Property*[2];
	appl.nodes=new Node*[4];
	appl.elements=new Element*[4];//Δεσμεύω 4 δείκτες για την κλάση Element.Πρεπει για κάθε ενα στοιχείο να δεσμεύσω χωριστά μνήμη
	///////////////////////////////////
	//////////Ορισμός Κόμβων///////////
	//////////////////////////////////
	Node ena(1,0.,0.);
	ena.SetFix(1,1);
	ena.SetFix(2,1);
	ena.SetExist(3,0);
	ena.ElementAdd(1);
	ena.ElementAdd(3);
	ena.info();
	waitt();
	Node dio(2,40.,0.);
	dio.SetFix(2,1);
	dio.SetExist(3,0);
	dio.SetLoad(1,20000.0);
	dio.ElementAdd(1);
	dio.ElementAdd(2);
	dio.info();
	waitt();
	Node tria(3,40.,30.);
	tria.SetExist(3,0);
	tria.SetLoad(2,-25000.0);
	tria.ElementAdd(2);
	tria.ElementAdd(3);
	tria.ElementAdd(4);
	tria.info();
	waitt();
	Node tess(4,0,30.);
	tess.SetFix(1,1);
	tess.SetFix(2,1);
	tess.SetExist(3,0);
	tess.ElementAdd(4);
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
	TrussMaterial mat(29.5*pow(10,6));
	appl.Ylika[0]=&mat;
	SpringMaterial spr(29.5*pow(10,6)/30);
	appl.Ylika[1]=&spr;
	//////////////////////////////////////
	////////Τέλος Δημιουργίας Υλικού//////
	//////////////////////////////////////
	//////////Δημιουργία Ιδιότητας////////
	//////////////////////////////////////
	TrussProperty prop(1);
	appl.Diatomes[0]=&prop;
	//////////////////////////////////////
	///////Τέλος Δημιουργίας Ιδιοτήτων////
	//////////////////////////////////////
	/////////Δημιουργία Στοιχείων/////////
	
	///////////Ορισμός Στοιχείων//////////
	Truss2D t1(&appl,1,1,2);//,t2(&appl,2,3,2);
	Truss2D t3(&appl,3,1,3),t4(&appl,4,3,4);
	Spring2D t2(&appl,2,3,2);
	///////////////////////////////////////
	/////////Θέτω Υλικά και Ιδιότητες//////
	t1.SetMaterial(1);
	t1.SetProperty(1);
	t2.SetMaterial(2);
	t2.SetProperty(1);
	t3.SetProperty(1);
	t3.SetMaterial(1);
	t4.SetMaterial(1);
	t4.SetProperty(1);
	///////////////////////////////////////
	/////////Τα τοποθετώ στον πίνακα///////
	appl.elements[1]=&t2;
	appl.elements[0]=&t1;
	appl.elements[2]=&t3;
	appl.elements[3]=&t4;
	appl.Initialize();
	cout<<"FORTIA "<<"\n"<<appl.AssembleLoads()<<"\n";
	waitt();
//	cout<<"\n"<<appl.elements[0]->CreateStiffnessMatrix();
	cout<<"\n"<<appl.elements[1]->CreateStiffnessMatrix();
//	cout<<"\n"<<appl.Ylika[appl.elements[1]->GetMaterial()-1]->C(appl.elements[1],1);
	flush(cout);
	waitt();
	cout<<"\n"<<appl.AssembleStiffnessMatrix();
	appl.ElasticSolve(1);
//	cout<<"\n"<<appl.ElasticSolve(1);
	cout<<"\n"<<(appl.elements[3]->CreateStiffnessMatrix())*(appl.elements[3]->GetDisplacementOfNodes());
	cout<<"tora \n";
	ColumnVector X=appl.ReactionAtNode(3);
	cout<<"\n"<<t3.GetStress(1,true,appl.Method)<<"\n";
	cout<<"\n"<<t3.InternalForce(true,appl.Method)<<"\n";
	cout<<(t3.CreateStiffnessMatrix()) * (t3.GetDisplacementOfNodes())<<"\n";

}

void waitt()
{
	int t;
	cout<<"Press any key to continue "<<"\n";
	flush(cout);
	t=_getch();
}