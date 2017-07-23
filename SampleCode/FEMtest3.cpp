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
	appl.numberofelements=3;
	//Δέσμευση Απαραίτητης Μνήμης
	appl.Ylika=new Material*[2];//Δεσμεύω 2 δείκτες για την κλάση Material
	appl.Diatomes=new Property*[2];
	appl.nodes=new Node*[4];
	appl.elements=new Element*[4];//Δεσμεύω 4 δείκτες για την κλάση Element.Πρεπει για κάθε ενα στοιχείο να δεσμεύσω χωριστά μνήμη
	///////////////////////////////////
	//////////Ορισμός Κόμβων///////////
	///////////////////////////////////
	Node ena(1,0.,0.);
	ena.ElementAdd(1);
	ena.ElementAdd(2);
	ena.SetLoad(1,10.0);
	ena.info();
	waitt();
	Node dio(2,5.,0.);
	dio.ElementAdd(1);
	dio.ElementAdd(3);
	dio.info();
	waitt();
	Node tria(3,0.,-8.);
	tria.SetFix(1,1);
	tria.SetFix(2,1);
	tria.SetFix(3,1);
	tria.ElementAdd(2);
	tria.info();
	waitt();
	Node tess(4,5.,-8.);
	tess.SetFix(1,1);
	tess.SetFix(2,1);
	tess.SetFix(3,1);
	tess.ElementAdd(3);
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
	TrussMaterial mat(21*pow(10,7));
	appl.Ylika[0]=&mat;
	BeamProperty prop(0.015,0.0045);
	appl.Diatomes[0]=&prop;
	Beam2D first(&appl,1,1,2);
	first.SetMaterial(1);
	first.SetProperty(1);
	Beam2D second(&appl,2,3,1);
	second.SetMaterial(1);
	second.SetProperty(1);
	Beam2D trd(&appl,3,4,2);
	trd.SetMaterial(1);
	trd.SetProperty(1);
	appl.elements[0]=&first;
	appl.elements[1]=&second;
	appl.elements[2]=&trd;
//	cout<<"\n"<<first.CreateStiffnessMatrix();
//	cout<<"\n"<<second.CreateStiffnessMatrix();
	appl.Initialize();
//	cout<<"\n"<<appl.AssembleLoads();
//	cout<<"\n"<<appl.AssembleStiffnessMatrix()(2,2);
	appl.ElasticSolve(1);
	cout<<"\n"<<second.GetDisplacementOfNodes();
//	cout<<"\n"<<second.ElementReaction();
//	cout<<"\n"<<appl.systemband();
//	cout<<"\n"<<appl.ReactionAtNode(1);
//	cout<<"\n"<<appl.ReactionAtNode(2);
//	cout<<"\n"<<appl.ReactionAtNode(3);
//	cout<<"\n"<<appl.ReactionAtNode(4);
}


void waitt()
{
	int t;
	cout<<"Press any key to continue "<<"\n";
	flush(cout);
	t=_getch();
}