#ifdef _DEFINITIONS_
	#include <math.h>
	#include <vector>
	#include <string>

	#define PI   3.1415926535

	typedef std::vector<int>    intvector;
	typedef std::vector<double> doublevector;
	//typedef BandMatrix StiffnessMatrix;
#endif

#ifdef _MATRICES_
	#include "newmatap.h"                // need matrix applications
	#include "newmat.h"
	typedef SquareMatrix StiffnessMatrix;
#endif

#ifdef _MATRIXIO_
	#include "Newmatio.h" // for testing
#endif


#ifdef _IOHEADERS_
	#include <fstream.h>
	#include "iostream.h" // for testing
#endif

#ifdef _CONTOLPANEL_
	#include "ControlCenter.h"
#endif

#ifdef _PARAMETERSPANEL_
	#include "ParametersCenter.h"
#endif

#ifdef _CONSTRAINS_
	#include "Constrains.h"
#endif

#ifdef _FEMIO_
	#include "Proccessing.h"
#endif

#ifdef _STAGES_
	#include "StageCommand.h"
#endif

#ifdef _EXTRAFUNCTIONS_
	#include "Statheres.h"
#endif

#ifdef _STRINGS_
	#include "Str.h"
#endif



#ifdef _NODE_
	#include "Node.h"
#endif

#ifdef _ELEMENT_
	#include "Element.h"
#endif


#ifdef _MATERIAL_
	#include "Material.h"
	#undef _MATERIAL_ 
#endif

#ifdef _PROPERTY_
	#include "Property.h"
#endif


#ifdef _FEM_
	#include "FEM.h"
#endif



#ifdef _ELEMENTS_
	#include "Element.h"
	#include "IsoparametricElement.h"
	#include "FourNode2DElasticityElement.h"
	#include "EightNode2DElasticityElement.h"
	#include "Truss2D.h"
	#include "Beam2D.h"
	#include "Truss2D.h"
	#include "Spring2D.h"
	#include "SpringXYZ2D.h"
#endif

#ifdef _MATERIALS_
	#include "TrussMaterial.h"
	#include "TrussMaterialS.h"
	#include "TrussPolyonimial.h"
	#include "SpringMaterial.h"
	#include "SpringMaterialS.h"
	#include "Veccio2DT.h"
	#include "Veccio2DS.h"
	#include "Veccio2DModS.h"
	#include "PlaneStress.h"
#endif

#ifdef _PROPERTIES_
	#include "TrussProperty.h"
	#include "BeamProperty.h"
	#include "PlaneElementProperty.h"
#endif


#ifdef _CONSTANTS_
	#define PENALTY_METHOD       1
	#define LAGRANGE_MULTIPLIERS 2
	#define AUGMENTED_LAGRANGE   3
#endif
