#ifndef _STAGECOMMAND_
#define _STAGECOMMAND_

#define _DEFINITIONS_
#include "FemIncludes.h"


struct StageCommand
{
	int       fromnode;
	int       tonode;
	int       step;
	int       dof;//degree of freedom
	double    load;
	int       fromstage;
	int       tostage;
	int       keep; // to keep loads after the tostage....
};

struct FixCommand
{
	int startingstage;
	int endingstage;
	doublevector coeff;
	intvector dofs; // Kombos, bathmos eleutherias, kombos .... kok....
	double result;
	int shiftingtimes; // poses fores tha epanalabw tin eksisosi me allous deiktes...
	int shiftvalue;
};

#endif