#include "FastExp.cpp"
// Already included #include "mex.h"

//	Static memory for quickly evaluating the exp() via table.
static double* FastExp_Table=null;

// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void Exemplar_DeleteStaticMemory(void)
{
	if( FastExp_Table!=null )
		delete FastExp_Table;
}

// This function initializes static memory if necessary.
static void Exemplar_EnsureStaticMemory()
{
	if( FastExp_Table==null )
	{
		FastExp_Table = FastExp_CreateTable();
		mexAtExit(Exemplar_DeleteStaticMemory);
	}
}