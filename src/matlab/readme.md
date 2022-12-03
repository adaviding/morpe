# matlab
**This code is filthy.**  You probably don't want to use it.

## status
This code worked in the years 2008-2010, and I have no idea if it works with modern Matlab version.  I have not used or tested this version since 2010.  I do not intend to work on this version in the future.

## summary
This is the original implementation of MoRPE.  I invented MoRPE almost by accident because nothing else seemed to fit the data I collected for my dissertation.  I tried many different things before MoRPE emerged from the chaos, and so this code is very chaotic.  There is plenty of dead code here, and the C code is almost impossible to understand (this was my first time writing anything in C or C++).

## known issues
* Not actively maintained or tested.
* API is quirky, poorly documented, and difficult to use.
* The optimization routine has a known error for cases where the input data is not centered around the origin.
  * To fix this error, compute the median of each category, and center your data around the average median. 

## build and usage
In order for the C files to be callable from Matlab, each file needs to be "mexed" from the Matlab command line.  Here I have already mexed the C files for 32-bit and 64-bit Windows, and so the software will work from a 32-bit or 64-bit Windows machine.  If you want this software to work on a different platform, you might need to mex the files from that platform (but it might work as is... I'm not sure).  For more information, see the Matlab documentation for `mex`.

Matlab Script for mexing files (if you need to):

```
setenv('VS120COMNTOOLS','C:\Program Files (x86)\Microsoft Visual Studio 12.0\Common7\Tools\');
mex Mcl_Accuracy.cpp;
mex Mcl_ConditionalEntropy.cpp;
mex Mcl_ForceMonotonic.cpp;
mex Mcl_MapDv.cpp;
mex Mcl_Poly_CalcDv.cpp;
mex Mcl_Poly_Coeff.cpp;
mex Mcl_Poly_Init.cpp;
mex Mcl_QuantizeDv.cpp;
mex Mcl_QuickSort.cpp;
mex Mcl_RandRotate.cpp;
mex Mcl_RangeLimit.cpp;
```

In this folder, you will find:
* A number of cpp-files.
	* Any cpp-file prefixed by "Mcl_" needs to be mexed (if not already mexed).
	* Any other cpp-file does not need to be mexed.  These files are used in `#include` directives.
* A number of m-files beginning with the `Mcl_` prefix.  These files are the implementation MoRPE (as MCL stands for MoRPE Classifier Library).
* A number of m-files beginning with the `sTest_` prefix.  These are script files that test parts of the MCL, or otherwise demonstrate its usage.
* This example shows how to train a MoRPE classifier and then classify data.
	* sMcl_TemplateForAnalysis.m
* This code also tests the implementation of a MoRPE classifier, and so it should be useful.
	* sTest_Mcl_Poly.m
