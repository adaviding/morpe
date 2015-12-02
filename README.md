#MoRPE:  A Probabilistic Classifier
The classification problem is commonly encountered when a finite sample of data is leveraged to determine the probabilistic relationship between a category label  *c*  and a multivariate coordinate  **x**   for an entire population.  Solving this problem requires approximating the optimal classifier, a function of  **x**   that evaluates the conditional probability  *p*(*c*¦**x**)  in an optimal way.  The paper *Ing_Morpe.docx* introduces MoRPE, a method for approximating optimal classifiers.  MoRPE is a machine learning method for probabilistic classification based on **Mo**notonic **R**egression of a **P**olynomial **E**xpansion.  It is conceptually related to Fisher’s Quadratic Discriminant (Fisher 1936) and Kernel Machines.  MoRPE has the ability to approximate an optimal classifier with remarkable precision in common scenarios.

As an algorithm, MoRPE is one of the best classifiers available in the public domain (perhaps *the* best).  Unfortunately, I have not had time to create a stable release for others to use.  The versions posted here are incomplete, difficult to use, and/or buggy.

#Academic Citation
### [Release 0.1](https://github.com/adaviding/Morpe/releases/tag/0.1) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13235.svg)](http://dx.doi.org/10.5281/zenodo.13235)

Ing, AD (2014) MoRPE:  A Probabilistic Classifier.  *GitHub*.  doi:10.5281/zenodo.13235.  Retrieved from https://github.com/adaviding/Morpe/releases/tag/0.1.

#Project Status
Work on the open source version is currently halted because I am currently developing MoRPE for a private entity.

The Matlab codebase is functioning, but it is quirky and difficult for others to use.  (See the Matlab notes below.)  The optimization routine has a known error for cases where the input data is not centered around the origin.  To fix this error, compute the median of each category, and center your data around the average median.

The C# codebase is halfway complete.  I have stopped working on this version because of a conflict with my employment contract.  (See the C# notes below.)

#Recommended Usage
MoRPE should be used in scenarios when the categories are partially overlapping (such that perfect classification performance is impossible).  If the categories are fully separated (such that perfect classification performance is possible), then you should use a different classifier such as a Support Vector Machine (hard margin) or a Decision Tree.

Like all other classifiers, the MoRPE classifier should not be applied blindly to a given data set.  Before MoRPE is applied, the analyst must design a feature space that minimizes category fragmentation (as discussed in the paper).  This means that the analyst must spend time visualizing the data in many possible feature spaces, and carefully select the set of features that appears to minimize fragmentation.  The analyst doesn't need to be perfect, but they should make a decent effort and then MoRPE will handle the rest (as it is designed to do).

MoRPE is intended for scenarios where you have at least a few hundred samples per category, and where the number of categories is relatively small (between 2 and 10).  MoRPE works best for 2-category problems (where it has the fewest free parameters for a constant polynomial rank).  If MoRPE has H parameters for the 2-category problem, it has M*H parameters for the M-category problem where M > 2 (for a constant polynomial rank).

MoRPE approximates the optimal classifier when category fragmentation is low, and when irrelevant dimensions are removed from the feature space (as discussed in the paper).

#Software Releases
I am planning to provide this software in C# so others can use it.  Keep in mind... I can only work on this project in my spare time.

##C# 
This version is halfway complete.  It basically lacks an optimization routine, but is otherwise complete.

I started a C# codebase many years after finishing the Matlab codebase.  I was halfway finished before I got hired to complete the work privately for my employer.  My employment contract prevents me from publishing that work in the public domain.  So for now, the C# version is halfway complete (without an optimization routine).  In the future, I may rewrite the optimization routine in my spare time, and post a complete version.

##Java
I am currently working on a private Java version for my  employer.  This version works on a GPU-enabled cluster.  I cannot release that code in the public domain.  Contact me if you want to get your hands on that version.  (It would need to go through my employer, and I don't know what kind of terms they would negotiate.)

##Matlab
This is the original implementation of MoRPE.  **This code is filthy.**  I invented MoRPE almost by accident because nothing else seemed to fit the data I collected for my dissertation.  I tried many different things before MoRPE emerged from the chaos, and so this code is very chaotic.  There is plenty of dead code here, and the C code is almost impossible to understand (this was my first time writing anything in C or C++).

In order for the C files to be callable from Matlab, each file needs to be "mexed" from the Matlab command line.  Here I have already mexed the C files for 32-bit and 64-bit Windows, and so the software will work from a 32-bit or 64-bit Windows machine.  If you want this software to work on a different platform, you might need to mex the files from that platform (but it might work as is... I'm not sure).  For more information, see the Matlab documentation for `mex`.

Matlab Script for mexing files (if you need to):

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

#Legal stuff

##License
The legal license file is attached (its Apache 2.0).  The owner, licensor, and copyright holder is Almon David Ing, PhD.

##Not patented
I have allowed the patent application to lapse.  This work is no longer patented.
