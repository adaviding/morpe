#MoRPE:  A Probabilistic Classifier
The classification problem is commonly encountered when a finite sample of data is leveraged to determine the probabilistic relationship between a category label  *c*  and a multivariate coordinate  **x**   for an entire population.  Solving this problem requires approximating the optimal classifier, a function of  **x**   that evaluates the conditional probability  *p*(*c*¦**x**)  in an optimal way.  This paper introduces MoRPE, a new method for approximating optimal classifiers.  MoRPE is a machine learning method for probabilistic classification based on **Mo**notonic **R**egression of a **P**olynomial **E**xpansion.  It is conceptually related to Fisher’s Quadratic Discriminant (Fisher 1936).  MoRPE has the ability to approximate an optimal classifier with remarkable precision in common scenarios.

#Academic Citation
### [Release 0.1](https://github.com/adaviding/Morpe/releases/tag/0.1) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13235.svg)](http://dx.doi.org/10.5281/zenodo.13235)

Ing, AD (2014) MoRPE:  A Probabilistic Classifier.  *GitHub*.  doi:10.5281/zenodo.13235.  Retrieved from https://github.com/adaviding/Morpe/releases/tag/0.1.

#Not patented
I have allowed the patent application to lapse.  This work is no longer patented.

#Software Releases
This software is implemented in multiple programming languages.

##Matlab
This is the original implementation of MoRPE.  I used it to conduct analyses for my PhD Dissertation.  I finished writing this code in 2009.

This release is actually a combination of Matlab and C.  In order for the C files to be callable from Matlab, each file needs to be "mexed" from the Matlab command line.  Here I have already mexed the C files for 32-bit Windows, and so the software will work from a 32-bit Windows machine.  If you want this software to work on a different platform, you must mex the necessary C files.  For more information, see the Matlab documentation for `mex`.

In this folder, you will find:
* A number of cpp-files.
	* Any cpp-file prefixed by "Mcl_" needs to be mexed (if not already mexed).
	* Any other cpp-file does not need to be mexed.  These files are used in `#include` directives.
* A number of m-files beginning with the `Mcl_` prefix.  These files are the implementation MoRPE (as MCL stands for MoRPE Classifier Library).
* A number of m-files beginning with the `s_` prefix.  These are script files that test parts of the MCL, or otherwise demonstrate its usage.
* Two examples of using MoRPE to classify data.
	* s_TestMcl_Example1.m
	* s_TestMcl_Example2.m
	
##Csharp
Not ready

##Java
Not ready

##C++
Not ready