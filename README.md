#MoRPE:  A Probabilistic Classifier
The classification problem is commonly encountered when a finite sample of data is leveraged to determine the probabilistic relationship between a category label  *c*  and a multivariate coordinate  **x**   for an entire population.  Solving this problem requires approximating the optimal classifier, a function of  **x**   that evaluates the conditional probability  *p*(*c*¦**x**)  in an optimal way.  The paper *Ing_Morpe.docx* introduces MoRPE, a method for approximating optimal classifiers.  MoRPE is a machine learning method for probabilistic classification based on **Mo**notonic **R**egression of a **P**olynomial **E**xpansion.  It is conceptually related to Fisher’s Quadratic Discriminant (Fisher 1936) and Kernel Machines.  MoRPE has the ability to approximate an optimal classifier with remarkable precision in common scenarios.

#Academic Citation
### [Release 0.1](https://github.com/adaviding/Morpe/releases/tag/0.1) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13235.svg)](http://dx.doi.org/10.5281/zenodo.13235)

Ing, AD (2014) MoRPE:  A Probabilistic Classifier.  *GitHub*.  doi:10.5281/zenodo.13235.  Retrieved from https://github.com/adaviding/Morpe/releases/tag/0.1.

#Software Releases
I am planning to provide this software in multiple programming languages.  Keep in mind... I can only work on this project in my spare time, and I don't have much spare time.

##C++
Not ready.
	
##C# 
Not ready.  Expected to be ready on July 1, 2015.

##Java
Not ready.

##Matlab
This is the original implementation of MoRPE.  It was written in Matlab and C so that I could conduct analyses for my PhD Dissertation.  This code was not written for mass consumption.  It contains some dead code, and the C code is particularly ugly.  I will eventually clean this up AFTER I have cross-coded MoRPE to C++.

In order for the C files to be callable from Matlab, each file needs to be "mexed" from the Matlab command line.  Here I have already mexed the C files for 32-bit Windows, and so the software will work from a 32-bit Windows machine.  If you want this software to work on a different platform, you must mex the necessary C files.  For more information, see the Matlab documentation for `mex`.

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
