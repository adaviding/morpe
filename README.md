#MoRPE:  A Probabilistic Classifier
The classification problem is commonly encountered when a finite sample of data is leveraged to determine the probabilistic relationship between a category label  *c*  and a multivariate coordinate  **x**   for an entire population.  Solving this problem requires approximating the optimal classifier, a function of  **x**   that evaluates the conditional probability  *p*(*c*¦**x**)  in an optimal way.  The paper *Ing_Morpe.docx* introduces MoRPE, a method for approximating optimal classifiers.  MoRPE is a machine learning method for probabilistic classification based on **Mo**notonic **R**egression of a **P**olynomial **E**xpansion.  It is conceptually related to Fisher’s Quadratic Discriminant (Fisher 1936) and Kernel Machines.  MoRPE has the ability to approximate an optimal classifier with remarkable precision in common scenarios.

#Academic Citation
### [Release 0.1](https://github.com/adaviding/Morpe/releases/tag/0.1) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13235.svg)](http://dx.doi.org/10.5281/zenodo.13235)

Ing, AD (2014) MoRPE:  A Probabilistic Classifier.  *GitHub*.  doi:10.5281/zenodo.13235.  Retrieved from https://github.com/adaviding/Morpe/releases/tag/0.1.

#Project Status
The Matlab codebase is functioning, but it is not user friendly.  (See the Matlab notes below.)

I am slowly getting together a C# code base (in my spare time).  I am hoping to finish by the end of July.  (See the C# notes below.)  Once C# is done, I will create a Java version.

#Recommended Usage
MoRPE is one of the best classifiers available in the public domain (perhaps *the* best).  However, this does not mean that it can be applied blindly to any data set.  Before MoRPE is applied, the analyst must design a feature space that minimizes category fragmentation (as discussed in the paper).  This means that the analyst must spend time visualizing the data in many possible feature spaces, and carefully select the set of features that appears to minimize fragmentation.  The analyst doesn't need to be perfect, but they should make a decent effort and then MoRPE will handle the rest (as it is designed to do).

MoRPE is intended for scenarios where you have at least a few hundred samples per category, and where the number of categories is relatively small (between 2 and 10).  MoRPE works best for 2-category problems (where it has the fewest free parameters for a constant polynomial rank).  If MoRPE has H parameters for the 2-category problem, it has M*H parameters for the M-category problem where M > 2 (for a constant polynomial rank).

MoRPE approximates the optimal classifier when category fragmentation is low, and when irrelevant dimensions are removed from the feature space (as discussed in the paper).

#Software Releases
I am planning to provide this software in multiple programming languages.  Keep in mind... I can only work on this project in my spare time, and I don't have much spare time.  It's funny how "work" prevents us from creating stuff that is economically valuable.

##C++
Not ready.  (Not close.)
	
##C# 
Not ready.  Expected to be ready on August 1, 2015.

##Java
Not ready.

##Matlab
This is the original implementation of MoRPE.  **This code is filthy.**  I invented MoRPE almost by accident because nothing else seemed to fit the data I collected for my dissertation.  I tried many different things before MoRPE emerged from the chaos, and so this code is very chaotic.  There is plenty of dead code here, and the C code is almost impossible to understand (this was my first time writing anything in C or C++).

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
