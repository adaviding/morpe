# MoRPE:  A Probabilistic Classifier
The classification problem is commonly encountered when a finite sample of data is leveraged to determine the probabilistic relationship between a category label  *c*  and a multivariate coordinate  **x**   for an entire population.  Solving this problem requires approximating the optimal classifier, a function of  **x**   that evaluates the conditional probability  *p*(*c*¦**x**)  in an optimal way.  The paper *Ing_Morpe.docx* introduces MoRPE, a method for approximating optimal classifiers.  MoRPE is a machine learning method for probabilistic classification based on **Mo**notonic **R**egression of a **P**olynomial **E**xpansion.  It is conceptually related to Fisher’s Quadratic Discriminant (Fisher 1936) and Kernel Machines.  MoRPE has the ability to approximate an optimal classifier with remarkable precision in common scenarios.

## Academic Citation
### [Release 0.1](https://github.com/adaviding/Morpe/releases/tag/0.1) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13235.svg)](http://dx.doi.org/10.5281/zenodo.13235)

Ing, AD (2014) MoRPE:  A Probabilistic Classifier.  *GitHub*.  doi:10.5281/zenodo.13235.  Retrieved from https://github.com/adaviding/Morpe/releases/tag/0.1.

## Project Status
This project moves very slowly because I have very little time to work on it.

### Projects that will be supported moving forward
Here are a couple implementations which are roughly equivalent in terms of the outer API and inner numerical properties, but computational performance may differ. 

* The [C# codebase](src/csharp/readme.md) is usable.
* The [C++ codebase](src/cpp/readme.md) is just getting started.

### The original project is unsupported
This was the code I had when MoRPE was first invented.

* The [Matlab codebase](src/matlab/readme.md) is here for posterity, but it is filthy and you should try to avoid it.

## Recommended Usage
MoRPE should be used in scenarios when the categories are distributed in a real-numerical feature space with partial overlap (such that perfect classification performance is impossible).

Like all other classifiers, the MoRPE classifier should not be applied blindly to a given data set.  Before MoRPE is applied, the analyst must design a feature space that minimizes category fragmentation (as discussed in the paper).  This means that the analyst must spend time visualizing the data in many possible feature spaces, and carefully select the set of features that appears to minimize fragmentation.  The analyst doesn't need to be perfect, but they should make a decent effort and then MoRPE will handle the rest (as it is designed to do).

MoRPE is intended for scenarios where you have at least a few hundred samples per category, and where the number of categories is relatively small (between 2 and 10).  MoRPE works best for 2-category problems (where it has the fewest free parameters for a constant polynomial rank).  If MoRPE has H parameters for the 2-category problem, it has M*H parameters for the M-category problem where M > 2 (for a constant polynomial rank).

MoRPE approximates the optimal classifier when category fragmentation is low, and when irrelevant dimensions are removed from the feature space (as discussed in the paper).

## License
The [MIT License](license.md).
