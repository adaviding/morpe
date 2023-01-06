# c#

## status
This version is usable.

## summary
This is the modern implementation of MoRPE.  My plan is to keep the C# version in sync with other modern versions.

## roadmap
* More testing
  * We currently have minimal testing which just proves basic overall functionality of the training algorithm.
  * We need unit tests for all the encapsulated "units" of functionality.
  * We need many more scenarios for classifier training.
* Visualization
  * The goal is to generate "subspace" plots, similar to figures 7, 8, 9, and 10 from the paper [Ing_Morpe.docx](../../Ing_Morpe.docx).
  * The plots would be SVG, but would incorporate a bitmapped graphics file for each 2D cloud.
* Multiple classifier management
  * We would expand the `TrainingContext` class to be a container for many classifier variants.
    * All 1D classifiers (1 per dimension)
    * All 2D classifiers (1 per dimension pair)
    * All full-dimensionality classifiers.
    * Classifiers of various polynomial rank.
  * This supports two kinds of operations.
    * The cross-validation procedure, used to choose a best-fitting classifier.
    * The generation of model output for the "subspace" plots.
* Model selection via cross-validation
  * Code to train multiple classifier variants, to pick the best variant.
