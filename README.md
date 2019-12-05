# d_RPCA
Deformable Groupwise Image Registration using Low-Rank and Sparse Decomposition

* A new groupwise distance measure for non-parametric image registration based on Low-Rank- and Sparse-Decomposition, inspired by [Heber and Pock, 2014](https://link.springer.com/chapter/10.1007/978-3-319-10599-4_48)
* Regularization through Total Variation penalty
* Optimization is performed using iterative relinearization and the primal-dual hybrid gradient algorithm from [Chambolle and Pock, 2011](https://hal.archives-ouvertes.fr/hal-00490826/document)

#### Implementation in __MATLAB__. 

Simply run

> Demo/demo_mf_nn_registration_no_ref_ml.m

for a demo on synthetic image data.

As a method of comparison, a simple variance-based groupwise registration method as in [Metz et al., 2011](https://www.sciencedirect.com/science/article/abs/pii/S1361841510001155?via%3Dihub) is implemented in 

 > Demo/demo_var_registration_no_ref_ml.m

In order to run the demo script for optical flow estimation
 
 > Demo/demo_mf_nn_registration_fix_ref_ml.m

you need to request the benchmark data from [HCI](https://hci.iwr.uni-heidelberg.de/benchmarks/Challenging_Data_for_Stereo_and_Optical_Flow) and unpack it into

> Data/ChallengingSequences

###### Contact for questions: haase [at] mic.uni-luebeck.de
