# d_RPCA
Deformable Groupwise Image Registration using Low-Rank and Sparse Decomposition

* A new groupwise distance measure for non-parametric image registration based on low-rank and sparse decompositions, inspired by [Heber and Pock, 2014](https://link.springer.com/chapter/10.1007/978-3-319-10599-4_48)
* Regularization through Total Variation or Curvature Penalty
* Optimization is performed using iterative relinearization and the primal-dual hybrid gradient algorithm from [Chambolle and Pock, 2011](https://hal.archives-ouvertes.fr/hal-00490826/document)

#### Implementation in __MATLAB__. 

Add your local copy of this repository to your MATLAB path. Then, simply run

> Demo/demo_mf_nn_tv_registration_no_ref_ml.m

for a demo on synthetic image data with TV regularization.

As a method of comparison, a simple variance-based groupwise registration method as in [Metz et al., 2011](https://www.sciencedirect.com/science/article/abs/pii/S1361841510001155?via%3Dihub) is implemented as well. A demo is found in

 > Demo/demo_var_tv_registration_no_ref_ml.m
