# An active-set method for Convex QP with piecewise linear terms

This is a MATLAB implementation of an active-set method suitable for convex quadratic programming problems with piecewise linear terms of the following form:

$$ \min_{x \in \mathbb{R}^n}\  c^\top x + \frac{1}{2} x^\top Q x + \sum_{i=1}^l \max((Cx+d)_i,0) + ||Dx||_1,\qquad \text{s.t., }Ax = b,\ x \in [a_l,a_u],$$

where $Q \in \mathbb{R}^{n\times n}$ is a positive semidefinite matrix, $D \in \mathbb{R}^{n\times n}$ is a diagonal "weight" matrix,
$C \in \mathbb{R}^{l\times n}$, $d \in \mathbb{R}^l$ (form the piecewise-linear terms in the objective) and $A \in \mathbb{R}^{m\times n}$, $b \in \mathbb{R}^m$
are the data associated with the linear constraints of the problem. 

The code is based on an accompanying research paper, in which we derive an appropriate Proximal Method of Multipliers (PMM) combined with
a SemiSmooth Newton method (SSN), the associated linear systems of which are solved via a preconditioned Krylob subspace solver. This repo is 
dedicated towards the reproducibility of the numerical results presented in the accompanying paper; however, the associated code is well-commented 
and is intented to be used as a template for research purposes.



**NOTE:** The datasets utilized by the code are not provided (due to licensing purposes). Nonetheless, we have included readme.txt files in
the appropriate places where the user should include the relevant datasets. Within these readme files, one can find the relevant links for
downloading the associated datasets. 




