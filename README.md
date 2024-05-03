# An active-set method for Convex QP with piecewise linear terms

This is a MATLAB implementation of an active-set method suitable for convex quadratic programming problems with piecewise linear terms of the following form:

$$ \min_{x \in \mathbb{R}^n}\  c^\top x + \frac{1}{2} x^\top Q x + \sum_{i=1}^l \max((Cx+d)_i,0) + \|Dx\|_1,\qquad \text{s.t., }Ax = b,\ x \in [a_l,a_u],$$

where $Q \in \mathbb{R}^{n\times n}$ is a positive semidefinite matrix, $D \in \matbbb{R}^{n\times n}$ is a diagonal "weight" matrix,
and $C \in \mathbb{R}^{l\times n}$, $d \in \mathbb{R}^l$. 

The code is based on an accompanying research paper, in which we derive an appropriate Proximal Method of Multipliers (PMM) combined with
a SemiSmooth Newton method (SSN), the associated linear systems of which are solved via a preconditioned Krylob subspace solver. This repo is 
dedicated towards the reproducibility of the numerical results presented in the accompanying paper; however, the associated code is well-commented 
and is intented to be used as a template for research purposes.




