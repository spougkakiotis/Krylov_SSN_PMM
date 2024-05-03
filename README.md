# An active-set method for Convex QP with piecewise linear terms

This is a MATLAB implementation of an active-set method suitable for convex quadratic programming problems with piecewise linear terms of the following form:

$$ \min_{x \in \mathbb{R}^n}\  c^\top x + \frac{1}{2} x^\top Q x + \sum_{i=1}^l \max\{(Cx+d)_i,0\},\qquad \text{s.t. }Ax = b,\ x \in [a_l,a_u],$$

where $Q \in \mathbb{R}^{n\times n}$ is a positive semidefinite matrix, and $(\cdot)_+ \equiv \max\{\cdot,0\}$.
