# An active-set method for Convex QP with piecewise linear terms

This is a MATLAB implementation of an active-set method suitable for convex quadratic programming problems with piecewise linear terms of the following form:

$$ \min_{x \in \mathbb{R}^n}\  c^\top x + \frac{1}{2} x^\top Q x + \sum_{i=1}^l \max((Cx+d)_i,0) + ||Dx||_1,\qquad \text{s.t., }Ax = b,\ x \in [a_l,a_u],$$

where $Q \in \mathbb{R}^{n\times n}$ is a positive semidefinite matrix, $c \in \mathbb{R}^n$, $D \in \mathbb{R}^{n\times n}$ is a diagonal "weight" matrix,
$C \in \mathbb{R}^{l\times n}$ and $d \in \mathbb{R}^l$ form the piecewise-linear (max) terms in the objective, and $A \in \mathbb{R}^{m\times n}$, $b \in \mathbb{R}^m$
are the data associated with the linear constraints of the problem. The decision vector $x$ is restricted to the box $[a_l,a_u]$, with $a_l\leq a_u$ (noting that 
$a_l$ can have $-\infty$ entries, while $a_u$ can have $+\infty$ ones). 

<br/>

The code is based on an accompanying research paper, in which we derive an appropriate Proximal Method of Multipliers (PMM) combined with
a SemiSmooth Newton method (SSN), the associated linear systems of which are solved via a preconditioned Krylob subspace solver. This repo is 
dedicated towards the reproducibility of the numerical results presented in the accompanying paper; however, the associated code is well-commented 
and is intented to be used as a template for research purposes.


<br/>
The core file containing the basic active-set method is <strong>SSN_PMM.m</strong>. <br/>

---
<strong> Input: </strong> A MATLAB struct <code> pb_struct </code>, with the following fields:

<ul>
  <li> <code>.Q</code> -> the (sparse) coefficient matrix of the quadratic in the objective (if empty, provide as <code>sparse(n,n)</code>) </li>
  <li> <code>.A</code> -> the (sparse) linear equalities coefficient matrix (if empty, provide as <code>sparse(0,n)</code>) </li>
  <li> <code>.C</code> -> the (sparse) coefficient matrix appearing within the $max(\cdot,0)$ term in the objective (if empty, provide as <code>sparse(0,n)</code>) </li>
  <li> <code>.C</code> -> the (sparse) coefficient matrix appearing within the $max(\cdot,0)$ term in the objective (if empty, provide as <code>sparse(0,n)</code>) </li>
  <li> <code>.c</code> -> the linear coefficients of the objective function (<em>Default</em>: all zeros) </li>
  <li> <code>.b</code> -> the right hand side of the linear equalities (Default: all zeros) (<em>Default</em>: all zeros) </li>
  <li> <code>.d</code> -> the constant displacement in $max(\cdot,0)$ terms in the objective (<em>Default</em>: all zeros) </li>
  <li> <code>.lb</code> -> the lower bound vector on the primal variables $x$ (<em>Default</em>: all $-\infty$) </li>
  <li> <code>.ub</code> -> the upper bound vector on the primal variables $x$ (<em>Default</em>: all $+\infty$) </li>
  <li> <code>.D</code> -> the weight vector for possible ell-1 norm in the objective (<em>Default</em>: all zeros) </li>
  <li> <code>.tol</code> -> the specified tolerance for termination (<em>Default</em>: $10^{-4}$) </li>
  <li> <code>.maxit</code> -> the maximum allowed number of PMM iterations (<em>Default</em>: 200)  </li>
  <li> <code>.pl</code> -> Possible choices <ul>
      <li>0 for not printing intermediate iterates</li>
      <li>1 for printing only PMM iterates (<em>Default</em>)</li>
      <li>2 for additionally printing SNN iterates</li>
      <li>3 for additionally printing Krylov iterate info</li>
    </ul> </li>
  <li> <code>.p_fid</code> -> file ID to print output (<em>Default</em>: 1 (prints on workspace)).  </li>

</ul>

---

---
<strong> Output: </strong> A MATLAB struct <code> pb_struct </code>, with the following fields:

<ul>
  <li> <code>.Q</code> -> the (sparse) coefficient matrix of the quadratic in the objective (if empty, provide as <code>sparse(n,n)</code>) </li>
  <li> <code>.A</code> -> the (sparse) linear equalities coefficient matrix (if empty, provide as <code>sparse(0,n)</code>) </li>
  <li> <code>.C</code> -> the (sparse) coefficient matrix appearing within the $max(\cdot,0)$ term in the objective (if empty, provide as <code>sparse(0,n)</code>) </li>
  <li> <code>.C</code> -> the (sparse) coefficient matrix appearing within the $max(\cdot,0)$ term in the objective (if empty, provide as <code>sparse(0,n)</code>) </li>
  <li> <code>.c</code> -> the linear coefficients of the objective function (<em>Default</em>: all zeros) </li>
  <li> <code>.b</code> -> the right hand side of the linear equalities (Default: all zeros) (<em>Default</em>: all zeros) </li>
  <li> <code>.d</code> -> the constant displacement in $max(\cdot,0)$ terms in the objective (<em>Default</em>: all zeros) </li>
  <li> <code>.lb</code> -> the lower bound vector on the primal variables $x$ (<em>Default</em>: all $-\infty$) </li>
  <li> <code>.ub</code> -> the upper bound vector on the primal variables $x$ (<em>Default</em>: all $+\infty$) </li>
  <li> <code>.D</code> -> the weight vector for possible ell-1 norm in the objective (<em>Default</em>: all zeros) </li>
  <li> <code>.tol</code> -> the specified tolerance for termination (<em>Default</em>: $10^{-4}$) </li>
  <li> <code>.maxit</code> -> the maximum allowed number of PMM iterations (<em>Default</em>: 200)  </li>
  <li> <code>.pl</code> -> Possible choices <ul>
      <li>0 for not printing intermediate iterates</li>
      <li>1 for printing only PMM iterates (<em>Default</em>)</li>
      <li>2 for additionally printing SNN iterates</li>
      <li>3 for additionally printing Krylov iterate info</li>
    </ul> </li>
  <li> <code>.p_fid</code> -> file ID to print output (<em>Default</em>: 1 (prints on workspace)).  </li>

</ul>

---


<br/><br/><br/><br/>

---
> **_NOTE:_** The datasets utilized by the code are not provided (due to licensing reasons). Nonetheless, we have included readme.txt files in
the appropriate places where the user should include the relevant datasets. Within these readme files, one can find the relevant links for
downloading the associated datasets. 
---



