# entr_ineq
MATLAB package for calculating entropic inequalities.  
[PACKAGE IS NOT COMPLETE. WE'RE CURRENTLY UPLOADING IT.]

The main goal of the package is to calculate entropic inequalities from
basic Shannon or von Neumann inequalities for a given set of random
variables and possible independence constraints among these variables.
These constraints can be drawn from graph of a causal model (Bayesian
Network or Markov Random Field).   

The main modules of the package are:  
-- Generation of basic inequalitites for a set of variables (e.g. Shannon,
von Neuman),  
-- d-separation algorithm for directed acyclic graphs of Bayesian Networks
and a similar algorithm for graphs of Markov Random Fields,   
-- Calculation of intersection of a cone given by set of inequalitites
and a cone given by a set of equations,  
-- Fourier-Motzkin elimination algorithm for projections of polytopic
structures (cones or polytopes).  

This packege contains additional programs which can be usefull for the
computations (e.g. generation of data processing inequalities, Graham
algorithm for determinig cyclicity of hypergraphs).  

SEE FULL DESCRIPTION IN THE MANUAL FILE (WORK IN PROGRESS).
