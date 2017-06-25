# entr_ineq
MATLAB package for calculating entropic inequalities.
[PACKAGE IS NOT COMPLETE. I'M CURRENTLY UPLOADING IT.]

The main goal of the package is to calculate entropic inequalities for a
given set of random variables and possible independence constraints (IC) 
among these variables. IC can be encoded in Markov Random Field (MRF),
Directed Acyclic Graph (DAG) or given explicitly. 

The main modules of the package are:
-- Generation of basic inequalitites for a set of variables (e.g. Shannon,
von Neuman),
-- d-separation algorithm for DAGs and similar for MRFs giving sets of
independence constraints associated with causal structure,
-- Calculation of intersection of a cone given by set of inequalitites
and a cone given by a set of equations,
-- Fourier-Motzkin elimination algorithm for projections of polytopic
structures (cones or polytopes).

This packege contains additional programs which can be usefull for the
computations (e.g. generation of data processing inequalities, Graham
algorithm for determinig cyclicity of hypergraphs).

SEE FULL DESCRIPTION IN THE MANUAL FILE (WORK IN PROGRESS).
