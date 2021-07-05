# tensor-ring-decomposition

Code accompanying the paper "On algorithms for and computing with the tensor ring decomposition", Numerical Linear Algebra with Applications 27.3 (2020): e2289. If you find the code useful, please consider citing the paper.

Arxiv version available at: https://arxiv.org/abs/1807.02513

README:
=======

Files marked 'ex_...' contain the numerical experiments in the paper.
The file 'ex_tr_wave_implicit_cp.m' requires Oseledet et al's TT-toolbox
in the current path (https://github.com/oseledets/TT-Toolbox).

The remaining files are standalone with no installation required.


Main files:
-----------
normTR.m
 - computes norm of tensor in TR-format.

 roundingTR.m
 - rounds a TR-representation with a bound on the Frobenius norm.

TRdecomp.m
  - computes TR-representation of a tensor.

TRcheckall.m
  - computes TR-representation of a tensor with low storage cost after
    an exhaustive search.

TRprealt2.m
  - computes TR-representation of a tensor with low storage cost using
    a heuristic approach.

CP2TR.m
  - comverts canonical decomposition into TR-format.

CP2TRfaster.m
  - comverts canonical decomposition into TR-format. Computations sped up by rounding subsystems.

CP2TRcheckall.m
  - comverts canonical decomposition into TR-format with low storage cost after
    an exhaustive search

CP2TRcheckallfaster.m
  - comverts canonical decomposition into TR-format with low storage cost after
    an exhaustive search. Computations sped up by rounding subsystems.

insert_gen_edge.m
  - inserts an edge between indices in a TR-representation.

TR2TT.m
  - converts TR-representation into TT-format.

TT2TR.m
  - converts TT-representation into TR-format.

fullTR.m
  - converts tensor in TR-format into full format.

check_all_edge.m
  - inserts exactly one edge between the first index of a TR-representation and a different index. The index resulting in the lowest storage cost is chosen.
