# BSc_thesis
Optimal Power Flow problems aim to minimize the power generating cost
while satisfying network constraints. Here we develop our own Sequential Linear
Programming and Sequential Quadratic Programming solvers from rst principles
in MATLAB. We use a MATLAB-AMPL interface to obtain the AC OPF model
data from AMPL. These solvers will aim to converge to optimal solutions for different
sized Optimal Power Flow problems. We nd SQP to be more reliable than
SLP. However SLP converges faster for cases where it obtains a solution. We conclude
that the solutions provided by the SLP solver are not optimal and our SQP
solver only nds optimal solutions for two cases. We suggest further alterations
to the algorithms we have written to improve optimality of the solutions.
