# Multi-Configuration many particle systems : Exemplification of the numerical implementation

## The problem studied

Here it is presented an numerical implementation based on C code for the
indexing problem of multi-configurational(truncated *Fock* states basis)
method for many-particle systems. A Perfect hashing function is already a
well known for both fermions and bosons, as cited in the main article,
though the literature lacks on specific details about the time required to
perform some fundamental calculations of relevant quantities, such as those
explored here, the one and two-body density matrices and the action of
Hamiltonian operator in a many-particle state represented in the Fock
basis. Additionally, it is explored some mapping structures to describe the
action of creation and annihilation operators in the Fock states, directly
linking the Fock states using a memory access instead of function calls.
This boost the performance though demands more memory. Finally, a parallel
implementation is studied both on CPU and GPU, in order to reveal the
scalability of the numerical problem.

## Hierarchy of files

To read the codes and understand the module, it is suggested to read the
file 'article.pdf' to get started on the problem. For those experienced
with the multi-configuration method the reading may be skipped. Among the
header files, the 'configurationsMap.h'
implement the basic functions, such as data types used, the hashing table
setup and the mappings used in improvements. All other headers depends on
this one. The name of the headers are intended to be very suggestive to
indicate what physical quantity is being calculated and a resonable
description is given on the files. To illustrate the functionality of the
basic functions on '*configurationMap.h*' one can compile and run the
'*demonstrateFockMap.c*' file, where the specific instructions are given in
the file. The main program used to collect the results given in the article
is 'performanceTest.c', which details the time required by different
implementations using the improving mapping structures. Furthermore, the
cuda folder implement specifically the most improved hamiltonian(all mappings
of particle creation/annihilation) acting on the Fock basis to explore the
scalability of the problem on GPUs.
