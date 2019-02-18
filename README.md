# BLQH
Exact diagonalization of bilayer quantum Hall (BLQH) systems in the torus geometry

## Description
This program performs exact diagonalization and calculate the excitation spectrum of the bilayer/monolayer quantum Hall
system (BLQH/MLQH) in the torus geometry and in the lowest-Landau-level (LLL). The bilayer quantum Hall system contains two layers of two-dimensional electron gas (2DEG), with intralayer and interlayer Coulomb interaction. Including the spin degrees of freedom, the noninteracting hamiltonian with bias voltage, Zeeman splitting energy and interlayer tunneling (symmetric-antisymmetric energy gap) can be included. Exact diagonalizing the full hamiltonian matrix, we could get all the eigenstates and analyze the results. To shrink the Hilbert space, and simplify the calculation, we utilize the magnetic translation operators in x and y directions, and then apply a three-pass Lanczos algorithm to calculate the lowest eigenstate at each K-point.


## Prerequisites
Environment:linux

Libraries:

* icc,mkl: https://software.intel.com/en-us/parallel-studio-xe

* g++,liblapack-dev,libblas-dev

* GNU Scientific Library(GSL): https://www.gnu.org/software/gsl/

## Installing
ubuntu:

>sudo apt-get install g++ g++-multilib build-essential liblapack-dev libblas-dev libgsl-dev
>
>make blqh

with icc,mkl installed:
>make blqh_mkl
## Usage
<pre><code>./blqh -h  
./blqh: option requires an argument -- 'h'   
Usage: ./blqh [Options]  
Options:  
  -l                       nLL  
  -n                       nphi  
  -e                       Total No. of electrons in upper-layer  
  -u                       No. of electrons in upper-layer  
  -k                       kx in upper-layer  
  -S                       Delta_SAS: tunnelling amplitude  
  -v                       Delta_V: bias voltage  
  -j                       total J in upper-layer or down-layer  
  -g                       gamma=lx/ly  aspect ratio  
  -d                       interlayer distance  
  -m                       Lambda  
  -t                       nthread  
  -s                       seed  
Default: (l,n,u,d,lambda) = (4,4,2,1,200)  
</code></pre>

## Running the tests
We calculate 8-electrons with 4-electrons in the upper-layer and 4-electrons in  
the down layer at layer distance d/l=0.0. The ground state energy is at K-point (j,k)=(4,4), we use four threads to calculate the hamiltonian.
> ./blqh -n 8 -e 8 -u 4 -j 4 -k 4 -t 4 -d 0.0
<pre><code>----------- ED results --------------
nHilbert: =79
E_gs:= -0.636707
----------- Lanczos results ---------
E_gs:= -0.6367069257
# ground state wave function
  71 :   |  _  1  2  3  4  _  _  _)|  0  _  _  _  _  5  6  7)   000000001110000100011110  0.3380617019
   6 :   |  0  _  _  _  4  _  6  7)|  _  1  2  3  _  5  _  _)   000000000010111011010001  0.3380617019
  48 :   |  0  _  _  _  4  5  6  _)|  _  1  2  3  _  _  _  7)   000000001000111001110001  0.3380617019
  29 :   |  0  _  2  _  _  5  _  7)|  _  1  _  3  4  _  6  _)   000000000101101010100101  0.3380617019
  20 :   |  0  _  _  _  4  5  _  7)|  _  1  2  3  _  _  6  _)   000000000100111010110001  0.3380617019
  10 :   |  0  _  _  3  _  _  6  7)|  _  1  2  _  4  5  _  _)   000000000011011011001001  0.3380617019
  54 :   |  0  _  2  _  _  5  6  _)|  _  1  _  3  4  _  _  7)   000000001001101001100101  0.3380617019
  25 :   |  0  _  _  3  _  5  _  7)|  _  1  2  _  4  _  6  _)   000000000101011010101001  0.3380617019
  66 :   |  0  1  _  _  4  5  _  _)|  _  _  2  3  _  _  6  7)   000000001100110000110011  0.2390457219
  57 :   |  0  _  2  _  4  _  6  _)|  _  1  _  3  _  5  _  7)   000000001010101001010101  0.1690308509
  19 :   |  0  _  _  _  4  _  6  7)|  0  _  2  3  _  _  6  _)   000000000100110111010001  1.306913788e-13
</code></pre>

The ground state energies of ED and Lanczos approaches are the same within accuracy of 1E-5.   
The first column of the wave function (WF) is the ID of the basis, the |\*\*)|\*\*) columns are
electrons which occupy the orbitals ('_' stands for not occupied,left for upper-layer and right for down-layer). the next column is
the electron occupaptions in binary representation ('1' for occupied and '0' for empty), and the last column are the coefficients in the many-body basis.  
As we can see, at d/l=0.0, all the basis with nonzero coefficients are Haplerin "111 state".

At distance d/l=3.0, we could see that the all the basis with nonzero coefficients are two decoupled composite Fermi liquids.
<pre><code>----------- ED results --------------
nHilbert: =79
E_gs:= -1.23588
----------- Lanczos results ---------
E_gs:= -1.235882219
# ground state wave function
  66 :   |  0  1  _  _  4  5  _  _)|  _  _  2  3  _  _  6  7)   000000001100110000110011  0.5388040068
  52 :   |  0  _  _  3  4  _  _  7)|  0  _  _  3  4  _  _  7)   000000001001100110011001  0.530402778
  51 :   |  _  1  _  3  4  _  6  _)|  0  _  _  3  4  _  _  7)   000000001001100101011010  0.3033881081
  33 :   |  0  _  2  _  _  5  _  7)|  _  1  2  _  _  5  6  _)   000000000110011010100101  0.3033881081
   5 :   |  _  _  2  3  _  _  6  7)|  0  _  2  3  _  5  _  _)   000000000010110111001100  0.3033881081
  68 :   |  0  1  _  _  4  5  _  _)|  _  1  _  _  4  _  6  7)   000000001101001000110011  0.3033881081
  29 :   |  0  _  2  _  _  5  _  7)|  _  1  _  3  4  _  6  _)   000000000101101010100101  0.1269396409
  55 :   |  0  _  _  3  _  5  6  _)|  0  _  2  _  _  5  _  7)   000000001010010101101001  0.1215007838
  35 :   |  0  _  2  _  _  5  _  7)|  0  _  _  3  _  5  6  _)   000000000110100110100101  0.1215007838
  34 :   |  0  _  _  3  _  5  6  _)|  0  _  _  3  _  5  6  _)   000000000110100101101001  0.1173494879
   3 :   |  0  1  _  _  4  5  _  _)|  _  1  2  3  4  _  _  _)   000000000001111000110011  0.009704390136
</code></pre>
